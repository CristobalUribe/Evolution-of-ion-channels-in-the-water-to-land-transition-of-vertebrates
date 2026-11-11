#!/usr/bin/env python3
"""Build audited codeml branch-site inputs from GUIDANCE3 codon alignments.

The script creates independent alternative and null directories, converts each
validated FASTA alignment to PAML sequential format, copies the matching marked
tree, regenerates codeml.ctl from the previous alternative template, and checks
all model-defining settings.  It never copies previous likelihood outputs.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import sys
import tempfile
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Iterable


STOP_CODONS = {"TAA", "TAG", "TGA"}
IUPAC_OR_MISSING = set("ACGTRYSWKMBDHVN?-")
CONTROL_KEYS = (
    "seqfile",
    "treefile",
    "outfile",
    "model",
    "NSsites",
    "fix_omega",
    "omega",
    "cleandata",
)


@dataclass
class Result:
    og: str
    status: str
    error: str = ""
    n_sequences: int | str = ""
    nt_columns: int | str = ""
    codon_columns: int | str = ""
    complete_codon_columns: int | str = ""
    deleted_by_cleandata_columns: int | str = ""
    gap_codon_cells: int | str = ""
    nnn_codon_cells: int | str = ""
    other_ambiguous_codon_cells: int | str = ""
    tree_tips: int | str = ""
    foreground_marks: int | str = ""
    seqfile: str = ""
    alternative_dir: str = ""
    null_dir: str = ""


def natural_key(value: str) -> tuple[int, str]:
    match = re.fullmatch(r"OG(\d+)", value)
    return (int(match.group(1)), value) if match else (sys.maxsize, value)


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="\n",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        handle.write(text)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []

    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence).upper()))
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"empty FASTA header at {path}:{line_number}")
                sequence = []
            else:
                if header is None:
                    raise ValueError(
                        f"sequence before first FASTA header at {path}:{line_number}"
                    )
                sequence.append("".join(line.split()))

    if header is not None:
        records.append((header, "".join(sequence).upper()))
    if not records:
        raise ValueError(f"empty FASTA: {path}")
    if any(not sequence for _, sequence in records):
        raise ValueError(f"one or more empty sequences in {path}")
    return records


def validate_alignment(records: list[tuple[str, str]]) -> dict[str, int]:
    identifiers = [header.split()[0] for header, _ in records]
    if any(header != header.split()[0] for header, _ in records):
        raise ValueError("FASTA identifiers must not contain whitespace")
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("duplicate sequence identifiers")

    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise ValueError(f"unequal alignment lengths: {sorted(lengths)}")
    nt_length = lengths.pop()
    if nt_length % 3:
        raise ValueError(f"alignment length is not divisible by three: {nt_length}")

    partial_gap_codons = 0
    stop_codons = 0
    invalid_characters: set[str] = set()
    gap_cells = 0
    nnn_cells = 0
    other_ambiguous_cells = 0
    codon_count = nt_length // 3
    complete_columns = 0

    for _, sequence in records:
        invalid_characters.update(set(sequence) - IUPAC_OR_MISSING)
        for start in range(0, nt_length, 3):
            codon = sequence[start : start + 3]
            if "-" in codon and codon != "---":
                partial_gap_codons += 1
            if codon in STOP_CODONS:
                stop_codons += 1
            if codon == "---":
                gap_cells += 1
            elif codon == "NNN":
                nnn_cells += 1
            elif any(base not in "ACGT" for base in codon):
                other_ambiguous_cells += 1

    for column in range(codon_count):
        codons = [
            sequence[column * 3 : column * 3 + 3]
            for _, sequence in records
        ]
        if all(set(codon) <= set("ACGT") and codon not in STOP_CODONS for codon in codons):
            complete_columns += 1

    if invalid_characters:
        raise ValueError(
            "invalid alignment characters: " + "".join(sorted(invalid_characters))
        )
    if partial_gap_codons:
        raise ValueError(f"{partial_gap_codons} partial-gap codons")
    if stop_codons:
        raise ValueError(f"{stop_codons} unmasked stop codons")
    if complete_columns == 0:
        raise ValueError("cleandata=1 would leave zero codon columns")

    return {
        "n_sequences": len(records),
        "nt_columns": nt_length,
        "codon_columns": codon_count,
        "complete_codon_columns": complete_columns,
        "deleted_by_cleandata_columns": codon_count - complete_columns,
        "gap_codon_cells": gap_cells,
        "nnn_codon_cells": nnn_cells,
        "other_ambiguous_codon_cells": other_ambiguous_cells,
    }


def render_paml(records: Iterable[tuple[str, str]]) -> str:
    materialized = list(records)
    nt_length = len(materialized[0][1])
    lines = [f"{len(materialized)} {nt_length}"]
    for header, sequence in materialized:
        lines.append(header)
        lines.extend(sequence[start : start + 60] for start in range(0, len(sequence), 60))
    return "\n".join(lines) + "\n"


def read_tree(path: Path) -> tuple[str, set[str], int, int]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"empty tree file: {path}")

    declared_tips = -1
    declared_trees = -1
    if re.fullmatch(r"\d+\s+\d+", lines[0]):
        declared_tips, declared_trees = map(int, lines[0].split())
        newick = "".join(lines[1:])
    else:
        newick = "".join(lines)

    if not newick.endswith(";"):
        raise ValueError(f"Newick tree lacks final semicolon: {path}")
    tips = set(re.findall(r"(?:^|[(,])\s*([^():,;\s]+)", newick))
    if not tips:
        raise ValueError(f"no tree tips detected: {path}")
    foreground_marks = len(re.findall(r"(?:#|\$)1\b", newick))

    if declared_tips >= 0 and declared_tips != len(tips):
        raise ValueError(
            f"tree header declares {declared_tips} tips but {len(tips)} were parsed"
        )
    if declared_trees >= 0 and declared_trees != 1:
        raise ValueError(f"tree header declares {declared_trees} trees; expected 1")
    if foreground_marks != 1:
        raise ValueError(
            f"expected exactly one foreground #1/$1 mark; found {foreground_marks}"
        )
    return "\n".join(lines) + "\n", tips, foreground_marks, len(tips)


def set_control_value(text: str, key: str, value: str) -> str:
    pattern = re.compile(
        rf"^(?P<indent>\s*){re.escape(key)}\s*=.*$",
        flags=re.MULTILINE | re.IGNORECASE,
    )
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise ValueError(
            f"expected one active {key}= line in codeml.ctl template; found {len(matches)}"
        )
    return pattern.sub(lambda match: f"{match.group('indent')}{key} = {value}", text)


def render_control(
    template: str,
    *,
    seqfile: str,
    outfile: str,
    fix_omega: int,
    omega: str,
    cleandata: int,
) -> str:
    values = {
        "seqfile": seqfile,
        "treefile": "arbol.tree",
        "outfile": outfile,
        "model": "2",
        "NSsites": "2",
        "fix_omega": str(fix_omega),
        "omega": omega,
        "cleandata": str(cleandata),
    }
    rendered = template
    for key in CONTROL_KEYS:
        rendered = set_control_value(rendered, key, values[key])
    return rendered.rstrip() + "\n"


def parse_control(path: Path) -> dict[str, str]:
    observed: dict[str, str] = {}
    key_lookup = {key.lower(): key for key in CONTROL_KEYS}
    pattern = re.compile(r"^\s*([A-Za-z_]+)\s*=\s*([^\s*]+)")
    for raw in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(raw)
        if not match:
            continue
        normalized = match.group(1).lower()
        if normalized in key_lookup:
            canonical = key_lookup[normalized]
            if canonical in observed:
                raise ValueError(f"duplicate active {canonical}= in {path}")
            observed[canonical] = match.group(2)
    missing = [key for key in CONTROL_KEYS if key not in observed]
    if missing:
        raise ValueError(f"missing control keys in {path}: {', '.join(missing)}")
    return observed


def verify_control(path: Path, expected: dict[str, str]) -> None:
    observed = parse_control(path)
    differences = [
        f"{key}={observed.get(key)!r}, expected {value!r}"
        for key, value in expected.items()
        if observed.get(key) != value
    ]
    if differences:
        raise ValueError(f"invalid {path}: " + "; ".join(differences))


def save_result(path: Path, result: Result) -> None:
    atomic_text(path, json.dumps(asdict(result), indent=2, ensure_ascii=False) + "\n")


def load_result(path: Path) -> Result:
    return Result(**json.loads(path.read_text(encoding="utf-8")))


def write_summary(path: Path, results: list[Result]) -> None:
    names = [field.name for field in fields(Result)]
    rows = ["\t".join(names)]
    for result in sorted(results, key=lambda item: natural_key(item.og)):
        data = asdict(result)
        rows.append(
            "\t".join(
                str(data[name]).replace("\t", " ").replace("\r", " ").replace("\n", " ")
                for name in names
            )
        )
    atomic_text(path, "\n".join(rows) + "\n")


def process_og(
    og: str,
    codon_root: Path,
    template_root: Path,
    output_root: Path,
    cleandata: int,
) -> Result:
    source_dir = codon_root / og
    marker = source_dir / "CODON_FILTER_OK"
    if not marker.is_file():
        raise ValueError(f"missing CODON_FILTER_OK: {marker}")
    fasta = source_dir / f"{og}_BestMSA_filtered_codon.fasta"
    if not fasta.is_file() or fasta.stat().st_size == 0:
        raise ValueError(f"missing or empty filtered codon FASTA: {fasta}")

    records = read_fasta(fasta)
    metrics = validate_alignment(records)
    ids = {header for header, _ in records}

    old_alt = template_root / "alternative" / og
    old_null = template_root / "null" / og
    template_ctl_path = old_alt / "codeml.ctl"
    alt_tree_path = old_alt / "arbol.tree"
    null_tree_path = old_null / "arbol.tree"
    for required in (template_ctl_path, alt_tree_path, null_tree_path):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"missing or empty template input: {required}")

    alt_tree_text, tree_tips, foreground_marks, n_tree_tips = read_tree(alt_tree_path)
    null_tree_text, null_tips, null_marks, _ = read_tree(null_tree_path)
    if alt_tree_text != null_tree_text or tree_tips != null_tips or foreground_marks != null_marks:
        raise ValueError("previous alternative and null trees differ")
    if ids != tree_tips:
        missing_tree = sorted(ids - tree_tips)
        extra_tree = sorted(tree_tips - ids)
        details = []
        if missing_tree:
            details.append("alignment-only IDs: " + ", ".join(missing_tree))
        if extra_tree:
            details.append("tree-only IDs: " + ", ".join(extra_tree))
        raise ValueError("tree/alignment ID mismatch; " + "; ".join(details))

    alt_dir = output_root / "alternative" / og
    null_dir = output_root / "null" / og
    alt_dir.mkdir(parents=True, exist_ok=True)
    null_dir.mkdir(parents=True, exist_ok=True)

    seqfile = f"{og}_G3_LINSI.paml"
    outfile = f"{og}.output"
    paml_text = render_paml(records)
    template = template_ctl_path.read_text(encoding="utf-8")
    alt_ctl = render_control(
        template,
        seqfile=seqfile,
        outfile=outfile,
        fix_omega=0,
        omega=".4",
        cleandata=cleandata,
    )
    null_ctl = render_control(
        template,
        seqfile=seqfile,
        outfile=outfile,
        fix_omega=1,
        omega="1",
        cleandata=cleandata,
    )

    for directory, control in ((alt_dir, alt_ctl), (null_dir, null_ctl)):
        atomic_text(directory / seqfile, paml_text)
        atomic_text(directory / "arbol.tree", alt_tree_text)
        atomic_text(directory / "codeml.ctl", control)

    common_expected = {
        "seqfile": seqfile,
        "treefile": "arbol.tree",
        "outfile": outfile,
        "model": "2",
        "NSsites": "2",
        "cleandata": str(cleandata),
    }
    verify_control(
        alt_dir / "codeml.ctl",
        {**common_expected, "fix_omega": "0", "omega": ".4"},
    )
    verify_control(
        null_dir / "codeml.ctl",
        {**common_expected, "fix_omega": "1", "omega": "1"},
    )

    result = Result(
        og=og,
        status="OK",
        tree_tips=n_tree_tips,
        foreground_marks=foreground_marks,
        seqfile=seqfile,
        alternative_dir=str(alt_dir),
        null_dir=str(null_dir),
        **metrics,
    )
    status_dir = output_root / "status"
    save_result(status_dir / f"{og}.json", result)
    marker_text = (
        f"OK\tog={og}\tsequences={result.n_sequences}\t"
        f"codons={result.codon_columns}\tcomplete={result.complete_codon_columns}\t"
        f"cleandata={cleandata}\n"
    )
    atomic_text(alt_dir / "PAML_INPUT_OK", marker_text)
    atomic_text(null_dir / "PAML_INPUT_OK", marker_text)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build and audit GUIDANCE3 branch-site codeml inputs."
    )
    parser.add_argument("--codon-root", required=True, type=Path)
    parser.add_argument("--template-paml-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--cleandata", choices=(0, 1), type=int, default=1)
    parser.add_argument("--og", action="append", help="Process only this OG; repeatable.")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rewrite OGs already marked PAML_INPUT_OK. Do not use after codeml starts.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    codon_root = args.codon_root.resolve()
    template_root = args.template_paml_root.resolve()
    output_root = args.output_root.resolve()
    status_dir = output_root / "status"
    summary_path = output_root / "guidance3_paml_input_summary.tsv"
    run_log_path = output_root / "prepare_guidance3_paml.log"
    ready_marker = output_root / "PAML_INPUTS_READY"

    for directory in (codon_root, template_root):
        if not directory.is_dir():
            raise SystemExit(f"ERROR: directory not found: {directory}")
    output_root.mkdir(parents=True, exist_ok=True)
    status_dir.mkdir(parents=True, exist_ok=True)
    ready_marker.unlink(missing_ok=True)

    available = sorted(
        [directory.name for directory in codon_root.glob("OG*") if directory.is_dir()],
        key=natural_key,
    )
    selected = sorted(args.og or available, key=natural_key)
    if not selected:
        raise SystemExit(f"ERROR: no OG directories found in {codon_root}")

    results: list[Result] = []

    def log(message: str) -> None:
        print(message, flush=True)
        with run_log_path.open("a", encoding="utf-8") as handle:
            handle.write(message + "\n")

    log(
        f"[START] OGs={len(selected)} cleandata={args.cleandata} "
        f"output={output_root}"
    )

    for index, og in enumerate(selected, 1):
        alt_marker = output_root / "alternative" / og / "PAML_INPUT_OK"
        null_marker = output_root / "null" / og / "PAML_INPUT_OK"
        result_path = status_dir / f"{og}.json"

        if (
            not args.force
            and alt_marker.is_file()
            and null_marker.is_file()
            and result_path.is_file()
        ):
            result = load_result(result_path)
            results.append(result)
            log(f"[SKIP {index}/{len(selected)}] {og} already prepared")
            write_summary(summary_path, results)
            continue

        try:
            result = process_og(
                og,
                codon_root,
                template_root,
                output_root,
                args.cleandata,
            )
            log(
                f"[OK {index}/{len(selected)}] {og} "
                f"sequences={result.n_sequences} codons={result.codon_columns} "
                f"complete={result.complete_codon_columns}"
            )
        except Exception as exc:
            alt_marker.unlink(missing_ok=True)
            null_marker.unlink(missing_ok=True)
            result = Result(og=og, status="ERROR", error=str(exc))
            save_result(result_path, result)
            log(f"[ERROR {index}/{len(selected)}] {og}: {exc}")

        results.append(result)
        write_summary(summary_path, results)

    complete = sum(result.status == "OK" for result in results)
    errors = len(results) - complete
    log(f"[FINISH] total={len(results)} complete={complete} errors={errors}")
    log(f"[SUMMARY] {summary_path}")
    if errors == 0:
        atomic_text(
            ready_marker,
            f"READY\tOGs={complete}\tmodels={complete * 2}\t"
            f"cleandata={args.cleandata}\n",
        )
    return int(errors > 0)


if __name__ == "__main__":
    raise SystemExit(main())
