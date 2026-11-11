#!/usr/bin/env python3
"""Prepare GUIDANCE3 BestMSA codon alignments for codeml.

For each ``OG*_guidance`` directory, the program:

1. renames BestMSA protein headers using the final ``[species]`` tag;
2. orders the corresponding unaligned CDS by that species identifier;
3. runs PAL2NAL on the complete (unfiltered) BestMSA;
4. removes codon triplets whose GUIDANCE3 residue-pair column score is below
   the requested cutoff, while retaining unscored columns; and
5. masks surviving stop codons as ``NNN`` and validates frame, gap structure,
   identifiers, and lengths.

If a raw ``cds_species.fasta`` record has a length that is not divisible by
three, the corresponding ungapped record from the previously generated
``full_codon.fasta`` is used as a documented fallback.

Successful OGs receive a ``CODON_FILTER_OK`` marker.  A resumable TSV summary
is rewritten after every OG.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Iterable


STOP_CODONS = {"TAA", "TAG", "TGA"}
PERL_ENV_VARS = {
    "PERL5LIB",
    "PERL5OPT",
    "PERL_LOCAL_LIB_ROOT",
    "PERL_MB_OPT",
    "PERL_MM_OPT",
}


@dataclass
class Result:
    og: str
    status: str
    error: str = ""
    n_sequences: int | str = ""
    aa_columns_original: int | str = ""
    scored_columns: int | str = ""
    score_ge_cutoff: int | str = ""
    score_lt_cutoff: int | str = ""
    unscored_columns_kept: int | str = ""
    aa_columns_final: int | str = ""
    nt_columns_final: int | str = ""
    partial_gap_codons: int | str = ""
    internal_stops: int | str = ""
    terminal_stops: int | str = ""
    stop_codons_masked: int | str = ""
    ambiguous_codons: int | str = ""
    empty_sequences: int | str = ""
    cds_fallback_sequences: int | str = ""
    bestmsa: str = ""
    column_scores: str = ""
    source_cds: str = ""
    filtered_codon_alignment: str = ""
    stop_mask_audit: str = ""


def natural_key(value: str) -> tuple[int, str]:
    match = re.fullmatch(r"OG(\d+)", value)
    if match:
        return int(match.group(1)), value
    return sys.maxsize, value


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


def write_fasta(records: Iterable[tuple[str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


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


def unique_glob(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise ValueError(
            f"expected exactly one {pattern!r} in {directory}; found {len(matches)}"
        )
    if not matches[0].is_file() or matches[0].stat().st_size == 0:
        raise ValueError(f"missing or empty input: {matches[0]}")
    return matches[0]


def aligned_length(records: list[tuple[str, str]], label: str) -> int:
    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise ValueError(f"{label} contains unequal sequence lengths: {sorted(lengths)}")
    return lengths.pop()


def parse_bestmsa(
    path: Path,
) -> tuple[list[tuple[str, str]], list[str]]:
    original = read_fasta(path)
    renamed: list[tuple[str, str]] = []
    species_order: list[str] = []

    for header, sequence in original:
        match = re.search(r"\[([^\[\]]+)\]\s*$", header)
        if not match:
            raise ValueError(f"cannot extract final [species] from BestMSA header: {header}")
        species = match.group(1).strip()
        if not species:
            raise ValueError(f"empty species tag in BestMSA header: {header}")
        species_order.append(species)
        renamed.append((species, sequence))

    duplicates = sorted({name for name in species_order if species_order.count(name) > 1})
    if duplicates:
        raise ValueError("duplicate BestMSA species: " + ", ".join(duplicates))
    aligned_length(renamed, "BestMSA")
    return renamed, species_order


def order_cds(
    path: Path,
    species_order: list[str],
    fallback_path: Path,
) -> tuple[list[tuple[str, str]], list[str]]:
    cds_by_species: dict[str, str] = {}
    for header, sequence in read_fasta(path):
        species = header.split()[0]
        if species in cds_by_species:
            raise ValueError(f"duplicate CDS species: {species}")
        if "-" in sequence or "." in sequence:
            raise ValueError(f"source CDS contains alignment gaps: {species}")
        cds_by_species[species] = sequence

    fallback_by_species: dict[str, str] = {}
    if fallback_path.is_file() and fallback_path.stat().st_size > 0:
        for header, sequence in read_fasta(fallback_path):
            species = header.split()[0]
            if species in fallback_by_species:
                raise ValueError(f"duplicate fallback CDS species: {species}")
            fallback_by_species[species] = sequence.replace("-", "").replace(".", "")

    expected = set(species_order)
    observed = set(cds_by_species)
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    if missing or extra:
        details: list[str] = []
        if missing:
            details.append("missing CDS: " + ", ".join(missing))
        if extra:
            details.append("extra CDS: " + ", ".join(extra))
        raise ValueError("; ".join(details))

    ordered: list[tuple[str, str]] = []
    fallback_species: list[str] = []

    for species in species_order:
        sequence = cds_by_species[species]
        if len(sequence) % 3:
            fallback = fallback_by_species.get(species, "")
            if not fallback:
                raise ValueError(
                    f"CDS length is not divisible by three and no fallback exists: "
                    f"{species} ({len(sequence)} nt)"
                )
            if len(fallback) % 3:
                raise ValueError(
                    f"fallback CDS length is not divisible by three: "
                    f"{species} ({len(fallback)} nt)"
                )
            sequence = fallback
            fallback_species.append(species)
        ordered.append((species, sequence))

    return ordered, fallback_species


def parse_scores(path: Path, alignment_length: int) -> dict[int, float]:
    scores: dict[int, float] = {}
    found_end = False

    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if line == "#END":
                found_end = True
                continue
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"invalid score row at {path}:{line_number}")
            try:
                position = int(parts[0])
                score = float(parts[1])
            except ValueError as exc:
                raise ValueError(
                    f"invalid score row at {path}:{line_number}: {line}"
                ) from exc
            if position in scores:
                raise ValueError(f"duplicate score for column {position}")
            if not 1 <= position <= alignment_length:
                raise ValueError(
                    f"score column {position} outside alignment length {alignment_length}"
                )
            if not math.isfinite(score):
                raise ValueError(f"non-finite score for column {position}")
            scores[position] = score

    if not found_end:
        raise ValueError(f"score file lacks final #END marker: {path}")
    return scores


def run_pal2nal(
    perl: Path,
    pal2nal: Path,
    protein_path: Path,
    cds_path: Path,
    output_path: Path,
    stderr_path: Path,
) -> None:
    environment = os.environ.copy()
    for variable in PERL_ENV_VARS:
        environment.pop(variable, None)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = output_path.with_suffix(output_path.suffix + ".tmp")

    with temporary_output.open("w", encoding="utf-8", newline="\n") as stdout, (
        stderr_path.open("w", encoding="utf-8", newline="\n")
    ) as stderr:
        completed = subprocess.run(
            [
                str(perl),
                str(pal2nal),
                str(protein_path),
                str(cds_path),
                "-output",
                "fasta",
            ],
            stdout=stdout,
            stderr=stderr,
            env=environment,
            check=False,
        )

    if completed.returncode != 0 or temporary_output.stat().st_size == 0:
        temporary_output.unlink(missing_ok=True)
        message = stderr_path.read_text(encoding="utf-8", errors="replace").strip()
        if len(message) > 1200:
            message = message[-1200:]
        raise RuntimeError(
            f"PAL2NAL exit={completed.returncode}; {message or 'empty output'}"
        )

    os.replace(temporary_output, output_path)


def validate_codon_alignment(
    protein: list[tuple[str, str]],
    codon: list[tuple[str, str]],
    *,
    reject_stops: bool = True,
) -> dict[str, int]:
    protein_ids = [header.split()[0] for header, _ in protein]
    codon_ids = [header.split()[0] for header, _ in codon]
    if protein_ids != codon_ids:
        raise ValueError("protein and codon identifiers/order differ")

    aa_length = aligned_length(protein, "protein alignment")
    nt_length = aligned_length(codon, "codon alignment")
    if nt_length != aa_length * 3:
        raise ValueError(
            f"codon alignment has {nt_length} nt; expected {aa_length * 3}"
        )

    partial_gaps = 0
    internal_stops = 0
    terminal_stops = 0
    ambiguous_codons = 0
    empty_sequences = 0

    for species, sequence in codon:
        nongap_codons: list[str] = []
        for start in range(0, len(sequence), 3):
            triplet = sequence[start : start + 3]
            if "-" in triplet and triplet != "---":
                partial_gaps += 1
            if triplet != "---" and "-" not in triplet:
                nongap_codons.append(triplet)
                if any(base not in "ACGT" for base in triplet):
                    ambiguous_codons += 1

        if not nongap_codons:
            empty_sequences += 1
            continue
        internal_stops += sum(codon_ in STOP_CODONS for codon_ in nongap_codons[:-1])
        terminal_stops += int(nongap_codons[-1] in STOP_CODONS)

    fatal = []
    if partial_gaps:
        fatal.append(f"{partial_gaps} partial-gap codons")
    if reject_stops and internal_stops:
        fatal.append(f"{internal_stops} internal stop codons")
    if reject_stops and terminal_stops:
        fatal.append(f"{terminal_stops} terminal stop codons")
    if empty_sequences:
        fatal.append(f"{empty_sequences} all-gap sequences")
    if fatal:
        raise ValueError("; ".join(fatal))

    return {
        "partial_gap_codons": partial_gaps,
        "internal_stops": internal_stops,
        "terminal_stops": terminal_stops,
        "ambiguous_codons": ambiguous_codons,
        "empty_sequences": empty_sequences,
    }


def mask_stop_codons(
    protein: list[tuple[str, str]],
    codon: list[tuple[str, str]],
) -> tuple[
    list[tuple[str, str]],
    list[tuple[str, str]],
    list[tuple[str, int, str, str]],
]:
    masked_protein: list[tuple[str, str]] = []
    masked_codon: list[tuple[str, str]] = []
    audit: list[tuple[str, int, str, str]] = []

    for (protein_id, protein_sequence), (codon_id, codon_sequence) in zip(
        protein, codon
    ):
        if protein_id != codon_id:
            raise ValueError("protein and codon identifiers differ while masking stops")

        protein_chars = list(protein_sequence)
        codon_chars = list(codon_sequence)

        for column in range(1, len(protein_sequence) + 1):
            start = (column - 1) * 3
            triplet = codon_sequence[start : start + 3]
            if triplet in STOP_CODONS:
                original_symbol = protein_chars[column - 1]
                codon_chars[start : start + 3] = list("NNN")
                protein_chars[column - 1] = "X"
                audit.append((protein_id, column, triplet, original_symbol))

        masked_protein.append((protein_id, "".join(protein_chars)))
        masked_codon.append((codon_id, "".join(codon_chars)))

    return masked_protein, masked_codon, audit


def write_stop_audit(
    path: Path, audit: list[tuple[str, int, str, str]]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["species", "filtered_column", "original_codon", "protein_symbol", "action"]
        )
        for species, column, codon, symbol in audit:
            writer.writerow([species, column, codon, symbol, "MASK_AS_NNN"])


def write_column_audit(
    path: Path, alignment_length: int, scores: dict[int, float], cutoff: float
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["original_column", "score", "action"])
        for position in range(1, alignment_length + 1):
            if position not in scores:
                writer.writerow([position, "NA", "KEEP_UNSCORED"])
            elif scores[position] < cutoff:
                writer.writerow([position, f"{scores[position]:.6f}", "REMOVE_LOW_SCORE"])
            else:
                writer.writerow([position, f"{scores[position]:.6f}", "KEEP"])


def write_summary(path: Path, results: list[Result]) -> None:
    column_names = [field.name for field in fields(Result)]
    rows: list[str] = []
    rows.append("\t".join(column_names))
    for result in sorted(results, key=lambda item: natural_key(item.og)):
        values = []
        data = asdict(result)
        for name in column_names:
            value = str(data[name]).replace("\t", " ").replace("\r", " ").replace("\n", " ")
            values.append(value)
        rows.append("\t".join(values))
    atomic_text(path, "\n".join(rows) + "\n")


def save_result_json(path: Path, result: Result) -> None:
    atomic_text(path, json.dumps(asdict(result), indent=2, ensure_ascii=False) + "\n")


def load_result_json(path: Path) -> Result:
    return Result(**json.loads(path.read_text(encoding="utf-8")))


def process_og(
    og: str,
    guidance_dir: Path,
    source_codon_root: Path,
    output_root: Path,
    cutoff: float,
    perl: Path,
    pal2nal: Path,
) -> Result:
    result = Result(og=og, status="ERROR")
    output_dir = output_root / og
    output_dir.mkdir(parents=True, exist_ok=True)

    bestmsa = unique_glob(guidance_dir, "*BestMSA.fasta")
    score_path = unique_glob(guidance_dir, "*BestMSA_res_pair_col.scr")
    cds_path = source_codon_root / og / "cds_species.fasta"
    fallback_cds_path = source_codon_root / og / "full_codon.fasta"
    if not cds_path.is_file() or cds_path.stat().st_size == 0:
        raise ValueError(f"missing or empty source CDS: {cds_path}")

    result.bestmsa = str(bestmsa)
    result.column_scores = str(score_path)
    result.source_cds = str(cds_path)

    protein_full, species_order = parse_bestmsa(bestmsa)
    cds_ordered, fallback_species = order_cds(
        cds_path, species_order, fallback_cds_path
    )
    aa_length = aligned_length(protein_full, "BestMSA")
    scores = parse_scores(score_path, aa_length)

    protein_full_path = output_dir / f"{og}_BestMSA_full_species.fasta"
    cds_ordered_path = output_dir / f"{og}_CDS_ordered.fasta"
    codon_full_path = output_dir / f"{og}_BestMSA_full_codon.fasta"
    pal2nal_stderr = output_dir / f"{og}_pal2nal.stderr"
    protein_filtered_path = output_dir / f"{og}_BestMSA_filtered.fasta"
    codon_filtered_path = output_dir / f"{og}_BestMSA_filtered_codon.fasta"
    column_audit_path = output_dir / f"{og}_column_filter.tsv"
    stop_audit_path = output_dir / f"{og}_stop_mask.tsv"

    write_fasta(protein_full, protein_full_path)
    write_fasta(cds_ordered, cds_ordered_path)
    run_pal2nal(
        perl,
        pal2nal,
        protein_full_path,
        cds_ordered_path,
        codon_full_path,
        pal2nal_stderr,
    )

    codon_full = read_fasta(codon_full_path)
    validate_codon_alignment(protein_full, codon_full, reject_stops=False)

    removed = {position for position, score in scores.items() if score < cutoff}
    kept = [position for position in range(1, aa_length + 1) if position not in removed]
    unscored = set(range(1, aa_length + 1)) - set(scores)

    protein_filtered = [
        (
            species,
            "".join(sequence[position - 1] for position in kept),
        )
        for species, sequence in protein_full
    ]
    codon_filtered = [
        (
            species,
            "".join(sequence[(position - 1) * 3 : position * 3] for position in kept),
        )
        for species, sequence in codon_full
    ]

    protein_filtered, codon_filtered, stop_audit = mask_stop_codons(
        protein_filtered, codon_filtered
    )
    metrics = validate_codon_alignment(
        protein_filtered, codon_filtered, reject_stops=True
    )
    write_fasta(protein_filtered, protein_filtered_path)
    write_fasta(codon_filtered, codon_filtered_path)
    write_column_audit(column_audit_path, aa_length, scores, cutoff)
    write_stop_audit(stop_audit_path, stop_audit)

    result.status = "OK"
    result.n_sequences = len(protein_full)
    result.aa_columns_original = aa_length
    result.scored_columns = len(scores)
    result.score_lt_cutoff = len(removed)
    result.score_ge_cutoff = len(scores) - len(removed)
    result.unscored_columns_kept = len(unscored)
    result.aa_columns_final = len(kept)
    result.nt_columns_final = len(kept) * 3
    result.partial_gap_codons = metrics["partial_gap_codons"]
    result.internal_stops = metrics["internal_stops"]
    result.terminal_stops = metrics["terminal_stops"]
    result.stop_codons_masked = len(stop_audit)
    result.ambiguous_codons = metrics["ambiguous_codons"]
    result.empty_sequences = metrics["empty_sequences"]
    result.cds_fallback_sequences = len(fallback_species)
    result.filtered_codon_alignment = str(codon_filtered_path)
    result.stop_mask_audit = str(stop_audit_path)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Project GUIDANCE3 BestMSAs to CDS with PAL2NAL, filter codon "
            "triplets by column score, and audit every OG."
        )
    )
    parser.add_argument("--guidance-root", required=True, type=Path)
    parser.add_argument("--source-codon-root", required=True, type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--cutoff", type=float, default=0.93)
    parser.add_argument("--pal2nal", type=Path, default=Path("/usr/bin/pal2nal.pl"))
    parser.add_argument("--perl", type=Path, default=Path("/usr/bin/perl"))
    parser.add_argument(
        "--og",
        action="append",
        help="Process only this OG; may be supplied more than once.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Reprocess OGs that already have CODON_FILTER_OK.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    guidance_root = args.guidance_root.resolve()
    source_codon_root = args.source_codon_root.resolve()
    output_root = (
        args.output_root.resolve()
        if args.output_root
        else guidance_root / "codon_results"
    )
    results_root = guidance_root / "results"

    if not 0 <= args.cutoff <= 1:
        raise SystemExit("ERROR: --cutoff must be between 0 and 1")
    for executable in (args.perl, args.pal2nal):
        if not executable.is_file():
            raise SystemExit(f"ERROR: executable/script not found: {executable}")
    if not results_root.is_dir():
        raise SystemExit(f"ERROR: results directory not found: {results_root}")
    if not source_codon_root.is_dir():
        raise SystemExit(f"ERROR: source codon root not found: {source_codon_root}")

    guidance_dirs: dict[str, Path] = {}
    for directory in results_root.glob("OG*_guidance"):
        if directory.is_dir():
            guidance_dirs[directory.name.removesuffix("_guidance")] = directory

    selected_ogs = sorted(args.og or guidance_dirs, key=natural_key)
    if not selected_ogs:
        raise SystemExit(f"ERROR: no OG*_guidance directories found in {results_root}")
    missing_dirs = [og for og in selected_ogs if og not in guidance_dirs]
    if missing_dirs:
        raise SystemExit("ERROR: GUIDANCE directories missing for: " + ", ".join(missing_dirs))

    output_root.mkdir(parents=True, exist_ok=True)
    summary_path = output_root / "guidance3_codon_summary.tsv"
    run_log_path = output_root / "guidance3_codon_run.log"
    run_results: list[Result] = []

    def log(message: str) -> None:
        print(message, flush=True)
        with run_log_path.open("a", encoding="utf-8") as handle:
            handle.write(message + "\n")

    log(
        f"[START] OGs={len(selected_ogs)} cutoff={args.cutoff} "
        f"output={output_root}"
    )

    for index, og in enumerate(selected_ogs, 1):
        output_dir = output_root / og
        marker = output_dir / "CODON_FILTER_OK"
        result_json = output_dir / "summary.json"

        if marker.is_file() and result_json.is_file() and not args.force:
            result = load_result_json(result_json)
            result.status = "OK"
            run_results.append(result)
            log(f"[SKIP {index}/{len(selected_ogs)}] {og} already complete")
            write_summary(summary_path, run_results)
            continue

        try:
            result = process_og(
                og=og,
                guidance_dir=guidance_dirs[og],
                source_codon_root=source_codon_root,
                output_root=output_root,
                cutoff=args.cutoff,
                perl=args.perl,
                pal2nal=args.pal2nal,
            )
            save_result_json(result_json, result)
            atomic_text(
                marker,
                f"OK\tog={og}\tcutoff={args.cutoff}\t"
                f"aa_final={result.aa_columns_final}\t"
                f"nt_final={result.nt_columns_final}\n",
            )
            log(
                f"[OK {index}/{len(selected_ogs)}] {og} "
                f"sequences={result.n_sequences} "
                f"aa={result.aa_columns_original}->{result.aa_columns_final}"
                f" fallbacks={result.cds_fallback_sequences}"
                f" stops_masked={result.stop_codons_masked}"
            )
        except Exception as exc:  # continue so all problematic OGs are reported
            marker.unlink(missing_ok=True)
            result = Result(og=og, status="ERROR", error=str(exc))
            save_result_json(result_json, result)
            log(f"[ERROR {index}/{len(selected_ogs)}] {og}: {exc}")

        run_results.append(result)
        write_summary(summary_path, run_results)

    completed = sum(result.status == "OK" for result in run_results)
    failed = len(run_results) - completed
    log(f"[FINISH] total={len(run_results)} complete={completed} errors={failed}")
    log(f"[SUMMARY] {summary_path}")
    return int(failed > 0)


if __name__ == "__main__":
    raise SystemExit(main())
