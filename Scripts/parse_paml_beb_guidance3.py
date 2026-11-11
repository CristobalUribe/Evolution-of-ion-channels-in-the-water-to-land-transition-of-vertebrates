#!/usr/bin/env python3
"""Extract MLE tables and starred BEB sites from the GUIDANCE3 PAML run.

By default the script reads only the alternative branch-site models stored in
``Guidance3_LINSI_all167/PAML_GUIDANCE3/alternative``.  This deliberately
excludes the null models.  The fixed list contains the 39 OGs that passed the
multiple-testing criterion in the GUIDANCE3 analysis.

The output is a semicolon-delimited CSV for straightforward import into Excel.
Its first line is ``sep=;``, which tells Excel to use semicolons even when the
computer's regional list separator is different.  Two intentionally empty
columns visually separate the identifiers, MLE table, and BEB results.  Each
OG uses a four-row MLE block (site class, proportion, background w, and
foreground w), matching the layout printed by PAML.  Each BEB cell contains at
most three sites, separated by `` | ``.

The reported BEB coordinates are positions in the GUIDANCE3-filtered PAML
alignment, exactly as printed by CODEML; they are not positions in an unaligned
reference protein or in the original unfiltered BestMSA.

Example
-------
python3 parse_paml_beb_guidance3.py
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path


GENES = OrderedDict(
    [
        ("OG6", "PKD1"),
        ("OG18", "PKDREJ"),
        ("OG148", "TRPM4"),
        ("OG153", "TRPA1"),
        ("OG563", "CHRND"),
        ("OG633", "GABRA6"),
        ("OG636", "KCNJ14"),
        ("OG353", "TPC3"),
        ("OG117", "GRIN2C"),
        ("OG510", "MCOLN1"),
        ("OG19", "PIEZO2"),
        ("OG286", "TRPV1"),
        ("OG805", "KCNJ8"),
        ("OG287", "TMEM63C"),
        ("OG240", "GRIA1"),
        ("OG613", "KCNA2"),
        ("OG646", "KCNK5"),
        ("OG238", "KCNQ2"),
        ("OG335", "TRPV3"),
        ("OG494", "SCNN1G"),
        ("OG527", "KCNK10"),
        ("OG163", "ANO3"),
        ("OG27", "CACNA1E"),
        ("OG154", "TRPM8"),
        ("OG1", "RYR2"),
        ("OG856", "KCNJ11"),
        ("OG560", "CHRNB2"),
        ("OG99", "TRPN1"),
        ("OG91", "GRIN2B"),
        ("OG642", "GABRR1"),
        ("OG671", "KCNA7"),
        ("OG370", "TRPV6"),
        ("OG296", "CLCN7"),
        ("OG841", "KCNJ10"),
        ("OG213", "HCN3"),
        ("OG22", "PIEZO1"),
        ("OG151", "ANO4"),
        ("OG273", "LRRC8D"),
        ("OG209", "ANO5"),
    ]
)

DEFAULT_PAML_ROOT = Path(
    "/home/cloyola/PAML4.10/Guidance3_LINSI_all167/PAML_GUIDANCE3"
)
DEFAULT_ALTERNATIVE_ROOT = DEFAULT_PAML_ROOT / "alternative"
DEFAULT_OUTPUT = DEFAULT_PAML_ROOT / "PAML_GUIDANCE3_selected_39_BEB.csv"

SITE_CLASSES = ("0", "1", "2a", "2b")
TABLE_ROWS = ("proportion", "background_w", "foreground_w")
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"

MLE_PATTERN = re.compile(
    rf"MLEs of dN/dS \(w\) for site classes \(K=4\).*?"
    rf"^\s*site class\s+0\s+1\s+2a\s+2b\s*$\s*"
    rf"^\s*proportion\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*$\s*"
    rf"^\s*background w\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*$\s*"
    rf"^\s*foreground w\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*$",
    flags=re.MULTILINE | re.DOTALL,
)

BEB_HEADER = re.compile(
    r"^Bayes Empirical Bayes \(BEB\) analysis\b", flags=re.MULTILINE
)

BEB_SITE = re.compile(
    rf"^\s*(\d+)\s+([A-Za-z*?\-])\s+({NUMBER})(\*{{1,2}})(?:\s+.*)?$"
)


def natural_og_key(og: str) -> int:
    return int(og[2:])


def identify_og(path: Path) -> str | None:
    """Identify an OG without confusing OG1 with OG19 or OG100."""
    match = re.match(r"^(OG\d+)(?:[^0-9].*)?$", path.stem, flags=re.IGNORECASE)
    if match:
        return match.group(1).upper()

    for part in reversed(path.parts):
        match = re.match(r"^(OG\d+)(?:_\d+)?$", part, flags=re.IGNORECASE)
        if match:
            return match.group(1).upper()
    return None


def extract_mle_table(text: str) -> dict[str, str] | None:
    """Return the last complete K=4 table, preserving PAML's printed precision."""
    matches = list(MLE_PATTERN.finditer(text))
    if not matches:
        return None

    values = matches[-1].groups()
    result: dict[str, str] = {}
    offset = 0
    for row_name in TABLE_ROWS:
        for site_class in SITE_CLASSES:
            result[f"{row_name}_{site_class}"] = values[offset]
            offset += 1
    return result


def extract_starred_beb_sites(text: str) -> list[str] | None:
    """Extract starred sites from BEB only; return None if no BEB block exists."""
    headers = list(BEB_HEADER.finditer(text))
    if not headers:
        return None

    beb_text = text[headers[-1].end() :]
    positive_header = re.search(
        r"^Positive sites for foreground lineages Prob\(w>1\):\s*$",
        beb_text,
        flags=re.MULTILINE,
    )
    if not positive_header:
        return []

    lines = beb_text[positive_header.end() :].splitlines()
    sites: list[str] = []
    started = False

    for line in lines:
        stripped = line.strip()
        if not stripped:
            if started:
                break
            continue

        # Every BEB result line begins with a position.  Once another section
        # begins (for example, "The grid"), the BEB table has ended.
        if not re.match(r"^\d+\s+", stripped):
            if started:
                break
            continue

        started = True
        match = BEB_SITE.match(line)
        if match:
            position, amino_acid, probability, stars = match.groups()
            sites.append(f"{position} {amino_acid.upper()} {probability}{stars}")

    return sites


def read_candidate(path: Path) -> tuple[str, dict[str, str] | None, list[str] | None]:
    text = path.read_text(encoding="utf-8", errors="replace")
    return text, extract_mle_table(text), extract_starred_beb_sites(text)


def choose_alternative_output(
    og: str, candidates: list[Path]
) -> tuple[Path | None, dict[str, str] | None, list[str] | None, list[str]]:
    """Prefer the complete output containing BEB; use mtime only to break ties."""
    warnings: list[str] = []
    parsed = []

    for path in candidates:
        try:
            _, table, beb_sites = read_candidate(path)
        except OSError as exc:
            warnings.append(f"{og}: no se pudo leer {path}: {exc}")
            continue

        rank = (
            int(beb_sites is not None),
            int(table is not None),
            path.stat().st_mtime,
        )
        parsed.append((rank, path, table, beb_sites))

    if not parsed:
        return None, None, None, warnings

    parsed.sort(key=lambda item: item[0], reverse=True)
    _, path, table, beb_sites = parsed[0]

    equally_valid = [
        item[1]
        for item in parsed
        if item[2] is not None and item[3] is not None
    ]
    if len(equally_valid) > 1:
        warnings.append(
            f"{og}: se encontraron {len(equally_valid)} outputs alternativos con BEB; "
            f"se uso el mas reciente: {path}"
        )

    return path, table, beb_sites, warnings


def collect_outputs(root: Path, selected_ogs: set[str]) -> dict[str, list[Path]]:
    candidates: dict[str, list[Path]] = defaultdict(list)

    paths = [root] if root.is_file() else root.rglob("*.output")
    for path in paths:
        if not path.is_file():
            continue
        og = identify_og(path)
        if og in selected_ogs:
            candidates[og].append(path)

    return candidates


def chunks(items: list[str], size: int) -> list[str]:
    # The CSV itself uses semicolons as column delimiters.  Pipes keep the three
    # sites together inside one cell without conflicting with that delimiter.
    return [" | ".join(items[i : i + size]) for i in range(0, len(items), size)]


def parse_only(value: str | None) -> list[str]:
    if not value:
        return list(GENES)

    requested = [item.strip().upper() for item in value.split(",") if item.strip()]
    unknown = [og for og in requested if og not in GENES]
    if unknown:
        raise ValueError("OG no incluido en la lista de 39: " + ", ".join(unknown))
    return requested


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Extrae la tabla MLE y los sitios BEB con asterisco de los 39 OG "
            "seleccionados positivamente en la corrida GUIDANCE3."
        )
    )
    parser.add_argument(
        "root",
        nargs="?",
        type=Path,
        default=DEFAULT_ALTERNATIVE_ROOT,
        help=(
            "Directorio alternative de PAML_GUIDANCE3; por defecto: "
            f"{DEFAULT_ALTERNATIVE_ROOT}"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"CSV de salida (por defecto: {DEFAULT_OUTPUT}).",
    )
    parser.add_argument(
        "--group-size",
        type=int,
        default=3,
        help="Numero de sitios BEB por celda/columna (por defecto: 3).",
    )
    parser.add_argument(
        "--only",
        help="Solo para pruebas: OG separados por coma, por ejemplo --only OG1,OG2.",
    )
    args = parser.parse_args()

    if not args.root.exists():
        parser.error(f"no existe: {args.root}")
    if args.root.is_dir() and (args.root / "alternative").is_dir():
        args.root = args.root / "alternative"
    if "null" in {part.lower() for part in args.root.parts}:
        parser.error(
            "la raiz apunta a un directorio null; use PAML_GUIDANCE3/alternative"
        )
    if args.group_size < 1:
        parser.error("--group-size debe ser al menos 1")

    try:
        og_order = parse_only(args.only)
    except ValueError as exc:
        parser.error(str(exc))

    candidates = collect_outputs(args.root, set(og_order))
    rows = []
    messages: list[str] = []

    for og in og_order:
        path, table, beb_sites, warnings = choose_alternative_output(
            og, candidates.get(og, [])
        )
        messages.extend(warnings)

        if path is None:
            status = "OUTPUT_NO_ENCONTRADO"
            messages.append(f"{og}: no se encontro ningun archivo .output")
        elif table is None:
            status = "TABLA_MLE_NO_ENCONTRADA"
            messages.append(f"{og}: falta la tabla MLE en {path}")
        elif beb_sites is None:
            status = "BEB_NO_ENCONTRADO"
            messages.append(f"{og}: el output elegido no contiene el bloque BEB: {path}")
        else:
            status = "OK"

        row = {
            "OMAGroup": og,
            "Gene_Symbol": GENES[og],
            **{
                f"{row_name}_{site_class}": (table or {}).get(
                    f"{row_name}_{site_class}", ""
                )
                for row_name in TABLE_ROWS
                for site_class in SITE_CLASSES
            },
            "BEB_starred_count": len(beb_sites or []),
            "Status": status,
        }
        row["_beb_groups"] = chunks(beb_sites or [], args.group_size)
        rows.append(row)

    max_groups = max((len(row["_beb_groups"]) for row in rows), default=0)
    beb_headers = [
        f"BEB_sites_{i * args.group_size + 1}-{(i + 1) * args.group_size}"
        for i in range(max_groups)
    ]
    headers = (
        [
            "OMAGroup",
            "Gene_Symbol",
            "",
            "MLEs of dN/dS (w) for site classes (K=4)",
            "",
            "",
            "",
            "",
            "",
            "BEB_starred_count",
        ]
        + beb_headers
        + ["Status"]
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8-sig", newline="") as handle:
        # Excel-specific separator declaration.  This makes a double-clicked
        # CSV open in separate columns even if Windows is configured to use a
        # comma as its regional list separator.
        handle.write("sep=;\n")
        writer = csv.writer(handle, delimiter=";")
        writer.writerow(headers)
        for row in rows:
            groups = row["_beb_groups"]
            padded_groups = groups + [""] * (max_groups - len(groups))

            # First row of each OG: identifiers, PAML site classes, and BEB.
            writer.writerow(
                [
                    row["OMAGroup"],
                    row["Gene_Symbol"],
                    "",
                    "site class",
                    *SITE_CLASSES,
                    "",
                    row["BEB_starred_count"],
                    *padded_groups,
                    row["Status"],
                ]
            )

            # Remaining three rows reproduce the MLE table printed by PAML.
            for row_name, display_name in (
                ("proportion", "proportion"),
                ("background_w", "background w"),
                ("foreground_w", "foreground w"),
            ):
                writer.writerow(
                    [
                        "",
                        "",
                        "",
                        display_name,
                        *[row[f"{row_name}_{site_class}"] for site_class in SITE_CLASSES],
                        "",
                        "",
                        *([""] * max_groups),
                        "",
                    ]
                )

    ok_count = sum(row["Status"] == "OK" for row in rows)
    print(f"Raiz alternativa analizada: {args.root}")
    print(f"CSV creado: {args.output}")
    print(f"OG procesados correctamente: {ok_count}/{len(rows)}")
    if messages:
        print("\nADVERTENCIAS:", file=sys.stderr)
        for message in messages:
            print(f"- {message}", file=sys.stderr)

    return 0 if ok_count == len(rows) else 2


if __name__ == "__main__":
    raise SystemExit(main())
