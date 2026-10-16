#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

BLAST_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "evalue",
    "bitscore",
    "qcovs",
]

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Screen an existing ASV FASTA against MBPD without rerunning DADA2."
    )
    parser.add_argument("--asv-fasta", required=True, type=Path)
    parser.add_argument("--abundance", type=Path, help="Optional CSV abundance table")
    parser.add_argument("--db-fasta", required=True, type=Path)
    parser.add_argument("--db-tax", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--identity", type=float, default=97.0)
    parser.add_argument("--coverage", type=float, default=90.0)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--max-target-seqs", type=int, default=50)
    return parser.parse_args()

def check_inputs(args: argparse.Namespace) -> None:
    for path in (args.asv_fasta, args.db_fasta, args.db_tax):
        if not path.is_file():
            raise FileNotFoundError(f"File not found: {path}")
    if args.abundance and not args.abundance.is_file():
        raise FileNotFoundError(f"File not found: {args.abundance}")
    if args.outdir.exists() and (not args.outdir.is_dir() or any(args.outdir.iterdir())):
        raise ValueError(f"Output path must be new or empty: {args.outdir}")
    for program in ("makeblastdb", "blastn"):
        if shutil.which(program) is None:
            raise RuntimeError(f"{program} is not available in PATH")
    if not 0 <= args.identity <= 100 or not 0 <= args.coverage <= 100:
        raise ValueError("identity and coverage must be between 0 and 100")
    if args.threads < 1 or args.max_target_seqs < 1:
        raise ValueError("threads and max-target-seqs must be positive")

def read_fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.startswith(">"):
                continue
            parts = line[1:].strip().split()
            if not parts:
                raise ValueError(f"Empty FASTA identifier at line {line_no}")
            seq_id = parts[0]
            if seq_id in seen:
                raise ValueError(f"Duplicate FASTA identifier: {seq_id}")
            seen.add(seq_id)
            ids.append(seq_id)
    if not ids:
        raise ValueError(f"No FASTA records found in {path}")
    return ids

def read_taxonomy(path: Path) -> dict[str, str]:
    taxonomy: dict[str, str] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            cells = line.rstrip("\r\n").split("\t")
            if not cells or not cells[0]:
                continue
            subject_id = cells[0].split()[0]
            annotation = ";".join(cell.strip() for cell in cells[1:] if cell.strip())
            taxonomy.setdefault(subject_id, annotation or "unannotated")
    return taxonomy

def build_database(db_fasta: Path, db_prefix: Path) -> None:
    required = [Path(f"{db_prefix}.{suffix}") for suffix in ("nhr", "nin", "nsq")]
    if all(path.exists() for path in required):
        return
    subprocess.run(
        [
            "makeblastdb",
            "-in",
            str(db_fasta),
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-out",
            str(db_prefix),
        ],
        check=True,
    )

def run_blast(args: argparse.Namespace, db_prefix: Path, raw_output: Path) -> None:
    outfmt = "6 " + " ".join(BLAST_FIELDS)
    subprocess.run(
        [
            "blastn",
            "-task",
            "megablast",
            "-query",
            str(args.asv_fasta),
            "-db",
            str(db_prefix),
            "-strand",
            "both",
            "-evalue",
            "1e-20",
            "-perc_identity",
            str(args.identity),
            "-qcov_hsp_perc",
            str(args.coverage),
            "-max_target_seqs",
            str(args.max_target_seqs),
            "-num_threads",
            str(args.threads),
            "-outfmt",
            outfmt,
            "-out",
            str(raw_output),
        ],
        check=True,
    )

def parse_hits(path: Path, taxonomy: dict[str, str]) -> dict[str, list[dict[str, object]]]:
    hits: dict[str, list[dict[str, object]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=BLAST_FIELDS)
        for row in reader:
            hit: dict[str, object] = dict(row)
            for field in ("pident", "evalue", "bitscore", "qcovs"):
                hit[field] = float(row[field])
            for field in ("length", "qlen"):
                hit[field] = int(row[field])
            hit["taxonomy"] = taxonomy.get(row["sseqid"], "taxonomy_not_found")
            hits[row["qseqid"]].append(hit)
    return hits

def write_all_hits(path: Path, hits: dict[str, list[dict[str, object]]]) -> None:
    fields = BLAST_FIELDS + ["taxonomy"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for query_id in sorted(hits):
            ordered = sorted(
                hits[query_id],
                key=lambda h: (-float(h["bitscore"]), -float(h["pident"]), -float(h["qcovs"])),
            )
            writer.writerows(ordered)

def summarize_hits(
    path: Path, asv_ids: list[str], hits: dict[str, list[dict[str, object]]]
) -> set[str]:
    candidate_ids: set[str] = set()
    fields = [
        "ASV_ID",
        "status",
        "best_identity",
        "best_query_coverage",
        "best_bitscore",
        "best_subject_count",
        "best_subject_ids",
        "best_taxonomy",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for asv_id in asv_ids:
            query_hits = hits.get(asv_id, [])
            if not query_hits:
                writer.writerow({"ASV_ID": asv_id, "status": "no_hit"})
                continue
            candidate_ids.add(asv_id)
            best_score = max(float(hit["bitscore"]) for hit in query_hits)
            best = [
                hit
                for hit in query_hits
                if abs(float(hit["bitscore"]) - best_score) < 1e-9
            ]
            annotations = sorted({str(hit["taxonomy"]) for hit in best})
            status = "unique_best_annotation" if len(annotations) == 1 else "ambiguous_best_hit"
            writer.writerow(
                {
                    "ASV_ID": asv_id,
                    "status": status,
                    "best_identity": max(float(hit["pident"]) for hit in best),
                    "best_query_coverage": max(float(hit["qcovs"]) for hit in best),
                    "best_bitscore": best_score,
                    "best_subject_count": len(best),
                    "best_subject_ids": ";".join(sorted({str(hit["sseqid"]) for hit in best})),
                    "best_taxonomy": " || ".join(annotations),
                }
            )
    return candidate_ids

def filter_abundance(source: Path, destination: Path, candidate_ids: set[str]) -> None:
    with source.open(newline="", encoding="utf-8-sig") as src, destination.open(
        "w", newline="", encoding="utf-8"
    ) as dst:
        reader = csv.reader(src)
        writer = csv.writer(dst)
        header = next(reader)
        writer.writerow(header)
        for row in reader:
            if row and row[0] in candidate_ids:
                writer.writerow(row)

def main() -> int:
    args = parse_args()
    try:
        check_inputs(args)
        asv_ids = read_fasta_ids(args.asv_fasta)
        taxonomy = read_taxonomy(args.db_tax)
        args.outdir.mkdir(parents=True, exist_ok=True)
        db_prefix = args.outdir / "mbpd_pathogen_blastdb"
        raw_output = args.outdir / "01_blast_all_hits.tsv"
        build_database(args.db_fasta, db_prefix)
        run_blast(args, db_prefix, raw_output)
        hits = parse_hits(raw_output, taxonomy)
        write_all_hits(args.outdir / "02_filtered_hits_with_taxonomy.tsv", hits)
        candidate_ids = summarize_hits(
            args.outdir / "03_asv_pathogen_candidate_summary.tsv", asv_ids, hits
        )
        if args.abundance:
            filter_abundance(
                args.abundance,
                args.outdir / "04_pathogen_candidate_abundance.csv",
                candidate_ids,
            )
        print(f"ASVs screened: {len(asv_ids)}")
        print(f"Candidate ASVs: {len(candidate_ids)}")
        print(f"Results: {args.outdir.resolve()}")
        return 0
    except (FileNotFoundError, RuntimeError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    raise SystemExit(main())
