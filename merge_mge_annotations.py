#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional
from urllib.parse import unquote

import pandas as pd


FA_EXTS = {".fa", ".fna", ".fasta", ".fas", ".fsa"}


def parse_gff_attributes(attr_str: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    if not attr_str:
        return attrs
    for item in attr_str.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k] = unquote(v)
        else:
            attrs[item] = ""
    return attrs



def read_table_auto(path: Path) -> pd.DataFrame:
    """Try to read a delimited table robustly."""
    for sep in ["\t", ",", None]:
        try:
            kwargs = {"comment": "#", "dtype": str, "keep_default_na": False}
            if sep is None:
                df = pd.read_csv(path, sep=None, engine="python", **kwargs)
            else:
                df = pd.read_csv(path, sep=sep, engine="python", **kwargs)
            if df.shape[1] >= 2:
                return df
        except Exception:
            continue
    raise ValueError(f"Could not parse table: {path}")



def normalize_mobileelementfinder(sample: str, sample_dir: Path) -> List[dict]:
    rows: List[dict] = []
    tool_dir = sample_dir / "mobileelementfinder"
    if not tool_dir.exists():
        return rows

    for gff in sorted(tool_dir.glob("*.gff")):
        with gff.open() as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 9:
                    continue
                seqid, source, feature_type, start, end, score, strand, phase, attrs = parts
                attrd = parse_gff_attributes(attrs)
                rows.append(
                    {
                        "sample": sample,
                        "tool": "MobileElementFinder",
                        "seqid": seqid,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "feature_type": feature_type,
                        "mge_name": attrd.get("Name") or attrd.get("ID") or feature_type,
                        "mge_class": attrd.get("type") or feature_type,
                        "annotation": attrd.get("Note") or attrd.get("note") or "",
                        "score": score,
                        "source_program": source,
                        "raw_file": str(gff),
                        "raw_attributes": json.dumps(attrd, ensure_ascii=False),
                    }
                )
    return rows



def normalize_integronfinder(sample: str, sample_dir: Path) -> List[dict]:
    rows: List[dict] = []
    tool_dir = sample_dir / "integronfinder"
    if not tool_dir.exists():
        return rows

    for f in sorted(tool_dir.rglob("*.integrons")):
        try:
            df = read_table_auto(f)
        except Exception as exc:
            sys.stderr.write(f"[WARN] Failed to parse {f}: {exc}\n")
            continue

        # Expected schema documented by IntegronFinder API.
        # We still tolerate partial schemas.
        for _, r in df.iterrows():
            row = {k: str(v) for k, v in r.to_dict().items()}
            seqid = row.get("ID_replicon") or row.get("replicon") or row.get("seqid") or ""
            start = row.get("pos_beg") or row.get("start") or row.get("begin") or ""
            end = row.get("pos_end") or row.get("end") or row.get("stop") or ""
            strand = row.get("strand", "")
            feature_type = row.get("element") or row.get("type_elt") or "integron_feature"
            mge_name = row.get("annotation") or row.get("element") or row.get("ID_integron") or feature_type
            mge_class = row.get("type") or row.get("type_elt") or feature_type
            rows.append(
                {
                    "sample": sample,
                    "tool": "IntegronFinder",
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "feature_type": feature_type,
                    "mge_name": mge_name,
                    "mge_class": mge_class,
                    "annotation": row.get("annotation", ""),
                    "score": row.get("evalue", ""),
                    "source_program": "integron_finder",
                    "raw_file": str(f),
                    "raw_attributes": json.dumps(row, ensure_ascii=False),
                }
            )
    return rows



def normalize_isescan(sample: str, sample_dir: Path) -> List[dict]:
    rows: List[dict] = []
    tool_dir = sample_dir / "isescan"
    if not tool_dir.exists():
        return rows

    for f in sorted(tool_dir.glob("*.tsv")):
        try:
            df = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
        except Exception as exc:
            sys.stderr.write(f"[WARN] Failed to parse {f}: {exc}\n")
            continue

        for _, r in df.iterrows():
            row = {k: str(v) for k, v in r.to_dict().items()}
            rows.append(
                {
                    "sample": sample,
                    "tool": "ISEScan",
                    "seqid": row.get("seqID", ""),
                    "start": row.get("isBegin", ""),
                    "end": row.get("isEnd", ""),
                    "strand": row.get("strand", ""),
                    "feature_type": "IS_element",
                    "mge_name": row.get("family") or row.get("cluster") or "IS_element",
                    "mge_class": row.get("family", ""),
                    "annotation": row.get("type", ""),
                    "score": row.get("E-value4copy") or row.get("E-value") or row.get("score") or "",
                    "source_program": "ISEScan",
                    "raw_file": str(f),
                    "raw_attributes": json.dumps(row, ensure_ascii=False),
                }
            )
    return rows



def sample_dirs(root: Path) -> Iterable[Path]:
    for p in sorted(root.iterdir()):
        if p.is_dir() and (p / "tables").exists():
            yield p



def main() -> int:
    ap = argparse.ArgumentParser(description="Merge and normalize MGE annotation outputs.")
    ap.add_argument("--root", required=True, help="Pipeline output root directory")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    combined: List[dict] = []
    manifest_rows: List[dict] = []

    for sdir in sample_dirs(root):
        sample = sdir.name
        tables_dir = sdir / "tables"
        tables_dir.mkdir(parents=True, exist_ok=True)

        me_rows = normalize_mobileelementfinder(sample, sdir)
        if_rows = normalize_integronfinder(sample, sdir)
        is_rows = normalize_isescan(sample, sdir)

        for name, rows in [
            ("mobileelementfinder.normalized.tsv", me_rows),
            ("integronfinder.normalized.tsv", if_rows),
            ("isescan.normalized.tsv", is_rows),
        ]:
            out = tables_dir / name
            df = pd.DataFrame(rows)
            if df.empty:
                df = pd.DataFrame(
                    columns=[
                        "sample", "tool", "seqid", "start", "end", "strand",
                        "feature_type", "mge_name", "mge_class", "annotation",
                        "score", "source_program", "raw_file", "raw_attributes",
                    ]
                )
            df.to_csv(out, sep="\t", index=False)
            manifest_rows.append({"sample": sample, "table": str(out), "n_rows": len(df)})
            combined.extend(rows)

    combined_dir = root / "combined"
    combined_dir.mkdir(parents=True, exist_ok=True)

    combined_df = pd.DataFrame(combined)
    if combined_df.empty:
        combined_df = pd.DataFrame(
            columns=[
                "sample", "tool", "seqid", "start", "end", "strand",
                "feature_type", "mge_name", "mge_class", "annotation",
                "score", "source_program", "raw_file", "raw_attributes",
            ]
        )
    combined_tsv = combined_dir / "mge_annotations.all.tsv"
    combined_df.to_csv(combined_tsv, sep="\t", index=False)

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_df.to_csv(combined_dir / "normalized_tables.manifest.tsv", sep="\t", index=False)

    # BED-like export for quick visualization in IGV / pyGenomeTracks / ggplot.
    bed_cols = ["seqid", "start", "end", "mge_name", "tool", "strand", "sample", "feature_type", "mge_class"]
    bed_df = combined_df.copy()
    if not bed_df.empty:
        bed_df = bed_df.assign(
            start=lambda d: pd.to_numeric(d["start"], errors="coerce").fillna(1).astype(int) - 1,
            end=lambda d: pd.to_numeric(d["end"], errors="coerce").fillna(1).astype(int),
        )
        bed_df = bed_df[bed_cols]
    else:
        bed_df = pd.DataFrame(columns=bed_cols)
    bed_df.to_csv(combined_dir / "mge_annotations.all.bedlike.tsv", sep="\t", index=False, header=False)

    print(f"Wrote: {combined_tsv}")
    print(f"Wrote: {combined_dir / 'mge_annotations.all.bedlike.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
