# MGEs-annotation-pipeline-for-short-sequences
# MGE annotation pipeline for GI FASTA files

This bundle contains:

- `install_mge_annotation_tools.sh` — install MobileElementFinder, IntegronFinder, and ISEScan into separate conda environments.
- `run_mge_annotation_pipeline.sh` — run all three tools on one FASTA file or a directory of FASTA files.
- `merge_mge_annotations.py` — normalize the raw outputs into per-tool TSVs and one merged TSV.

## Quick start

```bash
bash install_mge_annotation_tools.sh
bash run_mge_annotation_pipeline.sh -i /path/to/gi_fastas -o mge_out -t 4
```

## Output layout

For each input FASTA `sampleX.fasta`, the pipeline creates:

```text
mge_out/
  sampleX/
    input/
    logs/
    mobileelementfinder/
    integronfinder/
    isescan/
    tables/
      mobileelementfinder.normalized.tsv
      integronfinder.normalized.tsv
      isescan.normalized.tsv
  combined/
    mge_annotations.all.tsv
    mge_annotations.all.bedlike.tsv
    normalized_tables.manifest.tsv
```

## Notes

- IntegronFinder is forced to run with `--linear`, which is appropriate for GI fragments/contigs.
- The normalization step preserves raw outputs and writes a standardized merged table with columns such as:
  `sample`, `tool`, `seqid`, `start`, `end`, `strand`, `feature_type`, `mge_name`, `mge_class`, `annotation`.
- The BED-like file uses 0-based starts so it can be fed more easily into genome-browser style tools.
