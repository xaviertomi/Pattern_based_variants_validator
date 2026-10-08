bcftools 1.23.1, samtools 1.23.1, MEME Suite 5.5.7 (STREME, FIMO, Tomtom, fasta-get-markov, fasta-shuffle-letters), Python 3.12.10, and R 4.5.1

# Standalone motif workflow

Snakemake workflow for discovering sequence motifs around selected PacBio variants, testing transfer to their corresponding Illumina samples, freezing validated motif families, and annotating a separate multi-sample Illumina VCF. The workflow is standalone; it does not read or write `clone_vcf_filter/`.

A motif match is sequence evidence, not variant truth. `promising` means a reference-context match passed the configured workflow criteria. No match is `unclassified`, never a false variant call.

## Local software environment

The biological tools **must be loaded locally** from your site's modules or equivalent host environment. To use this mode, edit `config.yaml` and set:

```yaml
execution:
  use_containers: false
```

Before launching, load the site-specific modules that provide the exact versions on the first line, and confirm they are available on `PATH` on both the controller and compute nodes. Module names vary by site; this workflow does not guess or load biological-tool modules. Its launcher only loads `bioinfo-ifb` and `snakemake/9.4.0` when Snakemake 9.4.0 is not already available.

Host requirements: Snakemake 9.4.0, `snakemake-executor-plugin-cluster-generic` 1.0.9, Python 3.12.10, Bash, and an authorized SLURM environment (`sbatch`, `sacct`, `squeue`). In local-module mode, the required biological executables are `bcftools`, `samtools`, `streme`, `fimo`, `tomtom`, `fasta-get-markov`, `fasta-shuffle-letters`, and `Rscript`, at the versions listed above. Run the workflow's tool preflight after loading modules; it exercises the actual local tools and fails on missing or incompatible versions.

`run`, `pilot`, and `rerun` submit Snakemake rule jobs through the configured cluster-generic SLURM profile; local module loading does not make these production targets run serially on the login node. The default profile account/partition (`dedicated-cpu@cirad-long` / `cpu-dedicated`) are site-specific and require authorization. Adjust `profile/config.yaml` and `profile/controller.yaml` only to match your site's approved values.

## Configure inputs

Edit `config.yaml` before running:

- Set `reference_fasta`, `pacbio_vcf`, `illumina_vcf`, and `other_illumina` to readable input paths.
- Set `output_dir` and `scratch_dir` to writable locations outside `clone_vcf_filter/`. Relative paths are resolved from this `Motif/` directory.
- Replace `pacbio_samples` and `illumina_corresponding` with exact, case-sensitive VCF header names. The correspondence mapping must cover every configured PacBio sample; names are not inferred from file names.
- Add `sample_metadata` such as `clone_id` where available. Without sufficient metadata, external-sample independence is reported as unverified rather than assumed.
- Set `folds.grouping` for the reference contig naming scheme. When enabled, provide a regex with a capture group and choose `error` or `exclude` for unmatched contigs. When disabled, each reference contig is its own group. Five folds are generated from the staged reference index; do not provide BED files or a pre-existing index.
- Set `pilot.pacbio_sample` and `pilot.rotation` for the optional diagnostic pilot. The pilot is not a substitute for all five production rotations.

`selection.min_supported_rotations` defaults to a proposed, configurable three-of-five rule. Review all scientific thresholds in `config.yaml`; they are workflow settings, not calibrated biological truth thresholds.

## Frozen-library VCF exports and consensus

`vcf_outputs` exports promising original VCF records for every corresponding and external sample. Exports keep original INFO, FORMAT, GT, multiallelic rows and duplicates; motif labels are prioritization evidence, not variant truth. Corresponding annotations use their own namespace; external annotations remain separate.

Corresponding final annotations are `application/corresponding/{genotypes.tsv.gz,alt_annotations.tsv.gz,all_matches.tsv.gz,annotation.complete.json}` with counts and positions under `reports/corresponding/`. Their rows carry `role=illumina`, `selection_participant=true`, `evidence_context=selection_reuse`; external rows remain `role=other_illumina`, `selection_participant=false`, with the existing identity-check context under `application/`.

Configure role-scoped consensus groups by exact VCF sample names:

```yaml
vcf_outputs:
  consensus_fraction: 0.75
  groups:
    corresponding_a:
      role: illumina
      samples: [IL-library-A-1, IL-library-A-2, IL-library-A-3]
    external_cohort:
      role: other_illumina
      samples: [external-sample-1, external-sample-2]
```

`consensus_fraction` defaults to `0.75`; a group retains an allele at `ceil(fraction × configured_member_count)` supporting members. Every configured member stays in that fixed denominator, including missing, absent, or conflicting calls. Group membership is explicit, role-specific, may overlap, and external names are checked against the source header. An empty `groups` mapping produces sample exports and a global completion marker without consensus groups.

Run `annotate_corresponding` to apply the frozen motifs to selected corresponding samples; `annotate_other` annotates external samples only. `reports` retains its existing selection-then-external-report flow. After the freeze, `export_promising` requests all sample exports, `consensus_groups` requests configured groups, and `vcf_deliverables` requests all VCF outputs. `all` includes VCF delivery.

`vcf_outputs/datasets.json` records sample roles, sources, annotations and freeze identity; `groups.json` records configured membership, denominator and threshold. Each `vcf_outputs/{dataset_id}/` contains `promising.vcf.gz`, its `.csi` and `export.complete.json`. Each `vcf_outputs/groups/{group_id}/` contains `consensus.vcf.gz`, its `.csi`, `support.tsv.gz`, `member_evidence.tsv.gz`, `excluded_alleles.tsv.gz` and `complete.json`; `vcf_outputs/outputs.complete.json` covers all sample and group deliverables. Group reports retain original source-level motif context and distinguish support, missingness, duplicate conflicts and excluded unsupported alleles.

The application-only targets verify the existing freeze and do not rerun discovery or selection. `annotate_other` needs only the external source; `annotate_corresponding` needs the corresponding source; VCF delivery needs both. `all` and `reports` retain selection-before-application behavior.


## Run

From this directory, after loading the local tool modules and editing `config.yaml`:

```bash
# Non-executing checks; these do not submit jobs.
bash snakemake.sh dry preflight --configfile config.yaml
bash snakemake.sh dry all --configfile config.yaml
bash snakemake.sh dag all --configfile config.yaml > workflow_dag.dot

# Real local-tool preflight, then a diagnostic pilot.
bash snakemake.sh run preflight --configfile config.yaml
bash snakemake.sh dry pilot --configfile config.yaml
bash snakemake.sh pilot --configfile config.yaml

# Production: full five-rotation selection/freeze and external application.
bash snakemake.sh run all --configfile config.yaml
```

`run all` freezes the complete selection before applying the frozen library and producing reports. To run stages separately, use `bash snakemake.sh run freeze_library --configfile config.yaml`, followed by `bash snakemake.sh run annotate_other --configfile config.yaml` or `bash snakemake.sh run reports --configfile config.yaml`. Application-only verifies the frozen library and does not need the original selection VCFs. `rerun` forces recomputation; use it only when that is intended. `unlock` is for a stale Snakemake lock after confirming that no workflow is still running.

## Main outputs

The configured `output_dir` contains the generated five-fold partition and provenance manifests, per-sample genotype/context tables, all five rotations of discovery and transfer evidence, family and exact-representative validation, and the immutable frozen library under `library/`. External annotations are written under `application/`; summary tables, PDFs, and a generated report README are written under `reports/`. The principal deliverables are `library/freeze.json`, `library/final_motifs.meme`, and `application/genotypes.tsv.gz`.

Checksum-based provenance controls resume and invalidation. Do not edit generated folds, frozen artifacts, or provenance sidecars by hand; a mismatch causes validation failure or recomputation.
