# bcftools +rescue

A bcftools plugin that predicts whether protein-truncating variants (PTVs) are rescued by other variants on the same haplotype or alternative splice site usage.

For each haplotype in a phased VCF, `+rescue` scores splice changes with SpliceAI models. It tags PTVs (stop gained, frameshift, canonical splice donor or acceptor) that another variant on the haplotype rescues, tags variants that induce a PTV through exon skipping or intron retention, and annotates whether the transcript is predicted to undergo nonsense-mediated decay (NMD).

For analysis supporting Hampstead JE, et al. Incompletely penetrant protein-truncating variants confer developmental disorder risk through partial reductions in gene expression. *medRxiv* (2026). https://doi.org/10.64898/2026.10.22.26366026, see the dedicated repository [here](https://github.com/jhampstead/ptv_penetrance).

## Requirements

- A **phased** VCF, annotated with Ensembl VEP (`INFO/CSQ`) or `bcftools csq` (`INFO/BCSQ`)
- The reference genome (FASTA)
- A GFF3 annotation matching the VCF, in the format used by `bcftools csq`
- A directory of SpliceAI models (see [Models](#models))
- The [cpliceai](https://github.com/jhampstead/cpliceai) library (and its dependencies, such as ONNX Runtime), to build the plugin

## Models

The model files we distribute are the SpliceAI models ([Jaganathan et al., 2019](https://doi.org/10.1016/j.cell.2018.12.015)), converted and optimised for ONNX Runtime. They remain under Illumina's original license for the SpliceAI models, [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/): they may be used for academic and non-commercial purposes only. Any other use requires a commercial license from Illumina (see the [SpliceAI repository](https://github.com/Illumina/SpliceAI)).

## Build

Build bcftools as usual (see [INSTALL](../INSTALL)). The plugin is built as `plugins/rescue.so`:

```bash
make
```

If pkg-config can't find cpliceai, point the build to its installation:

```bash
make CPLICEAI_ROOT=/path/to/cpliceai CPLICEAI_DEPENDENCY_DIRS=/path/to/onnxruntime/lib
```

Then tell bcftools where the plugins are:

```bash
export BCFTOOLS_PLUGINS=/path/to/bcftools/plugins
```

## Usage

```bash
bcftools +rescue in.vcf.gz -m models/ -g annot.gff3.gz -f ref.fa -o out.vcf.gz -O z
```

| Option | Description |
|---|---|
| `-f, --fasta-ref FILE` | Reference genome (required) |
| `-g, --gff-annot FILE` | GFF3 annotation (required) |
| `-m, --models DIR` | SpliceAI models (required) |
| `-o, --output FILE` | Output file [stdout] |
| `-O, --output-type u\|b\|v\|z` | Output format [v] |
| `--annotation vep\|bcftools` | Tool that annotated the VCF [vep] |
| `--csq-tag TAG` | INFO tag with the consequences [CSQ for vep, BCSQ for bcftools] |
| `--exclude-unphased` | Ignore heterozygous variants with unknown phase (by default they are placed on both haplotypes) |
| `--reference-scores FILE` | Precomputed reference scores, to save time (optional, see [Reference scores](#reference-scores-optional)) |
| `--raw-scores FILE` | Write the raw SpliceAI scores to a TSV file (optional, useful for debugging) |

Run `bcftools +rescue -h` for all options.

## Reference scores (optional)

Instead of predicting the reference sequence while it runs, `+rescue` can read precomputed reference scores. Make them once per genome and annotation with `cpliceai_reference`, which comes with [cpliceai](https://github.com/jhampstead/cpliceai):

```bash
cpliceai_reference models/ ref.fa regions.tsv reference.bin
bcftools +rescue in.vcf.gz -m models/ -g annot.gff3.gz -f ref.fa --reference-scores reference.bin -o out.vcf.gz -O z
```

Use the same models and FASTA as for `+rescue`; the plugin refuses scores made from a different FASTA. The regions file lists the genes to score (`NAME`, `CHROM`, `STRAND`, `TX_START`, `TX_END`). Any gene in your GFF3 that the regions file doesn't match exactly (same name, coordinates and strand) is scored while running instead.

## Output

The plugin adds these INFO tags:

| Tag | Content |
|---|---|
| `SpliceAI` | SpliceAI delta scores for the variant alone |
| `SpliceAI_Haplotype` | SpliceAI delta scores from the reference to the full haplotype |
| `CS` | Change in splicing between the haplotype with and without the variant |
| `RESCUE` | Whether the PTV is predicted to be rescued |
| `RESCUE_PROB` | Probability of rescue |
| `RESCUE_TYPE` | Type of rescue, e.g. `complete_rescue`, `frame_restoring_indel`, `exon_skipping`, `intron_retention` |
| `PTV` | Newly predicted PTV consequence |
| `PTV_PROB` | Probability of the induced PTV |
| `PTV_TYPE` | Type of induced PTV |
| `PREDICTION_FLAGS` | Context of the prediction, including NMD (`nmd_target`, `long_exon`, `distal`, `start_proximal`, `non-coding`, `nmd_uncertain`), based on the NMD regions of [Blakes et al. (2026)](https://doi.org/10.1038/s41467-026-69983-z) |

## Citations

This plugin builds on the following work:

- **SpliceAI models:** Jaganathan K, et al. Predicting splicing from primary sequence with deep learning. *Cell* 176, 535–548 (2019). https://doi.org/10.1016/j.cell.2018.12.015
- **NMD regions:** Blakes AJM, et al. Regional nonsense constraint offers biological and clinical insights into genetic disease. *Nature Communications* 17, 3152 (2026). https://doi.org/10.1038/s41467-026-69983-z
