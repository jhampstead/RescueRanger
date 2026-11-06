# bcftools +rescue

A bcftools plugin that predicts whether protein-truncating variants (PTVs) are rescued by other variants on the same haplotype or alternative splice site usage.

For each haplotype in a phased VCF, `+rescue` scores splice changes with SpliceAI models. It tags PTVs (stop gained, frameshift, canonical splice donor or acceptor) that another variant on the haplotype rescues, tags variants that induce a PTV through exon skipping or intron retention, and annotates whether the transcript is predicted to undergo nonsense-mediated decay (NMD).

The predictions are made by `rescue`, the library of this repository, which can be run in two ways: as the bcftools plugin `bcftools +rescue`, or as a standalone library, without bcftools. See [Usage](#usage).

For analysis supporting Hampstead JE, et al. Incompletely penetrant protein-truncating variants confer developmental disorder risk through partial reductions in gene expression. *medRxiv* (2026). https://doi.org/10.64898/2026.10.22.26366026, see the dedicated repository [here](https://github.com/jhampstead/ptv_penetrance).

## Requirements

- A **phased** VCF, annotated with Ensembl VEP (`INFO/CSQ`) or `bcftools csq` (`INFO/BCSQ`). Records with more than one ALT allele are passed through untouched, so split them first, for example with `bcftools norm -m-`
- The reference genome (FASTA)
- A GFF3 annotation matching the VCF, in the format used by `bcftools csq`
- A directory of SpliceAI models (see [Models](#models))
- The `rescue` library of this repository, installed (see [Build](#build)). bcftools is only needed to run it as a plugin

## Models

The model files we distribute are the SpliceAI models ([Jaganathan et al., 2019](https://doi.org/10.1016/j.cell.2018.12.015)) converted and optimised for ONNX Runtime. They remain under Illumina's original license for the SpliceAI models [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/): they may be used for academic and non-commercial purposes only. Any other use requires a commercial license from Illumina (see the [SpliceAI repository](https://github.com/Illumina/SpliceAI)).

The ONNX Runtime-optimised models are in the `models/onnx` directory of [cpliceai](https://github.com/jhampstead/cpliceai). Pass that directory to
`-m` (plugin) or as the second argument to `rescue` (standalone library).

## Build

The `rescue` library of this repository is built first. It is all that is needed to run `rescue` as a standalone library:
installing it also installs the `rescue` program. To run it as a bcftools plugin, bcftools is then built against the
installed library.

### Dependencies

- [HTSlib](https://github.com/samtools/htslib)
- [cpliceai](https://github.com/jhampstead/cpliceai), with [ONNX Runtime](https://onnxruntime.ai), the library it runs the
  models with
- CMake 3.10 or later (the `cmake --install` below needs 3.15, `ctest --test-dir` 3.20), pkg-config and
  [Check](https://libcheck.github.io/check/) (for the unit tests) to build

Install HTSlib first (from your package manager or from source). Then install ONNX Runtime and cpliceai, here into `~/opt`:

```bash
curl -L -o onnxruntime.tgz https://github.com/microsoft/onnxruntime/releases/download/v1.28.0/onnxruntime-linux-x64-1.28.0.tgz
mkdir -p ~/opt/onnxruntime && tar -C ~/opt/onnxruntime --strip-components=1 -xzf onnxruntime.tgz

git clone https://github.com/jhampstead/cpliceai.git
cmake -S cpliceai -B cpliceai/build -DCMAKE_BUILD_TYPE=Release -DONNXRUNTIME_ROOT=$HOME/opt/onnxruntime -DCMAKE_INSTALL_PREFIX=$HOME/opt/cpliceai
cmake --build cpliceai/build
cmake --install cpliceai/build
```

Building cpliceai also needs [bats](https://github.com/bats-core/bats-core) and a recent C compiler (tested with GCC 13).
cpliceai installs no pkg-config file, so build the library with
`-DCPLICEAI_ROOT=$HOME/opt/cpliceai -DCPLICEAI_DEPENDENCY_DIRS=$HOME/opt/onnxruntime/lib` (see below). With Check from
Debian or Ubuntu packages, also add `-DCMAKE_MODULE_PATH=$PWD/cmake`.

bcftools itself is not needed to build the library.

### The library

```bash
cmake -S . -B build-lib -DCMAKE_BUILD_TYPE=Release
cmake --build build-lib
cmake --install build-lib            # or: sudo cmake --install build-lib
```

The prefix is `/usr/local` unless `-DCMAKE_INSTALL_PREFIX=...` is given. This installs

- `lib/librescue.so`
- `include/rescue/rescue.h` and `include/rescue/version.h` (generated, with the version of `CMakeLists.txt`)
- `lib/pkgconfig/rescue.pc`
- the program `bin/rescue`

To build against a cpliceai that pkg-config does not find, give its prefix and the directory of ONNX Runtime:

```bash
cmake -S . -B build-lib -DCPLICEAI_ROOT=/path/to/cpliceai -DCPLICEAI_DEPENDENCY_DIRS=/path/to/onnxruntime/lib
```

The directories are put in the RPATH of the library, so that the library finds them when it is loaded.

The unit tests are run with `ctest --test-dir build-lib` or `build-lib/tests/unit_test`.

### The plugin

The plugin is part of [bcftools](https://github.com/samtools/bcftools) but not of its release tarballs, so build bcftools
and HTSlib from GitHub (see bcftools' `INSTALL`). The plugin is built as
`plugins/rescue.so` when the library is found with pkg-config; otherwise it is left out with a notice:

```bash
git clone --recurse-submodules https://github.com/samtools/htslib.git
git clone https://github.com/samtools/bcftools.git
cd bcftools
make
```

If pkg-config can't find the library, point the build to its installation prefix:

```bash
make RESCUE_ROOT=/path/to/rescue
```

The directory of the library is put in the RPATH of the plugin. To tell bcftools where the plugin is, see
[Running as a bcftools plugin](#running-as-a-bcftools-plugin-bcftools-rescue).

## Usage

`bcftools +rescue` and the standalone `rescue` library run the same code: they take the same inputs, make the same predictions
and add the same INFO tags (see [Output](#output)). They differ in how the options are given and in the output formats.

### Running as a bcftools plugin (`bcftools +rescue`)

Tell bcftools where its plugins are, and check that `rescue` is among them:

```bash
export BCFTOOLS_PLUGINS=/path/to/bcftools/plugins
bcftools plugin -l
```

If `rescue` is not listed, `bcftools plugin -lvv` shows why it could not be loaded (for example, a library that was not found).
A plugin can also be given by its full path instead of its name: `bcftools +/path/to/bcftools/plugins/rescue.so ...`.

Then run it on a VCF:

```bash
bcftools +rescue in.vcf.gz -m models/ -g annot.gff3.gz -f ref.fa -o out.vcf.gz -O z
```

The options can come before or after the input file. Without an input file, the VCF is read from standard input, so the plugin
can be used in a pipeline, for example to run one chromosome at a time:

```bash
bcftools view -r chr17 in.vcf.gz -Ou | bcftools +rescue -m models/ -g annot.gff3.gz -f ref.fa -o chr17.vcf.gz -O z
```

| Option | Description |
|---|---|
| `-f, --fasta-ref FILE` | Reference genome (required) |
| `-g, --gff-annot FILE` | GFF3 annotation (required) |
| `-m, --models DIR` | SpliceAI models (required) |
| `-o, --output FILE` | Output file [stdout] |
| `-O, --output-type u\|b\|v\|z[0-9]` | Output format: `u`/`b` uncompressed/compressed BCF, `v`/`z` uncompressed/compressed VCF, `0-9` compression level [v] |
| `--annotation vep\|bcftools` | Tool that annotated the VCF [vep] |
| `--csq-tag TAG` | INFO tag with the consequences [CSQ for vep, BCSQ for bcftools] |
| `--af-tag TAG` | INFO tag with the allele frequency [INFO/AF] |
| `--af-threshold FLOAT` | Variants with an allele frequency at or below this trigger the splice analysis; more common variants are included in haplotypes but not tested individually [0.01] |
| `--gt-tag TAG` | FORMAT tag with the phased genotypes: `GT`, or a String field written like GT (`0\|1`) [FORMAT/GT] |
| `--exclude-unphased` | Ignore heterozygous variants with unknown phase (by default they are placed on both haplotypes) |
| `--prediction-radius INT` | Prediction window radius in bases; values other than 50 give scores not comparable to SpliceAI [50] |
| `--reference-scores FILE` | Precomputed reference scores, to save time (optional, see [Reference scores](#reference-scores-optional)) |
| `--raw-scores FILE` | Write the raw SpliceAI scores to a TSV file (optional, useful for debugging) |
| `--no-version` | Do not add the bcftools version and command line to the header |

Run `bcftools +rescue -h` for the help.

### Running as a standalone library

Without bcftools, `rescue` is run with the `rescue` program that is installed with the library, or called from your own code.

#### The `rescue` program

`cmake --install` puts the program in `<prefix>/bin/rescue` (`/usr/local/bin/rescue` by default). Before it is installed,
`cmake --build` writes it to `build/rescue` at the root of the repository, whatever the build directory is.

The program takes five arguments, in this order, followed by any options:

```
rescue <in.vcf> <models> <annot.gff3> <ref.fa> <out.vcf> [options]
```

| Argument | Description |
|---|---|
| `in.vcf` | Input VCF (or BCF), compressed or not |
| `models` | Directory with the SpliceAI models |
| `annot.gff3` | GFF3 annotation |
| `ref.fa` | Reference genome |
| `out.vcf` | Annotated VCF to write, or `-` for standard output |

The plugin example above, run standalone:

```bash
rescue in.vcf.gz models/ annot.gff3.gz ref.fa - | bgzip > out.vcf.gz
```

How it differs from the plugin:

- **Options go after the five arguments.** Anything before them is read as one of the arguments. An option the program doesn't
  know is skipped with a warning (`Ignoring invalid argument`) rather than stopping the run.
- **The output is always uncompressed VCF**, whatever the file name ends with. For a compressed VCF, write to `-` and pipe it to
  `bgzip`, as above, or run `bgzip` on the file afterwards.

| Option | Description |
|---|---|
| `--annotation vep\|bcftools` | Tool that annotated the VCF [vep] |
| `--csq-tag TAG` | INFO tag with the consequences [CSQ for vep, BCSQ for bcftools] |
| `--af-tag TAG` | INFO tag with the allele frequency [INFO/AF] |
| `--af-threshold FLOAT` | Variants with an allele frequency at or below this trigger the splice analysis; more common variants only take part in the haplotypes [0.01] |
| `--gt-tag TAG` | FORMAT tag with the phased genotypes: `GT`, or a String field written like GT (`0\|1`) [FORMAT/GT] |
| `--exclude-unphased` | Ignore heterozygous variants with unknown phase (by default they are placed on both haplotypes) |
| `--prediction-radius INT` | Prediction window radius in bases; values other than 50 give scores not comparable to SpliceAI [50] |
| `--reference-scores FILE` | Precomputed reference scores, to save time (optional, see [Reference scores](#reference-scores-optional)) |
| `--raw-scores FILE` | Write the raw SpliceAI scores to a TSV file (optional, useful for debugging) |

Run `rescue -h` for the help.

#### From your own code

The `rescue` library can be called from C (or C++). Its interface is `rescue/rescue.h`, the same one the plugin
and the `rescue` program are built on. The steps are those of the program: set the options, read the VCF header with HTSlib,
`rr_init` (which loads the models, the annotation and the reference, and adds the output INFO tags to the header), write the
header, then `rr_annotate` to read, annotate and write every record.

```c
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <rescue/rescue.h>

int main(void)
{
    RescueOptions opts;
    rr_opts_init(&opts);                    // the defaults; the required fields are NULL
    opts.model_dir   = "models/";
    opts.gff_fname   = "annot.gff3.gz";
    opts.fasta_fname = "ref.fa";
    // opts.annotation = RR_ANNOTATION_BCFTOOLS;  // for INFO/BCSQ from bcftools csq

    int ret = EXIT_FAILURE;
    htsFile *in = NULL, *out = NULL;
    bcf_hdr_t *hdr = NULL;
    Rescue *rr = NULL;

    if (!(in = hts_open("in.vcf.gz", "r")) || !(hdr = bcf_hdr_read(in))) goto done;
    if (!(rr = rr_init(&opts, hdr))) goto done;     // adds the output tags to hdr

    if (!(out = hts_open("out.vcf.gz", "wz")) || bcf_hdr_write(out, hdr) != 0) goto done;
    if (rr_annotate(rr, in, out) == 0) ret = EXIT_SUCCESS;

done:
    rr_destroy(rr);
    if (out && hts_close(out) != 0) ret = EXIT_FAILURE;
    if (hdr) bcf_hdr_destroy(hdr);
    if (in) hts_close(in);
    return ret;
}
```

Compile it with pkg-config. `rescue.pc` does not list HTSlib, so ask for it as well; the RPATH lets the program find the
library at run time when it is installed outside the standard directories:

```bash
cc -o annotate annotate.c $(pkg-config --cflags --libs rescue htslib) \
    -Wl,-rpath,$(pkg-config --variable=libdir rescue)
```

If ONNX Runtime is not in a standard library directory, the linker stops with `libonnxruntime.so.1, needed by libcpliceai.so.1,
not found`. Then also give it the directory of ONNX Runtime (the one given as `CPLICEAI_DEPENDENCY_DIRS`):
`-Wl,-rpath-link,/path/to/onnxruntime/lib`.

The header documents the rest of the interface: `rr_set_scores_output` to write the raw scores, `rr_set_log_level` for what is
logged to standard error, and `rr_version` for the version of the linked library (`RR_VERSION` is that of the header). Note that
the GFF reader, like that of bcftools, ends the program when the GFF cannot be parsed.

## Reference scores (optional)

Instead of predicting the reference sequence while it runs, `rescue` can read precomputed reference scores, both as the plugin and standalone. Make them once per genome and annotation with `cpliceai_reference`, which comes with [cpliceai](https://github.com/jhampstead/cpliceai):

```bash
cpliceai_reference models/ ref.fa regions.tsv reference.bin

# as a bcftools plugin
bcftools +rescue in.vcf.gz -m models/ -g annot.gff3.gz -f ref.fa --reference-scores reference.bin -o out.vcf.gz -O z

# standalone
rescue in.vcf.gz models/ annot.gff3.gz ref.fa out.vcf --reference-scores reference.bin
```

From your own code, set `opts.reference_scores` before `rr_init`.

Use the same models and FASTA as for the predictions; scores made from a different FASTA are refused. The regions file lists the genes to score (`NAME`, `CHROM`, `STRAND`, `TX_START`, `TX_END`). Any gene in your GFF3 that the regions file doesn't match exactly (same name, coordinates and strand) is scored while running instead.

## Output

Both ways of running add these INFO tags:

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

The raw SpliceAI scores, written with `--raw-scores`, are a TSV with the columns `CHR`, `POS`, `REF`, `ALT`, `GENE`,
`TRANSCRIPT`, `PHASE`, `OFFSET`, `SOURCE`, `TYPE`, `DONOR`, `ACCEPTOR` and `NEITHER`.

## Citations

`rescue` builds on the following work:

- **SpliceAI models:** Jaganathan K, et al. Predicting splicing from primary sequence with deep learning. *Cell* 176, 535–548 (2019). https://doi.org/10.1016/j.cell.2018.12.015
- **NMD regions:** Blakes AJM, et al. Regional nonsense constraint offers biological and clinical insights into genetic disease. *Nature Communications* 17, 3152 (2026). https://doi.org/10.1038/s41467-026-69983-z
