# Golden fixture

`make_golden.py` writes a small genome, an annotation and phased VEP-annotated VCFs. Each record exercises one item of the
"Golden-file checklist" in `REFACTOR_PLAN.md` or one NMD class. Run both builds on the same files and compare.

```
python3 make_golden.py data          # standard library only, deterministic
```

| file | content |
|---|---|
| `golden.fa`, `golden.fa.fai` | contig `chr1` |
| `golden.gff3` | one protein-coding transcript per gene, forward and reverse copies of each template |
| `golden.vcf` | all scenarios that every build is expected to handle; INFO `AF`, `CSQ` (`Allele\|Consequence\|Feature_type\|Feature`), one sample |
| `golden_bcsq.vcf` | the records of `golden.vcf` annotated as bcftools csq does (INFO/BCSQ). Run with `--annotation bcftools`; everything but the CSQ/BCSQ tag has to equal the result for `golden.vcf` |
| `golden.regions.tsv` | the gene regions for `cpliceai_reference` (0-based half-open), to make the file of reference scores |
| `golden_edge.vcf` | a gene whose first donor is 49 bases from the start of the contig (windows running off the contig). Separate because an older build may crash on it |
| `golden_broken.vcf` | two valid records and then a truncated line: the read failure |
| `manifest.tsv` | per record: copy, template, strand, haplotype, position, alleles, checklist items, expectation |

## Running

Two runs, both with the precomputed reference scores: the VEP-annotated records and the same records annotated by bcftools csq.
Apart from the CSQ/BCSQ tag the two outputs have to be equal. From the repository root:

```
# The six positionals come first, options after them. REF.bin is the reference scores of the golden genes (golden.regions.tsv).
# Every build writes build/rescue (src/CMakeLists.txt), so copy the binary after building it to keep two builds apart.
M=../cpliceai_haplotype/models/onnx; D=tests/golden/data
build/rescue $D/golden.vcf $M $D/golden.gff3 $D/golden.fa vep.vcf vep_scores.tsv --af-tag AF --reference-scores REF.bin
build/rescue $D/golden_bcsq.vcf $M $D/golden.gff3 $D/golden.fa bcsq.vcf bcsq_scores.tsv --af-tag AF --annotation bcftools --reference-scores REF.bin

# old: its models are another format, and the option for unphased variants is --max-unphased 0 where the new one has --exclude-unphased
```

Run the new build on `golden_edge.vcf` and `golden_broken.vcf` as well. The broken file has to give an error and a non-zero exit
status, after the first records were written. The edge file has to finish.

The path of the binary and the model directory depend on the build; the arguments are `variants model_dir gff fasta output scores_output`.
`--af-tag AF` is needed because the default `INFO/AF` is not what the records carry.

## Reading the manifest

`expect` is one of:

- `same`: no difference between old and new, apart from the NMD flags.
- `diff`: an expected difference; `items` are the numbers in the checklist that say which and why.
- `start_proximal`, `long_exon`, `distal`: the one NMD flag the unrescued or induced PTV of this record gets, the
  escape region its new stop is in (on the most probable isoform with a PTV). Per transcript and phase.
- `none`: a PTV whose new stop is outside every escape region, or that has no new stop (N8, `frameshift_921_no_stop`: a
  frameshift with no stop in the new frame), gets no NMD flag. For N8 the 3' end of the transcript is 52 bases of UTR; whether the
  scan continues into it decides if a stop is found there. If a flag is given, see `first_stop_codon` for what is scanned.
- `nmd-model`: N9. Whether the lost donor retains the intron is up to the model. When it does, the first-intron record
  gets no flag and the last-intron record `distal` (the intron lies in the distal region).
- `model`: the outcome depends on what the model predicts for the designed sequence; compare old and new, no outcome is asserted.
- `info`: only has to run (the edge file).

The NMD scenarios use the template with a 600 base exon: coding indices 0-149 are start-proximal, 250-448 the long-exon part, 800-947
the distal region (last 50 of the penultimate exon and the last exon); the rest is the target region. The stop positions
are chosen on both sides of each border (99/147 start-proximal, 150/201 target, 300/447 long exon, 450/600 target, 798/801 on the
distal border, 900/930 distal).

## What was not verified

The script could not be run when it was written. It checks its own invariants (REF of every record is the generated genome, the
reference transcripts have no premature stop, planted features do not overlap, no two records of one haplotype overlap), so a mistake in
it fails when it runs. Whether the GFF reader accepts the files, and the phases of the CDS lines, is not checked: if
`gff_parse` stops, the message names the line.
