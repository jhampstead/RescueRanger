#!/usr/bin/env python3
"""
Writes the golden-file fixture of rescue_ranger: a small genome, its annotation and phased, VEP-annotated VCFs whose records exercise
the changes listed in REFACTOR_PLAN.md ("Golden-file checklist") and the NMD flags, so that two builds (or the code before and after
a change) can be run on the same input and their output compared.

    python3 make_golden.py [output directory]        (default: ./data next to this script)

Nothing but the standard library is used. The output is deterministic: it is the same on every run.

Files
    golden.fa, golden.fa.fai   contig chr1
    golden.gff3                the annotation: protein-coding genes with one transcript each, as the vendored GFF reader expects
    golden.vcf                 the records of the scenarios that every build is expected to handle
    golden_edge.vcf            records at the start of the contig (splice sites within 50 bases of position 0, windows that run
                               off the contig); kept apart because an older build may crash on them
    golden_broken.vcf          valid records followed by a truncated line (the read failure, item 19)
    manifest.tsv               one row per record: scenario, checklist items, what is expected

Design
    Every scenario has a gene copy of its own, with its own random sequence, so that the haplotype of one scenario does not
    touch another. A copy has two haplotypes, so it hosts two scenarios (phase 0 is genotype 1|0, phase 1 is 0|1). A scenario with
    an unphased record (0/1) is on both haplotypes of its copy. Genes come in two templates (the standard one and the NMD one) and on
    both strands; the same scenario is run on both. Variants are described on the transcript (offsets in the direction of the
    transcript) and turned into VCF records on the + strand here, with REF read from the generated genome.

    Splice sites are real motifs (GTAAGT at the donor, a pyrimidine tract and AG at the acceptor), and some copies have a cryptic
    donor or acceptor planted in an intron, so that a lost canonical site has somewhere to go. What a model makes of the sequence is
    not known here: a scenario that depends on it is marked "model" in the manifest.

Run the program with --af-tag AF (the records carry INFO/AF), see README.md.
"""
import os
import sys

CONTIG = "chr1"
FLANK = 1500
RARE = 0.001      # At most the default threshold of --af-threshold (0.01): the record triggers the splice analysis
COMMON = 0.3      # A record that only takes part in a haplotype

COMP = str.maketrans("ACGT", "TGCA")
ROT = {"A": "C", "C": "G", "G": "T", "T": "A"}
STOP_CODONS = {"TAA", "TAG", "TGA"}

DONOR_MOTIF = "GTAAGT"
ACCEPTOR_MOTIF = "TTCCTTTCTTTCTTTCTTAG"          # 20 bases
DECOY_DONOR = "CAGGTAAGT"                          # The GT is 3 bases in
DECOY_ACCEPTOR = "TTCCTTTCTTTCTTTCTTAG"


def comp(s):
    return s.translate(COMP)


def revcomp(s):
    return comp(s)[::-1]


class Rng:
    """A generator that gives the same numbers on every Python version."""

    def __init__(self, seed):
        self.state = seed & 0xFFFFFFFF

    def next(self):
        self.state = (self.state * 1664525 + 1013904223) & 0xFFFFFFFF
        return self.state >> 8

    def base(self):
        return "ACGT"[self.next() & 3]

    def sense_codon(self):
        while True:
            codon = self.base() + self.base() + self.base()
            if codon not in STOP_CODONS:
                return codon


RNG = Rng(20240901)


class Template:
    """
    The structure of a gene in the direction of the transcript. exons is a list of (length, cds_offset, cds_length), where the
    offset is from the 5' end of the exon, and introns the lengths between them. constitutive has a flag for the coding region of
    every exon. Offsets t of a gene are bases from the 5' end of its transcript, whatever the strand.
    """

    def __init__(self, exons, introns, constitutive):
        self.exons = exons
        self.introns = introns
        self.constitutive = constitutive
        assert len(introns) == len(exons) - 1 and len(constitutive) == len(exons)

        self.exon_t, self.intron_t, self.cds_t = [], [], []
        t = 0
        for k, (length, cds_offset, cds_length) in enumerate(exons):
            self.exon_t.append((t, t + length))
            self.cds_t.extend(range(t + cds_offset, t + cds_offset + cds_length))
            t += length
            if k < len(introns):
                assert introns[k] >= 40, "an intron has to hold both motifs and room for a cryptic site"
                self.intron_t.append((t, t + introns[k]))
                t += introns[k]
        self.length = t
        assert len(self.cds_t) % 3 == 0, "the coding sequence has to be whole codons"


# The standard template: 3 exons, 273 coding bases. Exon 0 has a 5' UTR of 60 and the CDS of 90, exon 1 is all CDS (100), exon 2
# has the CDS of 83 and a 3' UTR of 67. Introns of 200. The coding regions of the first and the last exon are not constitutive.
STD = Template([(150, 60, 90), (100, 0, 100), (150, 0, 83)], [200, 200], [0, 1, 0])

# The NMD template: 948 coding bases. Exon 1 has 600 coding bases, so it is a long exon (the first 199 bases are more than 400 before
# its last base), the first 150 are start-proximal, and the last 50 of exon 1 with the whole of exon 2 are distal. In coding indices:
# 0-149 start-proximal, 250-448 long exon, 800-947 distal; the rest is the NMD target region.
NMD = Template([(300, 50, 250), (600, 0, 600), (150, 0, 98)], [200, 200], [1, 1, 1])

# A gene at the start of the contig: the donor of the first exon is within 50 bases of position 0
EDGE = Template([(45, 9, 36), (156, 0, 156)], [120], [1, 1])

# A gene of one exon, for the pair of genes that share a base
SINGLE = Template([(200, 10, 180)], [], [1])


def nmd_region(coding_index):
    """The NMD flag the first base of a stop codon at this coding index of the NMD template gives: its escape region, or none."""
    if coding_index < 150:
        return "start_proximal"
    if 250 <= coding_index < 449:
        return "long_exon"
    if coding_index >= 800:
        return "distal"
    return "none"


class Genome:
    def __init__(self, length):
        self.b = [RNG.base() for _ in range(length)]

    def get(self, i):
        assert 0 <= i < len(self.b), "position %d outside the contig" % i
        return self.b[i]

    def seq(self, i, n):
        assert 0 <= i and i + n <= len(self.b), "range %d+%d outside the contig" % (i, n)
        return "".join(self.b[i:i + n])

    def put(self, i, s):
        assert 0 <= i and i + len(s) <= len(self.b)
        self.b[i:i + len(s)] = list(s)


class Rec:
    """A VCF record on the + strand: 0-based position, REF and ALT."""

    def __init__(self, pos, ref, alt):
        self.pos, self.ref, self.alt = pos, ref, alt


class Gene:
    def __init__(self, genome, tpl, strand, name, start, number):
        self.genome, self.tpl, self.fwd, self.name, self.start = genome, tpl, strand == "+", name, start
        self.L = tpl.length
        self.gid = "ENSG%05d" % number
        self.tid = "ENST%05d" % number
        self.claims, self.edits = [], []

        n = len(tpl.cds_t)
        self.cds = list("ATG" + "".join(RNG.sense_codon() for _ in range(n // 3 - 2)) + "TAA")
        self.seq = [RNG.base() for _ in range(self.L)]

        for k in range(len(tpl.introns)):
            s, e = tpl.intron_t[k]
            self.claim(s, s + len(DONOR_MOTIF), "donor %d" % k)
            self.claim(e - len(ACCEPTOR_MOTIF), e, "acceptor %d" % k)

    # --- planting, before finish() ---
    def claim(self, t0, t1, tag):
        for (s, e, other) in self.claims:
            assert t1 <= s or t0 >= e, "%s: %s [%d, %d) overlaps %s [%d, %d)" % (self.name, tag, t0, t1, other, s, e)
        self.claims.append((t0, t1, tag))

    def ct(self, coding_index):
        return self.tpl.cds_t[coding_index]

    def plant_cds(self, coding_index, bases, tag="plant"):
        ts = [self.tpl.cds_t[coding_index + i] for i in range(len(bases))]
        assert all(ts[i + 1] == ts[i] + 1 for i in range(len(ts) - 1)), "%s: planted coding bases %d+ span an intron" % (self.name, coding_index)
        self.claim(ts[0], ts[-1] + 1, "%s at %d" % (tag, coding_index))
        for i, base in enumerate(bases):
            self.cds[coding_index + i] = base

    def plant_intron(self, k, offset, bases, tag="plant"):
        s, e = self.tpl.intron_t[k]
        t0 = s + offset
        assert t0 >= s and t0 + len(bases) <= e
        self.claim(t0, t0 + len(bases), "%s in intron %d" % (tag, k))
        self.edits.append((t0, bases))

    def decoy_donor(self, k, offset):
        """A cryptic donor whose GT is at offset + 3 of intron k: a lost canonical donor keeps offset + 3 bases of the intron."""
        self.plant_intron(k, offset, DECOY_DONOR, "cryptic donor")

    def decoy_acceptor(self, k, last_offset):
        """A cryptic acceptor whose AG ends at the offset last_offset of intron k."""
        self.plant_intron(k, last_offset - len(DECOY_ACCEPTOR) + 1, DECOY_ACCEPTOR, "cryptic acceptor")

    def finish(self):
        for i, t in enumerate(self.tpl.cds_t):
            self.seq[t] = self.cds[i]

        codons = ["".join(self.cds[i:i + 3]) for i in range(0, len(self.cds), 3)]
        assert codons[-1] in STOP_CODONS
        for i, codon in enumerate(codons[:-1]):
            assert codon not in STOP_CODONS, "%s: a stop codon at the codon %d of the reference" % (self.name, i)

        for k in range(len(self.tpl.introns)):
            s, e = self.tpl.intron_t[k]
            self.seq[s:s + len(DONOR_MOTIF)] = list(DONOR_MOTIF)
            self.seq[e - len(ACCEPTOR_MOTIF):e] = list(ACCEPTOR_MOTIF)
        for t0, bases in self.edits:
            self.seq[t0:t0 + len(bases)] = list(bases)

        s = "".join(self.seq)
        self.genome.put(self.start, s if self.fwd else revcomp(s))

    # --- positions ---
    def rel(self, t):
        """The offset from the start of the gene on the + strand of the base at t of the transcript."""
        return t if self.fwd else self.L - 1 - t

    def pos(self, t):
        return self.start + self.rel(t)

    def tbase(self, t):
        """The base at t in the direction of the transcript."""
        b = self.genome.get(self.pos(t))
        return b if self.fwd else comp(b)

    def intron_start(self, k):
        return self.tpl.intron_t[k][0]

    def exon_end(self, k):
        return self.tpl.exon_t[k][1] - 1

    # --- variants, in the direction of the transcript ---
    def snv(self, t, alt=None):
        p = self.pos(t)
        ref = self.genome.get(p)
        alt_t = alt if alt is not None else ROT[self.tbase(t)]
        alt_g = alt_t if self.fwd else comp(alt_t)
        assert alt_g != ref, "%s: the SNV at %d does not change the base" % (self.name, t)
        return Rec(p, ref, alt_g)

    def dele(self, t, k):
        """Deletes the k bases from t on in the direction of the transcript."""
        lo = min(self.pos(t), self.pos(t + k - 1))
        anchor = lo - 1
        ref = self.genome.seq(anchor, k + 1)
        return Rec(anchor, ref, ref[0])

    def ins(self, a, s):
        """Inserts s (in the direction of the transcript) between the bases a and a + 1."""
        if self.fwd:
            anchor, added = self.pos(a), s
        else:
            anchor, added = self.pos(a + 1), revcomp(s)
        base = self.genome.get(anchor)
        return Rec(anchor, base, base + added)

    def mnv(self, t, k):
        """Replaces k bases by their complements."""
        lo = min(self.pos(t), self.pos(t + k - 1))
        ref = self.genome.seq(lo, k)
        return Rec(lo, ref, comp(ref))

    def cplx(self, t, k, new_length):
        """Replaces k bases by new_length bases, none of them the first base of REF: a complex indel without an anchor."""
        assert new_length != k
        lo = min(self.pos(t), self.pos(t + k - 1))
        ref = self.genome.seq(lo, k)
        return Rec(lo, ref, comp(ref)[:new_length])

    # --- variants, on the + strand (offsets from the start of the gene) ---
    def snv_g(self, rel):
        p = self.start + rel
        ref = self.genome.get(p)
        return Rec(p, ref, ROT[ref])

    def dele_g(self, rel, k):
        """Deletes the k bases after the base rel."""
        ref = self.genome.seq(self.start + rel, k + 1)
        return Rec(self.start + rel, ref, ref[0])

    def ins_g(self, rel, s):
        base = self.genome.get(self.start + rel)
        return Rec(self.start + rel, base, base + s)


# ------------------------------------------------------------------------------------------------------------------------------
# Scenarios
# ------------------------------------------------------------------------------------------------------------------------------

def item(rec, csq=None, af=RARE, gt=None):
    """A record of a scenario: the VEP consequences for the transcript, INFO/AF (None for none), and a genotype for the unphased."""
    return (rec, csq, af, gt)


class Scn:
    """
    phase: 0 or 1, the haplotype of the records. items: the checklist items it touches. expect: "same" (no difference between
    builds expected), "diff" (an expected difference, see items), "start_proximal", "long_exon" or "distal" (the NMD flag the
    record gets), "none" (a PTV that gets no NMD flag), "nmd-model"
    (an NMD flag if the model predicts the splice event) or "info".
    """

    def __init__(self, label, phase, items, expect, desc, build):
        self.label, self.phase, self.items, self.expect, self.desc, self.build = label, phase, items, expect, desc, build


COPIES = []        # (name, template, strand, setup, scenarios)


def copy(name, tpl, strand, setup, scenarios):
    COPIES.append((name, tpl, strand, setup, scenarios))


def both(name, tpl, setup, scenarios):
    """The scenarios on a forward and on a reverse gene."""
    copy(name + "F", tpl, "+", setup, scenarios)
    copy(name + "R", tpl, "-", setup, scenarios)


def no_setup(g):
    pass


def donor_decoy_60(g):
    g.decoy_donor(0, 57)       # A lost donor keeps 60 bases of the first intron


def donor_decoy_61(g):
    g.decoy_donor(0, 58)       # 61 bases: with a deleted base in an exon the frame is restored


def acceptor_decoy_near(g):
    g.decoy_acceptor(0, 159)   # 40 bases before the canonical acceptor


# --- baseline: simple records that every build is expected to give the same ---
both("A1", STD, donor_decoy_60, [
    Scn("donor_snv", 0, "", "same", "SNV at the first base of the donor, a cryptic donor 60 bases in",
        lambda g: [item(g.snv(g.intron_start(0)), "splice_donor_variant")]),
    Scn("cds_snv", 1, "", "same", "missense SNV in the coding sequence of exon 1, rare",
        lambda g: [item(g.snv(g.ct(140)), "missense_variant")]),
])
both("A2", STD, acceptor_decoy_near, [
    Scn("acceptor_snv", 0, "", "same", "SNV at the last base of the acceptor, a cryptic acceptor 40 bases before it",
        lambda g: [item(g.snv(g.tpl.intron_t[0][1] - 1), "splice_acceptor_variant")]),
    Scn("utr_snv", 1, "", "same", "SNV in the 5' UTR, rare",
        lambda g: [item(g.snv(30), "5_prime_UTR_variant")]),
])
both("A3", STD, lambda g: g.plant_cds(93, "CAA", "stop codon"), [
    Scn("vep_stop", 0, "", "start_proximal", "stop_gained SNV (CAA to TAA) at the coding index 93: start-proximal",
        lambda g: [item(g.snv(g.ct(93), "T"), "stop_gained", None)]),
    Scn("intron_snv", 1, "", "same", "SNV deep in the first intron, rare",
        lambda g: [item(g.snv(250), "intron_variant")]),
])

# --- item 1: isoform positions with an indel in or next to an intron ---
both("B1", STD, donor_decoy_60, [
    Scn("donor_and_intron_del", 0, "1", "diff", "lost donor and a 3 base deletion in the retained part of the intron, one haplotype",
        lambda g: [item(g.snv(g.intron_start(0)), "splice_donor_variant"), item(g.dele(180, 3), "intron_variant", COMMON)]),
    Scn("exon_end_frameshift", 1, "1", "diff", "deletion of the 3rd last base of exon 1: a frameshift next to the intron",
        lambda g: [item(g.dele(447, 1), "frameshift_variant", None)]),
])

# --- item 2: partial indels across a border or the edge of the gene ---
both("B2", STD, no_setup, [
    Scn("del_over_cds_start", 0, "2,10", "diff", "deletion of 6 bases across the start of the coding sequence (2 in the UTR)",
        lambda g: [item(g.dele(58, 6), "start_lost&5_prime_UTR_variant")]),
    Scn("del_over_cds_end", 1, "2,10", "diff", "deletion of 8 bases across the end of the coding sequence (the stop codon and the UTR)",
        lambda g: [item(g.dele(731, 8), "stop_lost&3_prime_UTR_variant")]),
])
both("B3", STD, no_setup, [
    Scn("del_over_donor", 0, "2", "diff", "deletion of 7 bases across the exon end and the donor",
        lambda g: [item(g.dele(447, 7), "splice_donor_variant&coding_sequence_variant")]),
    Scn("del_over_acceptor", 1, "2", "diff", "deletion of 7 bases across the acceptor and the exon start",
        lambda g: [item(g.dele(346, 7), "splice_acceptor_variant&coding_sequence_variant")]),
])
both("B4", STD, no_setup, [
    Scn("del_over_gene_end", 0, "2,11", "diff", "deletion that begins inside the gene and runs 3 bases past its last base",
        lambda g: [item(g.dele_g(g.L - 3, 5), "3_prime_UTR_variant&downstream_gene_variant")]),
    Scn("del_over_gene_start", 1, "2,11", "diff", "deletion that begins 2 bases before the gene and runs 2 bases into it",
        lambda g: [item(g.dele_g(-3, 4), "5_prime_UTR_variant&upstream_gene_variant")]),
])

# --- item 3: regions in haplotype coordinates: MNVs, complex indels and insertions at borders ---
both("B5", STD, no_setup, [
    Scn("mnv_over_donor", 0, "3", "diff", "MNV of 4 bases across the exon end and the donor, and a frameshift deletion before it",
        lambda g: [item(g.dele(380, 1), "frameshift_variant", None), item(g.mnv(448, 4), "splice_donor_variant")]),
    Scn("complex_over_donor", 1, "3", "diff", "complex indel (5 bases to 2) across the exon end and the donor",
        lambda g: [item(g.cplx(447, 5, 2), "splice_donor_variant&coding_sequence_variant")]),
])
both("B6", STD, no_setup, [
    Scn("ins_at_acceptor_border", 0, "3", "diff", "insertion of 2 bases between the last intron base and the first exon base, and a frameshift before it",
        lambda g: [item(g.dele(120, 1), "frameshift_variant", None), item(g.ins(349, "AC"), "splice_acceptor_variant")]),
    Scn("ins_at_cds_utr_border", 1, "3", "diff", "insertion of 2 bases after the last coding base, before the 3' UTR, and a frameshift before it",
        lambda g: [item(g.dele(670, 1), "frameshift_variant", None), item(g.ins(732, "AC"), "3_prime_UTR_variant")]),
])

# --- item 4: an insertion on the last base of a UTR region (genomic coordinates) ---
both("B7", STD, no_setup, [
    Scn("ins_utr_edge_next_to_cds", 0, "4,5", "diff", "insertion of 2 bases after the UTR base next to the coding sequence (the lower UTR region of the + strand)",
        lambda g: [item(g.ins_g(66 if not g.fwd else 59, "AC"), "3_prime_UTR_variant" if not g.fwd else "5_prime_UTR_variant")]),
    Scn("ins_after_last_gene_base", 1, "4,11", "diff", "insertion of 2 bases after the last base of the gene on the + strand",
        lambda g: [item(g.ins_g(g.L - 1, "AC"), "downstream_gene_variant")]),
])

# --- item 5: variants next to a UTR ---
both("B8", STD, no_setup, [
    Scn("utr_and_cds_snvs", 0, "5", "diff", "SNVs on the last base of the 5' UTR and on the first coding base",
        lambda g: [item(g.snv(59), "5_prime_UTR_variant"), item(g.snv(60), "start_lost")]),
    Scn("del_anchored_in_utr", 1, "5", "diff", "deletion of the first coding base, anchored in the 5' UTR",
        lambda g: [item(g.dele(60, 1), "start_lost", None)]),
])
both("B9", STD, no_setup, [
    Scn("ins_in_3utr", 0, "5", "diff", "insertion of 2 bases in the first base of the 3' UTR",
        lambda g: [item(g.ins(733, "AC"), "3_prime_UTR_variant")]),
    Scn("snv_after_deletion", 1, "5", "diff", "an in-frame deletion of 30 coding bases and an SNV after it",
        lambda g: [item(g.dele(70, 30), "inframe_deletion", None), item(g.snv(140), "missense_variant")]),
])

# --- item 6: NON_CONSTITUTIVE and INFRAME flags at the edge of the coding sequence (exon 0 is not constitutive) ---
both("B10", STD, no_setup, [
    Scn("ins_at_exon_end", 0, "6", "diff", "insertion of 2 bases after the last base of the exon: the border with the intron",
        lambda g: [item(g.ins(149, "AC"), "frameshift_variant&splice_region_variant")]),
    Scn("del_into_intron", 1, "6", "diff", "deletion of the first intron base, anchored on the last exon base",
        lambda g: [item(g.dele(150, 1), "splice_donor_variant")]),
])
both("B11", STD, no_setup, [
    Scn("cplx_exon_to_intron", 0, "6", "diff", "complex indel (3 bases to 1) from the exon end into the intron",
        lambda g: [item(g.cplx(148, 3, 1), "splice_donor_variant&coding_sequence_variant")]),
    Scn("del_over_acceptor_anchor", 1, "6", "diff", "deletion of the last intron base and the first exon base, anchored in the intron",
        lambda g: [item(g.dele(349, 2), "splice_acceptor_variant&coding_sequence_variant")]),
])

# --- item 7: the reading frame trigger counts a deletion at its first deleted base ---
both("B12", STD, no_setup, [
    Scn("snv_then_del_into_intron", 0, "7", "diff", "SNV at the 3rd last exon base and a deletion of the first intron base (anchored in the exon)",
        lambda g: [item(g.snv(448), "missense_variant", COMMON), item(g.dele(450, 1), "splice_donor_variant")]),
    Scn("del_at_cds_start_then_snv", 1, "7", "diff", "deletion of the first coding base (anchored in the UTR) and an SNV 2 bases after it",
        lambda g: [item(g.dele(60, 1), "start_lost"), item(g.snv(62), "missense_variant", COMMON)]),
])

# --- item 8: frame-restoring pairing over the unique shifts ---
both("B13", STD, donor_decoy_61, [
    Scn("two_donor_snvs_and_frameshift", 0, "8", "diff", "two SNVs at the donor (the same splice site changes twice) and a frameshift deletion; a cryptic donor 61 bases in",
        lambda g: [item(g.snv(150), "splice_donor_variant"), item(g.snv(151), "splice_donor_variant"), item(g.dele(400, 1), "frameshift_variant", None)]),
    Scn("one_donor_snv_and_frameshift", 1, "8", "diff", "the same with one donor SNV",
        lambda g: [item(g.snv(150), "splice_donor_variant"), item(g.dele(400, 1), "frameshift_variant", None)]),
])

# --- item 9: splice windows from the reference bounds; cryptic acceptors 25 and 50 bases after the donor ---
def acceptors_near_donor(g):
    g.decoy_acceptor(0, 30)
    g.decoy_acceptor(0, 55)

both("B14", STD, acceptors_near_donor, [
    Scn("del_early_in_intron", 0, "9", "model", "acceptor lost, a deletion of 5 bases early in the intron; cryptic acceptors near the donor",
        lambda g: [item(g.dele(165, 5), "intron_variant", COMMON), item(g.snv(g.tpl.intron_t[0][1] - 1), "splice_acceptor_variant")]),
    Scn("del_late_in_intron", 1, "9", "model", "acceptor lost, a deletion of 4 bases near the acceptor",
        lambda g: [item(g.dele(330, 4), "intron_variant", COMMON), item(g.snv(g.tpl.intron_t[0][1] - 2), "splice_acceptor_variant")]),
])

# --- item 10: the stop scan at the ends of the coding sequence ---
both("B15", STD, no_setup, [
    Scn("del_first_coding_bases", 0, "10", "diff", "deletion of the first 6 coding bases",
        lambda g: [item(g.dele(60, 6), "start_lost")]),
    Scn("del_last_coding_bases", 1, "10", "diff", "deletion of the last 6 coding bases (the stop codon and the codon before)",
        lambda g: [item(g.dele(727, 6), "stop_lost")]),
])
both("B16", STD, no_setup, [
    Scn("long_ins_near_end", 0, "10", "diff", "insertion of 31 bases 12 bases before the end of the coding sequence",
        lambda g: [item(g.ins(720, "ACGT" * 7 + "ACG"), "frameshift_variant")]),
    Scn("long_ins_near_start", 1, "10", "diff", "insertion of 31 bases near the start of the coding sequence",
        lambda g: [item(g.ins(100, "ACGT" * 7 + "ACG"), "frameshift_variant")]),
])

# --- item 11: variants at the first and the last base of a gene ---
both("B17", STD, no_setup, [
    Scn("snv_first_base", 0, "11", "diff", "SNV on the first base of the gene on the + strand",
        lambda g: [item(g.snv_g(0), "5_prime_UTR_variant" if g.fwd else "3_prime_UTR_variant")]),
    Scn("snv_last_base", 1, "11", "diff", "SNV on the last base of the gene on the + strand",
        lambda g: [item(g.snv_g(g.L - 1), "3_prime_UTR_variant" if g.fwd else "5_prime_UTR_variant")]),
])
both("B18", STD, no_setup, [
    Scn("del_anchored_on_last_base", 0, "11", "diff", "deletion of 2 bases after the last base of the gene, anchored on that base",
        lambda g: [item(g.dele_g(g.L - 1, 2), "downstream_gene_variant")]),
    Scn("ins_before_first_base", 1, "11", "diff", "insertion of 2 bases before the first base of the gene",
        lambda g: [item(g.ins_g(-1, "AC"), "upstream_gene_variant")]),
])

# --- item 12, 13, 14: prediction windows (the default radius is 50); the focal record is rare ---
both("B19", STD, no_setup, [
    Scn("mnv_at_window_edge", 0, "12", "diff", "focal SNV and an MNV of 2 bases that the edge of its window (radius 50) runs through",
        lambda g: [item(g.snv(360), "splice_region_variant"), item(g.mnv(410, 2), "missense_variant", COMMON)]),
    Scn("cplx_at_window_edge", 1, "12,13", "diff", "focal SNV and a complex indel (4 bases to 2) that the edge of its window runs through",
        lambda g: [item(g.snv(361), "splice_region_variant"), item(g.cplx(410, 4, 2), "missense_variant", COMMON)]),
])
both("B20", STD, no_setup, [
    Scn("variants_next_to_gene_start", 0, "12", "diff", "focal SNV 12 bases into the gene and a deletion of 3 bases near the gene start",
        lambda g: [item(g.dele_g(2, 3), "5_prime_UTR_variant" if g.fwd else "3_prime_UTR_variant", COMMON), item(g.snv_g(12 if g.fwd else g.L - 13), "splice_region_variant")]),
    Scn("cplx_in_window", 1, "13", "diff", "focal SNV and a complex indel (3 bases to 2) well inside its window",
        lambda g: [item(g.snv(355), "splice_region_variant"), item(g.cplx(380, 3, 2), "missense_variant", COMMON)]),
])
both("B21", STD, no_setup, [
    Scn("del_before_window", 0, "14", "diff", "focal SNV and a deletion of 3 bases anchored on the base before its window",
        lambda g: [item(g.snv(400), "splice_region_variant"), item(g.dele_g(g.rel(400) - 51, 3), "missense_variant", COMMON)]),
])

# --- item 15: phasing flags; the unphased record is on both haplotypes ---
both("B22", STD, no_setup, [
    Scn("unphased_in_window", 0, "15,20", "diff", "focal donor SNV and an unphased heterozygous SNV 20 bases away",
        lambda g: [item(g.snv(g.intron_start(0)), "splice_donor_variant"), item(g.snv(g.intron_start(0) + 20), "intron_variant", COMMON, "0/1")]),
])
both("B23", STD, no_setup, [
    Scn("too_long_in_window", 0, "15", "diff", "focal donor SNV and a deletion of 60 bases (more than the radius) in its window",
        lambda g: [item(g.snv(g.intron_start(0)), "splice_donor_variant"), item(g.dele(175, 60), "intron_variant", COMMON)]),
])

# --- item 16: exactly 50 changed bases is kept ---
both("B24", STD, no_setup, [
    Scn("del_50", 0, "16", "diff", "deletion of exactly 50 bases in the second intron",
        lambda g: [item(g.dele(480, 50), "intron_variant")]),
    Scn("ins_50", 1, "16", "diff", "insertion of exactly 50 bases in the second intron",
        lambda g: [item(g.ins(560, "ACGTACGTAC" * 5), "intron_variant")]),
])

# --- item 18: the reading frame trigger: two variants that make a stop only together ---
def plant_cca_120(g):
    g.plant_cds(120, "CCA", "codon for a pair")

both("B25", STD, plant_cca_120, [
    Scn("pair_makes_stop", 0, "18", "diff", "two SNVs in one codon (CCA to TAA) on one haplotype: a stop only together",
        lambda g: [item(g.snv(g.ct(120), "T"), "missense_variant", COMMON), item(g.snv(g.ct(121), "A"), "missense_variant", COMMON)]),
])
both("B26", STD, plant_cca_120, [
    Scn("pair_on_phase_0", 0, "18", "same", "one of the two SNVs on this haplotype",
        lambda g: [item(g.snv(g.ct(120), "T"), "missense_variant", COMMON)]),
    Scn("pair_on_phase_1", 1, "18", "same", "the other SNV on the other haplotype: no stop on either",
        lambda g: [item(g.snv(g.ct(121), "A"), "missense_variant", COMMON)]),
])

# --- a record with two ALT alleles is skipped ---
def multi_alt(g):
    t = g.ct(140)
    ref = g.genome.get(g.pos(t))
    alt1, alt2 = ROT[ref], ROT[ROT[ref]]
    return [(Rec(g.pos(t), ref, alt1 + "," + alt2), "missense_variant", RARE, None)]

copy("B27F", STD, "+", no_setup, [Scn("multi_alt", 0, "multi-ALT", "diff", "a record with two ALT alleles: skipped with a warning", multi_alt)])


# --- NMD: where the new stop is (the NMD template). Stops are made by a substitution in a planted CAA codon at these coding indices,
# which are codon starts: the class is that of the first base of the stop codon ---
STOP_PAIRS = [(99, 147), (150, 201), (300, 447), (450, 600), (798, 801), (900, 930)]


def nmd_stop_copy(pair, number):
    first, second = pair

    def setup(g):
        g.plant_cds(first, "CAA", "stop")
        g.plant_cds(second, "CAA", "stop")

    both("N%d" % number, NMD, setup, [
        Scn("stop_%d" % first, 0, "NMD", nmd_region(first), "stop_gained SNV at the coding index %d (%s)" % (first, nmd_region(first)),
            lambda g, i=first: [item(g.snv(g.ct(i), "T"), "stop_gained", None)]),
        Scn("stop_%d" % second, 1, "NMD", nmd_region(second), "stop_gained SNV at the coding index %d (%s)" % (second, nmd_region(second)),
            lambda g, i=second: [item(g.snv(g.ct(i), "T"), "stop_gained", None)]),
    ])


for number, pair in enumerate(STOP_PAIRS, start=1):
    nmd_stop_copy(pair, number)


# A frameshift deletion of the coding base d (a codon start) with TAA planted at d + 1 to d + 3: the new stop is at d, and comes from
# the base d + 1
def frameshift_setup(d1, d2):
    def setup(g):
        g.plant_cds(d1 + 1, "TAA", "frameshift stop")
        g.plant_cds(d2 + 1, "TAA", "frameshift stop")
    return setup


both("N7", NMD, frameshift_setup(219, 60), [
    Scn("frameshift_219", 0, "NMD", nmd_region(220), "frameshift (deletion of the coding base 219): a new stop at 219, from the base 220",
        lambda g: [item(g.dele(g.ct(219), 1), "frameshift_variant", None)]),
    Scn("frameshift_60", 1, "NMD", nmd_region(61), "frameshift (deletion of the coding base 60): a new stop at 60, from the base 61",
        lambda g: [item(g.dele(g.ct(60), 1), "frameshift_variant", None)]),
])


# GCC repeated from the codon 912 to 944 so that no frame has a stop before the end of the coding sequence
def no_stop_setup(g):
    g.plant_cds(912, "GCC" * 11, "no stop run")
    g.plant_cds(274, "TAA", "frameshift stop")


both("N8", NMD, no_stop_setup, [
    Scn("frameshift_921_no_stop", 0, "NMD", "none", "frameshift near the end of the coding sequence with no stop in the new frame",
        lambda g: [item(g.dele(g.ct(921), 1), "frameshift_variant", None)]),
    Scn("frameshift_273", 1, "NMD", nmd_region(274), "frameshift (deletion of the coding base 273): a new stop at 273, in the long exon",
        lambda g: [item(g.dele(g.ct(273), 1), "frameshift_variant", None)]),
])


# A lost donor, with a stop in the retained intron in the frame of the isoform. Whether the intron is retained is up to the model
def retention_setup(g):
    # Strong canonical donors, so that losing one has a large effect: the consensus CAG|GTAAGTAT around the donor of both introns.
    # The coding bases are codons that cannot be a stop (CAG after the codon start, G before the next exon's codon)
    for last_coding_index, next_index in ((249, 250), (849, 850)):
        g.plant_cds(last_coding_index - 2, "CAG", "donor consensus")
        g.plant_cds(next_index, "G", "acceptor consensus")
    g.plant_intron(0, len(DONOR_MOTIF), "AT", "donor consensus")
    g.plant_intron(1, len(DONOR_MOTIF), "AT", "donor consensus")

    g.plant_intron(0, 50, "TAG", "stop in intron 0")     # The intron after exon 0 begins at the coding index 250: 250 + 50 is in frame
    g.plant_intron(1, 47, "TAG", "stop in intron 1")     # The intron after exon 1 begins at 850: 850 + 47 = 897 is in frame


both("N9", NMD, retention_setup, [
    Scn("retained_first_intron", 0, "NMD", "nmd-model", "lost donor of exon 0, TAG in its intron: no flag if the intron is retained (not the last intron)",
        lambda g: [item(g.snv(g.intron_start(0)), "splice_donor_variant"), item(g.snv(g.intron_start(0) + 1), "splice_donor_variant")]),
    Scn("retained_last_intron", 1, "NMD", "nmd-model", "lost donor of exon 1, TAG in its intron: distal if the intron is retained (the last intron)",
        lambda g: [item(g.snv(g.intron_start(1)), "splice_donor_variant"), item(g.snv(g.intron_start(1) + 1), "splice_donor_variant")]),
])

# ------------------------------------------------------------------------------------------------------------------------------
# Assembly
# ------------------------------------------------------------------------------------------------------------------------------


def build():
    # Layout: the gene at the start of the contig, the pair of genes, then a copy after another, with flanks between
    cursor = 5
    edge_start = cursor
    cursor += EDGE.length + FLANK
    pair_start = cursor
    cursor += 399 + FLANK
    starts = []
    for (name, tpl, strand, setup, scenarios) in COPIES:
        starts.append(cursor)
        cursor += tpl.length + FLANK
    total = cursor

    genome = Genome(total)
    number = 0

    genes = []   # Every gene, for the annotation

    def new_gene(tpl, strand, name, start):
        nonlocal number
        number += 1
        gene = Gene(genome, tpl, strand, name, start, number)
        genes.append(gene)
        return gene

    edge = new_gene(EDGE, "+", "GEDGE", edge_start)
    pair1 = new_gene(SINGLE, "+", "GPAIR1", pair_start)
    pair2 = new_gene(SINGLE, "+", "GPAIR2", pair_start + 199)       # Shares its first base with the last base of the first
    assert pair1.start + pair1.L - 1 == pair2.start

    copies = []
    for index, (name, tpl, strand, setup, scenarios) in enumerate(COPIES):
        gene = new_gene(tpl, strand, "G" + name, starts[index])
        if setup is not None:
            setup(gene)
        copies.append((name, gene, scenarios))

    for gene in genes:
        gene.finish()

    # A gene that shares a base with the previous one: the genome has one base, which both agree on
    assert genome.get(pair2.start) == (pair2.seq[0] if pair2.fwd else comp(pair2.seq[0]))
    return genome, genes, edge, pair1, pair2, copies, total


class Record:
    def __init__(self, scenario_id, copy_name, tpl_name, strand, phase, rec, csq, af, gt, items, expect, desc):
        self.scenario_id, self.copy_name, self.tpl_name, self.strand = scenario_id, copy_name, tpl_name, strand
        self.phase, self.rec, self.csq, self.af, self.gt = phase, rec, csq, af, gt
        self.items, self.expect, self.desc = items, expect, desc
        self.id = scenario_id


GT = {0: "1|0", 1: "0|1"}


def check_records(records):
    # The records of one haplotype of one copy do not overlap (the anchor of an indel counts as part of it)
    for copy_name in set(r.copy_name for r in records):
        for phase in (0, 1):
            on_phase = sorted((r for r in records if r.copy_name == copy_name and (r.gt == "0/1" or r.phase == phase or r.phase is None)), key=lambda r: r.rec.pos)
            for previous, current in zip(on_phase, on_phase[1:]):
                assert previous.rec.pos + len(previous.rec.ref) <= current.rec.pos, \
                    "%s and %s overlap on the haplotype %d" % (previous.id, current.id, phase)


def csq_field(rec, csq, tids):
    alts = rec.alt.split(",")
    entries = []
    for alt in alts:
        for tid in tids:
            entries.append("%s|%s|Transcript|%s" % (alt, csq, tid))
    return ",".join(entries)


def af_field(af, n_alts):
    return ",".join(("%g" % af) for _ in range(n_alts))


# The names of bcftools csq for the consequences of VEP that differ
BCFTOOLS_NAMES = {
    "frameshift_variant": "frameshift",
    "splice_donor_variant": "splice_donor",
    "splice_acceptor_variant": "splice_acceptor",
    "splice_region_variant": "splice_region",
    "missense_variant": "missense",
}


def bcsq_field(csq, tids):
    """The same consequences as bcftools csq writes them: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change."""
    consequence = "&".join(BCFTOOLS_NAMES.get(c, c) for c in csq.split("&"))
    return ",".join("%s|GENE|%s|protein_coding|+|.|." % (consequence, tid) for tid in tids)


def vcf_line(r, tids, bcsq=False):
    rec = r.rec
    info = []
    if r.af is not None:
        info.append("AF=" + af_field(r.af, len(rec.alt.split(","))))
    if r.csq is not None:
        info.append(("BCSQ=" + bcsq_field(r.csq, tids)) if bcsq else ("CSQ=" + csq_field(rec, r.csq, tids)))
    gt = r.gt if r.gt is not None else GT[r.phase]
    return "\t".join([CONTIG, str(rec.pos + 1), r.id, rec.ref, rec.alt, ".", ".", ";".join(info) if info else ".", "GT", gt])


def vcf_header(total, bcsq=False):
    if bcsq:
        annotation = ('##INFO=<ID=BCSQ,Number=.,Type=String,Description="Haplotype-aware consequence annotation from BCFtools/csq. '
                      'Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change">')
    else:
        annotation = '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature_type|Feature">'

    return "\n".join([
        "##fileformat=VCFv4.2",
        "##contig=<ID=%s,length=%d>" % (CONTIG, total),
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        annotation,
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    ]) + "\n"


def collect(genome, copies, pair1, pair2, edge):
    records, edge_records = [], []

    def add(target, scenario_id, copy_name, tpl_name, strand, phase, items, expect, desc, built):
        for n, (rec, csq, af, gt) in enumerate(built):
            suffix = "" if len(built) == 1 else ".%d" % (n + 1)
            target.append(Record(scenario_id + suffix, copy_name, tpl_name, strand, phase, rec, csq, af, gt, items, expect, desc))

    tpl_names = {id(STD): "standard", id(NMD): "nmd", id(EDGE): "edge", id(SINGLE): "single"}

    for (name, gene, scenarios) in copies:
        for scn in scenarios:
            scenario_id = "%s_p%d_%s" % (name, scn.phase, scn.label)
            add(records, scenario_id, name, tpl_names[id(gene.tpl)], "+" if gene.fwd else "-", scn.phase, scn.items, scn.expect, scn.desc, scn.build(gene))

    # Two genes that share one base (item 17). Both are forward. The first gene's MNV ends one base before the second begins
    def pair_build(gene_a, gene_b):
        shared = gene_a.start + gene_a.L - 1
        transcripts = [gene_a.tid, gene_b.tid]
        ids = []
        ids.append(("GPAIR_p0_mnv_before_next_gene", 0, [item(gene_a.mnv(197, 2), "3_prime_UTR_variant")], "an MNV whose REF ends one base before the second gene begins", transcripts[:1]))
        ids.append(("GPAIR_p0_snv_in_second_gene", 0, [item(gene_b.snv(gene_b.ct(90)), "missense_variant")], "SNV in the coding sequence of the second gene (a block of its own or the same block)", transcripts[1:]))
        ids.append(("GPAIR_p1_snv_on_shared_base", 1, [(Rec(shared, genome.get(shared), ROT[genome.get(shared)]), None, RARE, None)], "SNV on the one base both genes have", transcripts))
        ids.append(("GPAIR_p1_snv_in_first_gene", 1, [item(gene_a.snv(gene_a.ct(90)), "missense_variant")], "SNV in the coding sequence of the first gene", transcripts[:1]))
        return ids

    pair_tids = {}
    for (sid, phase, built, desc, tids) in pair_build(pair1, pair2):
        before = len(records)
        add(records, sid, "GPAIR", "single", "+", phase, "17", "diff", desc, built)
        for r in records[before:]:
            pair_tids[r.id] = tids
            if r.csq is None:
                r.csq = "3_prime_UTR_variant&5_prime_UTR_variant"
            elif len(tids) == 2:
                pass

    # An SNV that is in no gene
    intergenic_pos = edge.start + edge.L + 700
    ref = genome.get(intergenic_pos)
    records.append(Record("INTERGENIC_p0_snv", "none", "none", "+", 0, Rec(intergenic_pos, ref, ROT[ref]), None, RARE, None, "", "same", "SNV in no gene: written unchanged"))

    # The gene at the start of the contig
    edge_items = [
        ("GEDGE_p0_donor_snv", 0, "wraparound", "info", "SNV on the first donor base 49 bases from the start of the contig", [item(edge.snv(edge.intron_start(0)), "splice_donor_variant")]),
        ("GEDGE_p0_exon_snv", 0, "wraparound", "info", "SNV in the first exon", [item(edge.snv(edge.ct(10)), "missense_variant")]),
        ("GEDGE_p1_del_near_start", 1, "12", "info", "deletion of 2 bases 8 bases into the contig, before the gene begins", [item(edge.dele_g(2, 2), "upstream_gene_variant")]),
    ]
    for (sid, phase, items, expect, desc, built) in edge_items:
        add(edge_records, sid, "GEDGE", "edge", "+", phase, items, expect, desc, built)

    return records, edge_records, pair_tids


def gff_lines(gene):
    strand = "+" if gene.fwd else "-"
    L, start = gene.L, gene.start

    def genomic(t0, t1):
        """The 1-based inclusive range on the + strand of the transcript offsets [t0, t1)."""
        if gene.fwd:
            return start + t0 + 1, start + t1
        return start + L - t1 + 1, start + L - t0

    lines = []
    lo, hi = genomic(0, L)
    lines.append("\t".join([CONTIG, "golden", "gene", str(lo), str(hi), ".", strand, ".", "ID=gene:%s;Name=%s;biotype=protein_coding" % (gene.gid, gene.name)]))
    lines.append("\t".join([CONTIG, "golden", "mRNA", str(lo), str(hi), ".", strand, ".", "ID=transcript:%s;Parent=gene:%s;biotype=protein_coding" % (gene.tid, gene.gid)]))

    # The phase is counted in the direction of the transcript; the lines are written in ascending order like the annotation of VEP
    parts = []
    coding_so_far = 0
    for k, (length, cds_offset, cds_length) in enumerate(gene.tpl.exons):
        t0, t1 = gene.tpl.exon_t[k]
        lo, hi = genomic(t0, t1)
        parts.append((lo, 0, "\t".join([CONTIG, "golden", "exon", str(lo), str(hi), ".", strand, ".", "Parent=transcript:%s" % gene.tid])))
        if cds_length:
            lo, hi = genomic(t0 + cds_offset, t0 + cds_offset + cds_length)
            phase = (3 - coding_so_far % 3) % 3          # The bases to skip to reach the next codon, as the reader checks them
            coding_so_far += cds_length
            parts.append((lo, 1, "\t".join([CONTIG, "golden", "CDS", str(lo), str(hi), ".", strand, str(phase),
                                            "Parent=transcript:%s;constitutive=%d" % (gene.tid, gene.tpl.constitutive[k])])))
    lines.extend(line for (_, _, line) in sorted(parts))
    return lines


def write_fasta(directory, genome, total):
    sequence = "".join(genome.b)
    with open(os.path.join(directory, "golden.fa"), "w") as out:
        out.write(">%s\n" % CONTIG)
        for i in range(0, len(sequence), 60):
            out.write(sequence[i:i + 60] + "\n")
    with open(os.path.join(directory, "golden.fa.fai"), "w") as out:
        out.write("%s\t%d\t%d\t60\t61\n" % (CONTIG, total, len(">%s\n" % CONTIG)))


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    os.makedirs(directory, exist_ok=True)

    genome, genes, edge, pair1, pair2, copies, total = build()
    records, edge_records, pair_tids = collect(genome, copies, pair1, pair2, edge)

    tids = {}
    for (name, gene, scenarios) in copies:
        tids[name] = [gene.tid]
    tids["GEDGE"] = [edge.tid]
    tids["none"] = []

    def tids_of(r):
        if r.id in pair_tids:
            return pair_tids[r.id]
        return tids[r.copy_name]

    # The reference of every record has to be the genome
    for r in records + edge_records:
        assert genome.seq(r.rec.pos, len(r.rec.ref)) == r.rec.ref, "%s: REF is not the reference" % r.id
        assert all(c in "ACGT," for c in r.rec.alt) and r.rec.alt != "", r.id

    check_records(records)
    check_records(edge_records)

    records.sort(key=lambda r: r.rec.pos)
    edge_records.sort(key=lambda r: r.rec.pos)

    write_fasta(directory, genome, total)

    with open(os.path.join(directory, "golden.gff3"), "w") as out:
        out.write("##gff-version 3\n")
        for gene in sorted(genes, key=lambda g: g.start):
            out.write("\n".join(gff_lines(gene)) + "\n")

    # The gene regions for cpliceai_reference: 0-based half-open (BED), as src/gene_regions.c reads them. The extent of a gene is that of
    # its exons here, which is also what the annotation of rescue_ranger gives the gene
    with open(os.path.join(directory, "golden.regions.tsv"), "w") as out:
        out.write("NAME\tCHROM\tSTRAND\tTX_START\tTX_END\n")
        for gene in sorted(genes, key=lambda g: g.start):
            out.write("%s\t%s\t%s\t%d\t%d\n" % (gene.name, CONTIG, "+" if gene.fwd else "-", gene.start, gene.start + gene.L))

    header = vcf_header(total)
    with open(os.path.join(directory, "golden.vcf"), "w") as out:
        out.write(header)
        for r in records:
            out.write(vcf_line(r, tids_of(r)) + "\n")

    # The same records annotated like bcftools csq does (INFO/BCSQ); run with --annotation bcftools, the result is the same as for golden.vcf
    with open(os.path.join(directory, "golden_bcsq.vcf"), "w") as out:
        out.write(vcf_header(total, bcsq=True))
        for r in records:
            out.write(vcf_line(r, tids_of(r), bcsq=True) + "\n")

    with open(os.path.join(directory, "golden_edge.vcf"), "w") as out:
        out.write(header)
        for r in edge_records:
            out.write(vcf_line(r, tids_of(r)) + "\n")

    # The read failure: valid records, then a line that stops after the 4th column
    coding = [r for r in records if r.csq is not None][:2]
    with open(os.path.join(directory, "golden_broken.vcf"), "w") as out:
        out.write(header)
        for r in coding:
            out.write(vcf_line(r, tids_of(r)) + "\n")
        out.write("%s\t%d\t.\tA\n" % (CONTIG, coding[-1].rec.pos + 2000))

    with open(os.path.join(directory, "manifest.tsv"), "w") as out:
        out.write("id\tcopy\ttemplate\tstrand\tphase\tpos\tref\talt\tgt\taf\tcsq\titems\texpect\tdescription\n")
        for r in records + edge_records:
            gt = r.gt if r.gt is not None else GT[r.phase]
            out.write("\t".join([r.id, r.copy_name, r.tpl_name, r.strand, str(r.phase), str(r.rec.pos + 1), r.rec.ref, r.rec.alt, gt,
                                 "" if r.af is None else "%g" % r.af, r.csq or "", r.items, r.expect, r.desc]) + "\n")

    print("Wrote %s: %d genes, %d records (%d at the contig start), %d bases" % (directory, len(genes), len(records), len(edge_records), total))


if __name__ == "__main__":
    main()
