#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>

#include "bcftools/gff.h"
#include "klib/kvec.h"
#include "record.h"
#include "record_block.h"
#include "session.h"

#include "check_record_block.h"

/*
* classify_read: what the result of bcf_read means. 0 is a record, -1 without an error code is the end of the input, and every
* other negative result, or an error code on the record, is a failure. A read error used to be taken for the end of the input.
*/

START_TEST(check_classify_read__a_record) {
    ck_assert_int_eq(classify_read(0, 0), READ_OK);
} END_TEST

START_TEST(check_classify_read__the_end_of_the_input) {
    ck_assert_int_eq(classify_read(-1, 0), READ_END);
} END_TEST

START_TEST(check_classify_read__a_parse_error_sets_the_error_code) {
    ck_assert_int_eq(classify_read(-1, 1), READ_FAILED);
    ck_assert_int_eq(classify_read(-2, 4), READ_FAILED);
} END_TEST

// A read error of the file (a truncated compressed file for example) is a result below -1 and no error code on the record
START_TEST(check_classify_read__a_read_error_is_not_the_end) {
    ck_assert_int_eq(classify_read(-2, 0), READ_FAILED);
    ck_assert_int_eq(classify_read(-5, 0), READ_FAILED);
} END_TEST

/*
* RecordBuffer_fill: the next record into the buffer.
*/

START_TEST(check_record_buffer_fill__a_record_that_is_waiting_is_not_replaced) {
    RecordBuffer buffer = { .v = NULL, .unused = true, .failed = false };

    ck_assert_int_eq(RecordBuffer_fill(NULL, NULL, &buffer), READ_OK); // Nothing is read: the file is not touched
} END_TEST

START_TEST(check_record_buffer_fill__after_a_failure_nothing_is_read_again) {
    RecordBuffer buffer = { .v = NULL, .unused = false, .failed = true };

    ck_assert_int_eq(RecordBuffer_fill(NULL, NULL, &buffer), READ_FAILED);
    ck_assert_int_eq(RecordBuffer_fill(NULL, NULL, &buffer), READ_FAILED);
} END_TEST

// A VCF with the given text after the header, in a file of its own that is removed by close_vcf
typedef struct VcfFile {
    char path[64];
    htsFile *fp;
    bcf_hdr_t *hdr;
} VcfFile;

static VcfFile open_vcf(const char *records) {
    VcfFile f;
    strcpy(f.path, "check_record_block_XXXXXX");
    int fd = mkstemp(f.path);
    ck_assert_int_ge(fd, 0);

    FILE *out = fdopen(fd, "w");
    fputs("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", out);
    fputs(records, out);
    fclose(out);

    f.fp = hts_open(f.path, "r");
    ck_assert_ptr_nonnull(f.fp);
    f.hdr = bcf_hdr_read(f.fp);
    ck_assert_ptr_nonnull(f.hdr);

    return f;
}

static void close_vcf(VcfFile *f) {
    bcf_hdr_destroy(f->hdr);
    hts_close(f->fp);
    unlink(f->path);
}

START_TEST(check_record_buffer_fill__records_then_the_end) {
    VcfFile f = open_vcf("chr1\t10\t.\tA\tG\t.\t.\t.\nchr1\t20\t.\tC\tT\t.\t.\t.\n");
    RecordBuffer buffer = { .v = bcf_init(), .unused = false, .failed = false };

    ck_assert_int_eq(RecordBuffer_fill(f.fp, f.hdr, &buffer), READ_OK);
    ck_assert_int_eq(buffer.v->pos, 9);
    ck_assert(buffer.unused);

    buffer.unused = false; // Processed
    ck_assert_int_eq(RecordBuffer_fill(f.fp, f.hdr, &buffer), READ_OK);
    ck_assert_int_eq(buffer.v->pos, 19);

    buffer.unused = false;
    ck_assert_int_eq(RecordBuffer_fill(f.fp, f.hdr, &buffer), READ_END);
    ck_assert(!buffer.failed);
    ck_assert_int_eq(RecordBuffer_fill(f.fp, f.hdr, &buffer), READ_END); // And it stays the end

    bcf_destroy(buffer.v);
    close_vcf(&f);
} END_TEST

// A line that cannot be parsed (too few columns) is a failure, not the end, and the buffer stays failed
START_TEST(check_record_buffer_fill__an_unparsable_line_is_a_failure_that_stays) {
    VcfFile f = open_vcf("chr1\t10\t.\tA\tG\t.\t.\t.\nchr1\t20\t.\tC\n");
    RecordBuffer buffer = { .v = bcf_init(), .unused = false, .failed = false };

    ck_assert_int_eq(RecordBuffer_fill(f.fp, f.hdr, &buffer), READ_OK);
    buffer.unused = false;

    ck_assert_int_eq(RecordBuffer_fill(f.fp, f.hdr, &buffer), READ_FAILED);
    ck_assert(buffer.failed);
    ck_assert_int_eq(RecordBuffer_fill(f.fp, f.hdr, &buffer), READ_FAILED);

    bcf_destroy(buffer.v);
    close_vcf(&f);
} END_TEST

/*
* RecordBlock_next and the phase copies, on a small genome written to files: chr1 of 800 bases with two protein-coding genes on the
* forward strand, G1 at 101-300 (GFF, 1-based, so [100, 300) as a Range) and G2 at 501-700, each with one transcript (T1, T2) that
* is one exon with the CDS 111-290 and 511-690. The session holds what the loader reads (GFF, FASTA, VCF and the VEP format);
* the models and the output are not used by the loader.
*/
#define GFF_G1 "chr1\tsrc\tgene\t101\t300\t.\t+\t.\tID=gene:G1;Name=GENE1;biotype=protein_coding\n" \
               "chr1\tsrc\tmRNA\t101\t300\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding\n" \
               "chr1\tsrc\texon\t101\t300\t.\t+\t.\tParent=transcript:T1\n" \
               "chr1\tsrc\tCDS\t111\t290\t.\t+\t0\tParent=transcript:T1\n"
#define GFF_G2 "chr1\tsrc\tgene\t501\t700\t.\t+\t.\tID=gene:G2;Name=GENE2;biotype=protein_coding\n" \
               "chr1\tsrc\tmRNA\t501\t700\t.\t+\t.\tID=transcript:T2;Parent=gene:G2;biotype=protein_coding\n" \
               "chr1\tsrc\texon\t501\t700\t.\t+\t.\tParent=transcript:T2\n" \
               "chr1\tsrc\tCDS\t511\t690\t.\t+\t0\tParent=transcript:T2\n"

#define VCF_HEADER "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=800>\n" \
    "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|Feature_type|Feature\">\n" \
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

// Two records in G1, one between the genes, a record in G2 and a record in G2 with two ALT alleles
#define VCF_RECORDS "chr1\t151\t.\tA\tG\t.\t.\tCSQ=G|missense_variant|Transcript|T1\n" \
                    "chr1\t161\t.\tA\tT\t.\t.\tCSQ=T|missense_variant|Transcript|T1\n" \
                    "chr1\t401\t.\tA\tG\t.\t.\t.\n" \
                    "chr1\t551\t.\tA\tG\t.\t.\tCSQ=G|missense_variant|Transcript|T2\n" \
                    "chr1\t561\t.\tA\tC,T\t.\t.\tCSQ=C|missense_variant|Transcript|T2,T|missense_variant|Transcript|T2\n"

typedef struct Loader {
    char gff_path[64], fasta_path[64], vcf_path[64];
    Session s;
    RecordBuffer buffer;
} Loader;

static void write_temp_file(char *path, const char *prefix, const char *content) {
    snprintf(path, 64, "%s_XXXXXX", prefix);
    const int fd = mkstemp(path);
    ck_assert_int_ge(fd, 0);

    FILE *out = fdopen(fd, "w");
    ck_assert_ptr_nonnull(out);
    fputs(content, out);
    fclose(out);
}

// The Loader is on the heap because the GFF keeps a pointer to the path of its file
static Loader *open_loader(const char *vcf_records) {
    Loader *l = calloc(1, sizeof(Loader));

    char fasta[900] = ">chr1\n";
    const size_t header = strlen(fasta);
    memset(fasta + header, 'A', 800);
    strcpy(fasta + header + 800, "\n");

    char vcf[2048];
    snprintf(vcf, sizeof(vcf), "%s%s", VCF_HEADER, vcf_records);

    write_temp_file(l->gff_path, "check_loader_gff", "##gff-version 3\n" GFF_G1 GFF_G2);
    write_temp_file(l->fasta_path, "check_loader_fa", fasta);
    write_temp_file(l->vcf_path, "check_loader_vcf", vcf);

    l->s.gff = gff_init(l->gff_path);
    gff_parse(l->s.gff);
    l->s.fai = fai_load(l->fasta_path); // Writes the index next to the FASTA
    ck_assert_ptr_nonnull(l->s.fai);
    l->s.fp_in = hts_open(l->vcf_path, "r");
    ck_assert_ptr_nonnull(l->s.fp_in);
    l->s.hdr = bcf_hdr_read(l->s.fp_in);
    ck_assert_ptr_nonnull(l->s.hdr);
    l->s.vep_fmt = parse_csq_hdr_fmt(l->s.hdr, ANNOTATION_VEP, NULL);
    ck_assert_ptr_nonnull(l->s.vep_fmt);
    l->s.prediction_radius = 50;

    l->buffer = (RecordBuffer) { .v = bcf_init(), .unused = false, .failed = false };

    return l;
}

static void close_loader(Loader *l) {
    char fai_path[80];
    snprintf(fai_path, sizeof(fai_path), "%s.fai", l->fasta_path);

    bcf_destroy(l->buffer.v);
    Session_close(&l->s);
    unlink(l->gff_path); unlink(l->fasta_path); unlink(fai_path); unlink(l->vcf_path);
    free(l);
}

// The next block, which has to be there
static RecordBlock next_block(Loader *l) {
    RecordBlock block;
    ck_assert_int_eq(RecordBlock_next(&l->s, &l->buffer, &block), 0);
    return block;
}

START_TEST(check_record_block_next__records_in_one_gene_form_a_block_with_a_shared_annotation) {
    Loader *l = open_loader(VCF_RECORDS);
    RecordBlock block = next_block(l);

    ck_assert_int_eq(block.is_coding, 1);
    ck_assert_int_eq(kv_size(block.records), 2);
    ck_assert_ptr_nonnull(block.annotation);
    ck_assert_int_eq(block.annotation->n_records, 2);

    Record *first = kv_A(block.records, 0), *second = kv_A(block.records, 1);
    ck_assert_int_eq(first->v->pos, 150);   // 0-based
    ck_assert_int_eq(second->v->pos, 160);

    // The state that depends on the annotation only is there once, without the per-phase copies of the transcripts
    Record *records[2] = { first, second };
    for (int i = 0; i < 2; i++) {
        Record *r = records[i];
        ck_assert_ptr_nonnull(r->tscripts);
        ck_assert_int_eq(kv_size(*r->tscripts), 1);
        ck_assert_ptr_null(kv_A(*r->tscripts, 0));
        ck_assert_int_eq(kv_size(*r->genes), 1);
        ck_assert_ptr_nonnull(r->transcript_genes);
        ck_assert_ptr_nonnull(r->native_scores);
        ck_assert_ptr_nonnull(r->diff_scores);
        ck_assert_ptr_nonnull(r->haplotype_scores);
        ck_assert_ptr_nonnull(r->isoforms);
        ck_assert_ptr_nonnull(r->flags);
        ck_assert_int_eq(r->skipped, 0);
    }

    // The gene and the transcript are shared by the records of the block, not copied
    Gene *gene = kv_A(*first->genes, 0);
    ck_assert_str_eq(gene->name, "GENE1");
    ck_assert_int_eq(gene->coords.start, 100);   // Half-open
    ck_assert_int_eq(gene->coords.end, 300);
    ck_assert_int_eq(gene->seq.l, 200);
    ck_assert_ptr_eq(gene, kv_A(*second->genes, 0));
    ck_assert_ptr_eq(gene, kv_A(*first->transcript_genes, 0));
    ck_assert_ptr_nonnull(kv_A(block.annotation->records[0].tscripts, 0));
    ck_assert_ptr_eq(kv_A(block.annotation->records[0].tscripts, 0), kv_A(block.annotation->records[1].tscripts, 0));

    RecordBlock_destroy(block);
    close_loader(l);
} END_TEST

START_TEST(check_record_block_next__a_record_outside_every_gene_is_a_block_of_its_own) {
    Loader *l = open_loader(VCF_RECORDS);
    RecordBlock_destroy(next_block(l));

    RecordBlock block = next_block(l);

    ck_assert_int_eq(block.is_coding, 0);
    ck_assert_int_eq(kv_size(block.records), 1);
    ck_assert_ptr_null(block.annotation);
    ck_assert_int_eq(kv_A(block.records, 0)->v->pos, 400);
    ck_assert_ptr_null(kv_A(block.records, 0)->tscripts);   // Nothing to skip later: it has no coding state

    RecordBlock_destroy(block);
    close_loader(l);
} END_TEST

START_TEST(check_record_block_next__a_record_with_two_alt_alleles_is_skipped) {
    Loader *l = open_loader(VCF_RECORDS);
    RecordBlock_destroy(next_block(l));
    RecordBlock_destroy(next_block(l));

    RecordBlock block = next_block(l);

    ck_assert_int_eq(block.is_coding, 1);
    ck_assert_int_eq(kv_size(block.records), 2);

    // The record of G2 is set up like any other. The one with two ALT alleles has no coding state, so it takes no part in a
    // phase, and it is marked as skipped
    Record *plain = kv_A(block.records, 0), *multi = kv_A(block.records, 1);
    ck_assert_int_eq(plain->v->pos, 550);
    ck_assert_ptr_nonnull(plain->tscripts);
    ck_assert_int_eq(plain->skipped, 0);

    ck_assert_int_eq(multi->v->pos, 560);
    ck_assert_int_eq(multi->skipped, 1);
    ck_assert_ptr_null(multi->tscripts);
    ck_assert_ptr_null(multi->genes);
    ck_assert_ptr_null(multi->flags);
    ck_assert_ptr_null(multi->isoforms);

    // The phases leave it alone
    RecordBlock_begin_phase(&block);
    ck_assert_ptr_nonnull(kv_A(*plain->tscripts, 0));
    ck_assert_ptr_null(multi->tscripts);
    RecordBlock_end_phase(&block);

    RecordBlock_destroy(block);
    close_loader(l);
} END_TEST

START_TEST(check_record_block_next__after_the_last_block_it_is_the_end_of_the_input) {
    Loader *l = open_loader(VCF_RECORDS);
    for (int i = 0; i < 3; i++) RecordBlock_destroy(next_block(l));

    RecordBlock block;
    ck_assert_int_eq(RecordBlock_next(&l->s, &l->buffer, &block), 1);
    ck_assert_int_eq(RecordBlock_next(&l->s, &l->buffer, &block), 1);   // And it stays the end

    close_loader(l);
} END_TEST

// The records before a read failure come first, as a block, and the failure with the next call
START_TEST(check_record_block_next__the_records_before_a_read_failure_come_first) {
    Loader *l = open_loader("chr1\t151\t.\tA\tG\t.\t.\tCSQ=G|missense_variant|Transcript|T1\n"
                            "chr1\t161\t.\tA\n");
    RecordBlock block = next_block(l);

    ck_assert_int_eq(kv_size(block.records), 1);
    ck_assert_int_eq(kv_A(block.records, 0)->v->pos, 150);
    RecordBlock_destroy(block);

    ck_assert_int_eq(RecordBlock_next(&l->s, &l->buffer, &block), 2);
    ck_assert_int_eq(RecordBlock_next(&l->s, &l->buffer, &block), 2);   // The failure stays

    close_loader(l);
} END_TEST

START_TEST(check_record_block_phase__every_record_gets_its_own_copies_that_end_with_the_phase) {
    Loader *l = open_loader(VCF_RECORDS);
    RecordBlock block = next_block(l);
    Record *first = kv_A(block.records, 0), *second = kv_A(block.records, 1);
    Transcript *shared = kv_A(block.annotation->records[0].tscripts, 0);

    // Twice, to see that a phase starts from nothing
    for (int phase = 0; phase < 2; phase++) {
        RecordBlock_begin_phase(&block);

        Transcript *copy_of_first = kv_A(*first->tscripts, 0), *copy_of_second = kv_A(*second->tscripts, 0);
        ck_assert_ptr_nonnull(copy_of_first);
        ck_assert_ptr_nonnull(copy_of_second);
        ck_assert_ptr_ne(copy_of_first, shared);
        ck_assert_ptr_ne(copy_of_first, copy_of_second);

        RecordBlock_end_phase(&block);

        ck_assert_ptr_null(kv_A(*first->tscripts, 0));
        ck_assert_ptr_null(kv_A(*second->tscripts, 0));
        ck_assert_ptr_eq(kv_A(block.annotation->records[0].tscripts, 0), shared);   // The annotation is not touched
    }

    RecordBlock_destroy(block);
    close_loader(l);
} END_TEST

Suite *record_block_suite() {
    Suite *s = suite_create("Record block");

    TCase *tc_classify = tcase_create("classify_read");
    tcase_add_test(tc_classify, check_classify_read__a_record);
    tcase_add_test(tc_classify, check_classify_read__the_end_of_the_input);
    tcase_add_test(tc_classify, check_classify_read__a_parse_error_sets_the_error_code);
    tcase_add_test(tc_classify, check_classify_read__a_read_error_is_not_the_end);
    suite_add_tcase(s, tc_classify);

    TCase *tc_fill = tcase_create("RecordBuffer_fill");
    tcase_add_test(tc_fill, check_record_buffer_fill__a_record_that_is_waiting_is_not_replaced);
    tcase_add_test(tc_fill, check_record_buffer_fill__after_a_failure_nothing_is_read_again);
    tcase_add_test(tc_fill, check_record_buffer_fill__records_then_the_end);
    tcase_add_test(tc_fill, check_record_buffer_fill__an_unparsable_line_is_a_failure_that_stays);
    suite_add_tcase(s, tc_fill);

    TCase *tc_next = tcase_create("RecordBlock_next");
    tcase_add_test(tc_next, check_record_block_next__records_in_one_gene_form_a_block_with_a_shared_annotation);
    tcase_add_test(tc_next, check_record_block_next__a_record_outside_every_gene_is_a_block_of_its_own);
    tcase_add_test(tc_next, check_record_block_next__a_record_with_two_alt_alleles_is_skipped);
    tcase_add_test(tc_next, check_record_block_next__after_the_last_block_it_is_the_end_of_the_input);
    tcase_add_test(tc_next, check_record_block_next__the_records_before_a_read_failure_come_first);
    suite_add_tcase(s, tc_next);

    TCase *tc_phase = tcase_create("RecordBlock_phase");
    tcase_add_test(tc_phase, check_record_block_phase__every_record_gets_its_own_copies_that_end_with_the_phase);
    suite_add_tcase(s, tc_phase);

    return s;
}
