#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <klib/kvec.h>
#include <klib/kstring.h>

#include <time.h>
#include <zlib.h>
#include <htslib/hts.h>
#include <htslib/khash.h>
#include <htslib/khash_str2int.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
// #include <cpliceai/predict.h>
// #include <cpliceai/utils.h>

#include "bcftools/gff.h"
#include "bcftools/regidx.h"

#include "logging/log.h"
#include "region.h"

#define BUFFER_SIZE 2024
#define ID 0
#define TYPE 2
#define BEGIN 3
#define END 4
#define DIRECTION 6
#define PHASE 7
#define ATTRIBUTE 8

#define VCF_DELIM ','
#define CSQ_FIELD "CSQ"
#define HDR_DESC "Description"
#define VEP_DELIM '|'
#define VEP_DESC_PREFIX "Format: "
#define VEP_DESC_PREFIX_LEN 8
#define CONSEQUENCE_DELIM '&'

#define VEP_CONSEQUENCE "Consequence"
#define VEP_FEATURE_TYPE "Feature_type"
#define VEP_FEATURE "Feature"
#define VEP_FEATURE_TYPE_TRANSCRIPT "Transcript"

#define SPLICE_DONOR_VARIANT "splice_donor_variant"
#define SPLICE_ACCEPTOR_VARIANT "splice_acceptor_variant"
#define FRAMESHIFT_VARIANT "frameshift_variant"
#define STOP_GAINED_VARIANT "stop_gained"

#define N_CSQ_FIELDS 29
#define NUM_PHASES 2

#define MAX_FRAMESHIFT_RADIUS 80

#define CAI_DISTANCE 50
#define WIDTH 10101
#define VARIANT_POS 5050
#define SPLICEAI_THRESHOLD 0.8
#define SPLICEAI_RESCUE_THRESHOLD 0.2

#define SPLICE_TRIGGER 0
#define FRAMESHIFT_TRIGGER 1
#define STOP_TRIGGER 2

#define STOP_CODONS {"TAG", "TGA", "TAA"}
#define COMPLEMENTARY_STOP_CODONS {"ATC", "ACT", "ATT"}

#define SPLICEAI_HAP_TAG "SpliceAI_Haplotype"
#define SPLICEAI_HAP_DESC "##INFO=<ID=SpliceAI_Haplotype,Number=.,Type=String,Description=\"CpliceAIv1.0 haplotype variant annotation. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL\">"

#define CS_TAG "CS"
#define CS_DESC "##INFO=<ID=CS,Number=.,Type=String,Description=\"Change in splicing between haplotype ref and haplotype alt at positions DP_AG, DP_AL, DP_DG, DP_DL. Format: CS_AG|CS_AL|CS_DG|CS_DL|CS_SUM\">"

#define RESCUE_TAG "RESCUE"
#define RESCUE_DESC "##INFO=<ID=RESCUE,Number=.,Type=String,Description=\"Boolean rescue prediction\">"

#define RESCUE_TYPE_TAG "RESCUE_TYPE"
#define RESCUE_TYPE_DESC "##INFO=<ID=RESCUE_TYPE,Number=.,Type=String,Description=\"Rescue type\">"

#define PTV_TAG "PTV"
#define PTV_DESC "##INFO=<ID=PTV,Number=.,Type=String,Description=\"Newly predicted PTV consequence\">"

#define PTV_TYPE_TAG "PTV_TYPE"
#define PTV_TYPE_DESC "##INFO=<ID=PTV_TYPE,Number=.,Type=String,Description=\"Type of PTV\">"

#define PREDICTION_FLAGS_TAG "PREDICTION_FLAGS"
#define PREDICTION_FLAGS_DESC "##INFO=<ID=PREDICTION_FLAGS,Number=.,Type=String,Description=\"Additional flags for the context in which the prediction happened.\">"

#define INFO_TAGS { SPLICEAI_HAP_TAG, CS_TAG, RESCUE_TAG, RESCUE_TYPE_TAG, PTV_TAG, PTV_TYPE_TAG, PREDICTION_FLAGS_TAG }

typedef struct
{
    uint32_t start, end;
} Range;

enum Trigger {
    SPLICE,
    FRAMESHIFT,
    STOP
};

typedef kvec_t(bcf1_t *) Variants;

typedef struct Gene {
    uint32_t gff_id;
    Range coords;
    kstring_t seq;
    uint32_t strand;
} Gene;

typedef kvec_t(Gene *) Genes;

typedef struct ConsequencePredictions {
    kstring_t **hap_score;
    kstring_t **splice_change;
    kstring_t **is_rescue;
    kstring_t **rescue_type;
    kstring_t **is_ptv;
    kstring_t **ptv_type;
    kstring_t **flags;
    size_t num_consequences;
} ConsequencePredictions;

typedef struct Transcript {
    int rid;
    gf_tscript_t *definition;
} Transcript;

typedef kvec_t(Transcript) Transcripts;

typedef kvec_t(gf_tscript_t *) gf_tscripts;

typedef struct TscriptCsqs {
    char *tid;
    int num_consequences;
    char **consequences;
} TscriptCsqs;

typedef kvec_t(TscriptCsqs) VepCsqs;

typedef struct Record {
    bcf1_t *v;
    VepCsqs *vep;
    enum Trigger triggers[3];
    ConsequencePredictions *c;
    gf_tscripts *tscripts;
    Genes *genes;
} Record;

typedef kvec_t(Record *) Records;

typedef struct RecordBlock {
    int is_coding;
    Records *records;
} RecordBlock;

typedef kvec_t(Records *) PhasedRecords;

typedef struct HdrFormat {
    char *data_; // Only used for storing data, access this data should be done through id2str
    char **id2str;
    kh_str2int_t *str2id;
} HdrFormat;



KHASH_MAP_INIT_STR(trCsqHash, char **);


void init_consequence_predictions(ConsequencePredictions **dst, int num_consequences) {
    ConsequencePredictions *csq_preds = malloc(sizeof(ConsequencePredictions));
    csq_preds->num_consequences = num_consequences;
    csq_preds->hap_score = malloc(num_consequences * sizeof(kstring_t *));
    csq_preds->splice_change = malloc(num_consequences * sizeof(kstring_t *));
    csq_preds->is_rescue = malloc(num_consequences * sizeof(kstring_t *));
    csq_preds->rescue_type = malloc(num_consequences * sizeof(kstring_t *));
    csq_preds->is_ptv = malloc(num_consequences * sizeof(kstring_t *));
    csq_preds->ptv_type = malloc(num_consequences * sizeof(kstring_t *));
    csq_preds->flags = malloc(num_consequences * sizeof(kstring_t *));


    for (int i = 0; i < num_consequences; i++) {
        csq_preds->hap_score[i] = calloc(1, sizeof(kstring_t));
        csq_preds->splice_change[i] = calloc(1, sizeof(kstring_t));
        csq_preds->is_rescue[i] = calloc(1, sizeof(kstring_t));
        csq_preds->rescue_type[i] = calloc(1, sizeof(kstring_t));
        csq_preds->is_ptv[i] = calloc(1, sizeof(kstring_t));
        csq_preds->ptv_type[i] = calloc(1, sizeof(kstring_t));
        csq_preds->flags[i] = calloc(1, sizeof(kstring_t));
    }

    *dst = csq_preds;
}

void destroy_consequence_predictions(ConsequencePredictions *csq_preds) {
    for (int i = 0; i < csq_preds->num_consequences; i++) {
        if (csq_preds->hap_score[i]->m > 0) free(csq_preds->hap_score[i]->s);
        free(csq_preds->hap_score[i]);

        if (csq_preds->splice_change[i]->m > 0) free(csq_preds->splice_change[i]->s);
        free(csq_preds->splice_change[i]);

        if (csq_preds->is_rescue[i]->m > 0) free(csq_preds->is_rescue[i]->s);
        free(csq_preds->is_rescue[i]);

        if (csq_preds->rescue_type[i]->m > 0) free(csq_preds->rescue_type[i]->s);
        free(csq_preds->rescue_type[i]);

        if (csq_preds->is_ptv[i]->m > 0) free(csq_preds->is_ptv[i]->s);
        free(csq_preds->is_ptv[i]);

        if (csq_preds->ptv_type[i]->m > 0) free(csq_preds->ptv_type[i]->s);
        free(csq_preds->ptv_type[i]);

        if (csq_preds->flags[i]->m > 0) free(csq_preds->flags[i]->s);
        free(csq_preds->flags[i]);
    }

    free(csq_preds->hap_score);
    free(csq_preds->splice_change);
    free(csq_preds->is_rescue);
    free(csq_preds->rescue_type);
    free(csq_preds->is_ptv);
    free(csq_preds->ptv_type);
    free(csq_preds->flags);

    free(csq_preds);
}

/*
 *
 * WARN: Does not add last field if it is empty (e.g., 28 fields instead of 29 for VEP)
 */
int *str_split(kstring_t *s, char delim, int *num_splits) {
    int i = 0;
    kvec_t(int) field_idxs; kv_init(field_idxs);
    while (i < s->l) {
        kv_push(int, field_idxs, i);

        while(*(s->s + i) && *(s->s + i) != delim ) i++;

        s->s[i++] = '\0';
    }

    *num_splits = field_idxs.n;
    return field_idxs.a;
}

/*
*
* Usage: 
*   int id; khash_str2int_get(idHash, "Consequence", &id);
*   fields[id];
*/
HdrFormat *parse_csq_format(bcf_hdr_t *hdr) {
    char *hdr_info;
    if (bcf_hdr_id2int(hdr, BCF_DT_ID, CSQ_FIELD) < 0) return NULL;

    // Get CSQ Description
    bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, BCF_HL_INFO, NULL, CSQ_FIELD, NULL);
    int key = bcf_hrec_find_key(hrec, HDR_DESC);
    const char *description = strstr(hrec->vals[key], VEP_DESC_PREFIX) + VEP_DESC_PREFIX_LEN; // Skip "Format: "
    int len = strlen(description)-1;

    // Split CSQ description
    kstring_t s = {0}; kputsn(description, len, &s);
    int *field_idxs, n; field_idxs = ksplit(&s, VEP_DELIM, &n);

    khiter_t k;
    kh_str2int_t *idHash = kh_init_str2int();
    char **fields = malloc(n * sizeof(char *));
    int ret;
    for (int i = 0; i < n; i++) {
        k = kh_put_str2int(idHash, s.s + field_idxs[i], &ret);
        khash_str2int_set(idHash, s.s + field_idxs[i], i);
        fields[i] = s.s + field_idxs[i];
    }
    free(field_idxs);

    HdrFormat *hdr_fmt = malloc(sizeof(HdrFormat));
    *hdr_fmt = (HdrFormat) {s.s, fields, idHash};
    return hdr_fmt;
}

void destroy_hdr_fmt(HdrFormat *hdr_fmt) {
    kh_destroy_str2int(hdr_fmt->str2id);
    free(hdr_fmt->id2str);
    free(hdr_fmt->data_);
    free(hdr_fmt);
}

VepCsqs *parse_vep_transcript_consequences(bcf1_t *v, const HdrFormat *vep_fmt, bcf_hdr_t *hdr) {
    char *vep_str = NULL; int n_vep_values = 0;
    bcf_get_info_string(hdr, v, CSQ_FIELD, &vep_str, &n_vep_values);

    VepCsqs *transcripts = malloc(sizeof(VepCsqs)); kv_init(*transcripts);

    kstring_t s = {0}; kputs(vep_str, &s);
    int *field_idxs, n; field_idxs = str_split(&s, VCF_DELIM, &n);
    for (int i = 0; i < n; i++) {
        int flen = strlen(s.s + field_idxs[i]);
        kstring_t transcript_csq = { flen, flen+1, s.s + field_idxs[i] }; // WARN: NOT MALLOC'ED, SHOULD NOT BE FREED
        int *vep_field_idxs, vep_n; vep_field_idxs = str_split(&transcript_csq, VEP_DELIM, &vep_n);

        int feature_type_id; khash_str2int_get(vep_fmt->str2id, VEP_FEATURE_TYPE, &feature_type_id);
        if (strcmp(transcript_csq.s + vep_field_idxs[feature_type_id], VEP_FEATURE_TYPE_TRANSCRIPT) != 0) { free(vep_field_idxs); continue; }

        int feature_id; khash_str2int_get(vep_fmt->str2id, VEP_FEATURE, &feature_id);
        int consequences_id; khash_str2int_get(vep_fmt->str2id, VEP_CONSEQUENCE, &consequences_id);

        // Parse consequences for transcript, vep allows multiple with '&' separators
        int clen = strlen(transcript_csq.s + vep_field_idxs[consequences_id]);
        kstring_t consequences = { clen, clen+1, transcript_csq.s + vep_field_idxs[consequences_id] };
        int *consequence_field_idxs, consequence_n; consequence_field_idxs = str_split(&consequences, CONSEQUENCE_DELIM, &consequence_n);
        char **cons = malloc(consequence_n * sizeof(char *));
        for (int j = 0; j < consequence_n; j++) {
            cons[j] = strdup(consequences.s + consequence_field_idxs[j]);
        }

        // Store in hashmap
        char *tr_id = transcript_csq.s + vep_field_idxs[feature_id];
        TscriptCsqs t = {strdup(tr_id), consequence_n, cons};
        kv_push(TscriptCsqs, *transcripts, t);

        free(consequence_field_idxs);
        free(vep_field_idxs);
    }
    free(s.s);
    free(vep_str);
    free(field_idxs);

    return transcripts;
}

void destroy_vep_csqs(VepCsqs *vep_csqs) {
    for (int i = 0; i < vep_csqs->n; i++) {
        TscriptCsqs vep_csq_field = kv_A(*vep_csqs, i);

        for (int j = 0; j < vep_csq_field.num_consequences; j++) {
            free(vep_csq_field.consequences[j]);
        }
        free(vep_csq_field.consequences);
        free(vep_csq_field.tid);
    }

    kv_destroy(*vep_csqs);
    free(vep_csqs);
}

Record *init_record(bcf1_t *v, VepCsqs *csqs) {
    Record *r = malloc(sizeof(Record));

    bcf1_t *v_dup = bcf_dup(v);
    bcf_unpack(v_dup, BCF_UN_ALL);
    r->v = v_dup;

    r->triggers[0] = 0;
    r->triggers[1] = 0;
    r->triggers[2] = 0;

    r->c = NULL;
    r->genes = NULL;
    r->tscripts = NULL;
    r->vep = NULL;

    if (csqs == NULL) {
        return r;
    }

    r->vep = csqs;
    init_consequence_predictions(&r->c, csqs->n);

    return r;
}

void destroy_record(Record *r) {
    bcf_destroy(r->v);
    if (r->vep != NULL) {
        destroy_consequence_predictions(r->c);
        destroy_vep_csqs(r->vep);
    }

    if (r->tscripts != NULL) {
        kv_destroy(*r->tscripts);
        free(r->tscripts);
    }

    if (r->genes != NULL) {
        kv_destroy(*r->genes);
        free(r->genes);
    }

    free(r);
}


// int has_overlapping_transcripts(bcf1_t *v, bcf_hdr_t *hdr, gff_t *gff) {
//     regidx_t *tscripts = gff_get(gff, idx_tscript);
//     regitr_t *itr = regitr_init(tscripts);
//
//     int has_overlap = regidx_overlap(tscripts, bcf_hdr_id2name(hdr, v->rid), v->pos, v->pos + v->rlen, itr);
//
//     regitr_destroy(itr);
//     return has_overlap;
// }

int find_gene_region_for_variant(bcf1_t *v, bcf_hdr_t *hdr, gff_t *gff, IRegion *region) {
    regidx_t *genes = gff_get(gff, idx_gene);
    regitr_t *itr = regitr_init(genes);
    const char *chr = bcf_hdr_id2name(hdr, v->rid);

    if (!regidx_overlap(genes, chr, v->pos, v->pos + v->rlen, itr)) { regitr_destroy(itr); return -1; };

    IRegion r = {0};
    r.rid = v->rid;
    while (regitr_overlap(itr)) {
        gf_gene_t *t = regitr_payload(itr, gf_gene_t *);
        if (t->type != GF_PROTEIN_CODING) continue;

        if (r.start == 0) r.start = t->beg;
        if (itr->end > r.end) r.end = t->end;
    }

    regitr_destroy(itr);
    if (r.start == 0 && r.end == 0) return -1; // No protein coding gene was found

    *region = r;

    return 0;
}

int find_gf_tscripts_by_iregion(IRegion iregion, bcf_hdr_t *hdr, gff_t *gff, gf_tscripts **out_gf_tscripts) {
    regidx_t *gff_tscripts = (regidx_t *) gff_get(gff, idx_tscript);
    regitr_t *itr = regitr_init(gff_tscripts);

    const char *chr = bcf_hdr_id2name(hdr, iregion.rid);
    if (!regidx_overlap(gff_tscripts, chr, iregion.start, iregion.end, itr)) { regitr_destroy(itr); return -1; }

    gf_tscripts *tscripts = malloc(sizeof(gf_tscripts));
    kv_init(*tscripts);

    while (regitr_overlap(itr)) {
        kv_push(gf_tscript_t *, *tscripts, regitr_payload(itr, gf_tscript_t *));
    }

    regitr_destroy(itr);
    *out_gf_tscripts = tscripts;

    return 0;
}

void match_gf_tscripts_to_vep_tscripts(gf_tscripts *tscripts, Record *r, gff_t *gff, gf_tscripts **out_gf_scripts) {
    gf_tscripts *record_tscripts = malloc(sizeof(gf_tscripts));
    kv_init(*record_tscripts);

    for (int i = 0; i < r->vep->n; i++) {
        int has_matched = 0;

        if (tscripts == NULL) { kv_push(gf_tscript_t *, *record_tscripts, NULL); continue; }
        for (int j = 0; j < tscripts->n; j++) {
            gf_tscript_t *t = kv_A(*tscripts, j);
            const char *tid = gff_id2string(gff, 0, t->id);

            if (strcmp(tid, kv_A(*r->vep, i).tid) != 0) continue;

            has_matched = 1;
            kv_push(gf_tscript_t *, *record_tscripts, t);
            break;
        }

        if (!has_matched) kv_push(gf_tscript_t *, *record_tscripts, NULL);
    }

    *out_gf_scripts = record_tscripts;
}

int find_gf_tscripts_for_record(Record *r, bcf_hdr_t *hdr, gff_t *gff, gf_tscripts **out_gf_tscripts) {
    IRegion region = variant_to_iregion(r->v);
    gf_tscripts *tscripts = NULL;
    find_gf_tscripts_by_iregion(region, hdr, gff, &tscripts);

    gf_tscripts *matched_tscripts; match_gf_tscripts_to_vep_tscripts(tscripts, r, gff, &matched_tscripts);
    if (tscripts != NULL) { kv_destroy(*tscripts); free(tscripts); }

    *out_gf_tscripts = matched_tscripts;
    return 0;
}

// int find_iregion_transcripts(IRegion *region, bcf_hdr_t *hdr, gff_t *gff, faidx_t *fai, Transcripts **out_transcripts) {
//     regidx_t *tscripts = (regidx_t *) gff_get(gff, idx_tscript);
//     regitr_t *itr = regitr_init(tscripts);
//
//     const char *chr = bcf_hdr_id2name(hdr, region->rid);
//     if (!regidx_overlap(tscripts, chr, region->start, region->end, itr)) return -1;
//
//     Transcripts *transcripts = malloc(sizeof(Transcripts));
//     kv_init(*transcripts);
//     int slen;
//     while (regitr_overlap(itr)) {
//         gf_tscript_t *tr = regitr_payload(itr, gf_tscript_t *);
//         // char *seq = faidx_fetch_seq(fai, chr, itr->beg, itr->end, &slen);
//         // seq[slen] = '\0';
//         // kstring_t s = { slen, slen+1, seq };
//
//         
//
//         Transcript t = { bcf_hdr_name2id(hdr, itr->seq), tr};
//         kv_push(Transcript, *transcripts, t);
//     }
//     regitr_destroy(itr);
//
//     *out_transcripts = transcripts;
//     return 0;
// }

// IRegion transcript_to_iregion(Transcript *t) {
//     return (IRegion) { t->rid, t->definition->beg, t->definition->end };
// }

// void predict_for_seq(Model *models, const char *seq, int seq_len, int strand, float *predictions[], int *num_predictions) {
//     float *ohe; one_hot_encode(seq, seq_len, &ohe);
//
//     if (strand == STRAND_REV) reverse_encoding(ohe, seq_len * ENCODING_SIZE);
//
//     predict(models, seq_len * ENCODING_SIZE, 1, &ohe, num_predictions, predictions);
//     free(ohe);
//
//     if (strand == STRAND_REV) reverse_prediction(*predictions, *num_predictions, 3);
// }

void get_best_insertion_score(float *insert_start, int length, float *best) {
    for (int i = 0; i < length; i+=3) {
        for (int j = 0; j < 3; j++) {
            if (insert_start[i + j] > best[j]) best[j] = insert_start[i+j];
        }
    }
}

void unpack_variants(Records *vs) {
    for (int i = 0; i < vs->n; i++) {
        bcf_unpack(kv_A(*vs, i)->v, BCF_UN_ALL);
    }
}

PhasedRecords *split_records_on_phase(Records *vs, bcf_hdr_t *hdr) {
    PhasedRecords *pvs = malloc(sizeof(PhasedRecords)); kv_init(*pvs); 
    for (int i = 0; i < NUM_PHASES; i++) {
        Records *tmp = malloc(sizeof(Records)); kv_init(*tmp);
        kv_push(Records *, *pvs, tmp);
    }

    for (int i = 0; i < vs->n; i++) {
        Record *r = kv_A(*vs, i);

        int *gt = NULL; int ngts = 0;
        bcf_get_genotypes(hdr, r->v, &gt, &ngts);

        if (gt[0] == gt[1]) { // Homozygous variant, equal opportunity phaser
            kv_push(Record *, *(pvs->a[0]), r);
            kv_push(Record *, *(pvs->a[1]), r);
            free(gt);
            continue;
        }

        // TODO: Add a CLI option to set this behaviour
        if (!bcf_gt_is_phased(gt[0]) && !bcf_gt_is_phased(gt[1])) { // No phasing for heterozygous variants, including on both for now
            kv_push(Record *, *(pvs->a[0]), r);
            kv_push(Record *, *(pvs->a[1]), r);
            log_debug("Adding %s:%li (heterozygous, no phase)", bcf_hdr_id2name(hdr, r->v->rid), r->v->pos); 
            free(gt);
            continue; 
        }

        //Something is phased, check what. Enforced split of multiallelic sites makes this easy
        if (gt[0] > gt[1]) kv_push(Record *, *(pvs->a[0]), r); // gt[0]=4, gt[1]=3, 1|0
        else kv_push(Record *, *(pvs->a[1]), r); // gt[0]=2, gt[1]=5, 0|1
        free(gt);
    }

    return pvs;
}

void remove_overlapping_variants(Records *vs, const bcf_hdr_t *hdr) {
    IRegion last_deletion = {-1, -1, -1};
    int *should_delete = calloc(vs->n, sizeof(int));
    for (int j = 0; j < vs->n; j++) {
        bcf1_t *v = kv_A(*vs, j)->v;

        // Skip variant if it overlaps a previous deletion
        IRegion r = variant_to_iregion(v);
        if (iregion_intersects(last_deletion, r)) {
            log_warn("Skipping on %s:%li, overlaps previous variant", bcf_hdr_id2name(hdr, v->rid), v->pos);
            should_delete[j] = 1;
        }
        if (strlen(v->d.allele[1]) > CAI_DISTANCE) {
            log_warn("Skipping on %s:%li, insertion is too long", bcf_hdr_id2name(hdr, v->rid), v->pos);
            should_delete[j] = 1;
        }
        if (v->rlen > CAI_DISTANCE) {
            log_warn("Skipping on %s:%li, deletion is too long", bcf_hdr_id2name(hdr, v->rid), v->pos);
            should_delete[j] = 1;
        }

        // Store region if deletion
        if (v->rlen > 1) { 
            last_deletion = r;
        }
    }

    int num_variants = vs->n;
    for (int i = 0, vs_i = 0; i < num_variants; i++, vs_i++) {
        if (!should_delete[i]) continue;

        int current_num_variants = vs->n;
        if (vs_i == current_num_variants-1) {
            kv_pop(*vs);
            vs_i--;
            continue;
        }

        for (int j = vs_i; j < vs->n-1; j++) {
            vs->a[j] = vs->a[j+1];
        }
        vs->n--; vs_i--;
    }
}

// int prediction_offset(bcf1_t *v, Variants *vs, int transcript_start, int transcript_end) {
//     int offset = 0;
//     for (int i = 0; i < vs->n; i++) {
//         bcf1_t *vv = kv_A(*vs, i);
//         if (vv->pos + vv->rlen <= v->pos - (CAI_DISTANCE+1) || vv->pos > v->pos + CAI_DISTANCE) continue;
//
//         offset += get_allele_offset_variant(vv, transcript_start, transcript_end);
//     }
//
//     return offset;
// }

int ranges_overlap(Range r1, Range r2) {
    return r1.start <= r2.end && r2.start <= r1.end;
}

int intersect_ranges(Range r1, Range r2, Range *intersection) {
    if (!ranges_overlap(r1, r2)) return 1;

    if (r1.start < r2.start) {
        if (r1.end > r2.end) {
            *intersection = (Range) {r2.start, r2.end};
        } else {
            *intersection = (Range) {r2.start, r1.end};
        }
    } else {
        if (r1.end > r2.end) {
            *intersection = (Range) {r1.start, r2.end};
        } else {
            *intersection = (Range) {r1.start, r1.end};
        }
    }

    return 0;
}

int range_len(Range r) {
    return r.end - r.start;
}

Range scale_range(Range r, int scale) {
    return (Range) {r.start * scale, r.end * scale};
}

void shift_range(Range *r, int value) {
    r->start += value;
    r->end += value;
}

float *match_predictions_to_reference(float *predictions, bcf1_t *v, Variants *vs, Transcript *t, size_t num_reference_predictions) {
    float *corrected_predictions = malloc(num_reference_predictions * sizeof(float));
    Range pred_range = {v->pos - (CAI_DISTANCE+1), v->pos + CAI_DISTANCE};

    int ref_idx = 0, hap_idx = 0;
    for (int l = 0; l < vs->n; l++) {
        bcf1_t *vvv = kv_A(*vs, l);

        Range ref;
        if (intersect_ranges((Range) {vvv->pos, vvv->pos + vvv->rlen}, pred_range, &ref) > 0) continue;
        shift_range(&ref, -pred_range.start);
        int ref_len = range_len(ref);
        if (ref_len == 0) continue;

        Range alt;
        int alt_len = 0;
        if (intersect_ranges((Range) {vvv->pos, vvv->pos + strlen(vvv->d.allele[1])}, pred_range, &alt) == 0) {
            shift_range(&alt, -pred_range.start);
            alt_len = range_len(alt);
        }

        Range ref_segment = scale_range(ref, 3);

        // Copy segments between variants/boundaries
        int segment_len = (ref_segment.start - ref_idx);
        memcpy(corrected_predictions + ref_idx, predictions + hap_idx, segment_len * sizeof(float));

        int offset = ref_len - alt_len; // Naive offset in end position vs ref and alt
        if (offset == 0) {
            memcpy(corrected_predictions + ref_segment.start, predictions + hap_idx + segment_len, range_len(ref_segment) * sizeof(float));
        } else if (offset < 0) {
            float best_scores[3] = {0};
            get_best_insertion_score(predictions + hap_idx + segment_len, range_len(alt), best_scores);
            memcpy(corrected_predictions + ref_segment.start, best_scores, 3 * sizeof(float));
         } else { // offset > 0
            memcpy(corrected_predictions + ref_segment.start, predictions + hap_idx + segment_len, 3 * sizeof(float));
            memset(corrected_predictions + (ref_segment.start + 3), 0, offset * 3 * sizeof(float));
        }

        hap_idx += segment_len + alt_len*3;
        ref_idx = ref_segment.end;
    }

    memcpy(corrected_predictions + ref_idx, predictions + hap_idx, (num_reference_predictions - ref_idx) * sizeof(float));

    return corrected_predictions;
}

// int64_t get_position_in_spliced_haplotype(bcf1_t *v, Variants *vs, Transcript *t) {
//     int positions_skipped = 0;
//     int previous_exon_end = t->definition->beg;
//
//     for (int i = 0; i < t->definition->ncds; i++) {
//         int64_t start = t->definition->cds[i]->beg;
//         if (start > v->pos) break;
//         int64_t end = start + t->definition->cds[i]->len;
//
//         int64_t start_in_hap = get_position_in_haplotype(start, vs, t->definition->beg, t->definition->end);
//         int64_t end_in_hap = get_position_in_haplotype(end, vs, t->definition->beg, t->definition->end);
//
//         positions_skipped += start_in_hap - previous_exon_end; // Size of intron between current and previous exon
//
//         if (iregion_intersects((IRegion) {-1, v->pos, v->pos + strlen(v->d.allele[1])}, (IRegion) {-1, start, end})) {
//             return get_position_in_haplotype(v->pos, vs, t->definition->beg, t->definition->end) - t->definition->beg - positions_skipped;
//         }
//
//         previous_exon_end = end_in_hap;
//     }
//
//     return -1;
// }

// void process_phased_variants(Records *vs, HdrFormat *vep_fmt, gff_t *gff, faidx_t *fai, bcf_hdr_t *hdr, Model *models) {
//     kvec_t(int64_t) frameshifts_rescued; kv_init(frameshifts_rescued);
//
//
//     for (int j = 0; j < vs->n; j++) {
//         Record *r = kv_A(*vs, j);
//         bcf1_t *vv = kv_A(*vs, j)->v;
//         VepTranscripts *transcript_csqs = get_vep_transcript_consequences(vv, vep_fmt, hdr);
//
//         // Filter transcripts for overlap with current variant
//         Region variant_region = iregion_to_region(variant_to_iregion(vv), hdr);
//         Transcripts *variant_transcripts; find_region_transcripts(&variant_region, hdr, gff, fai, &variant_transcripts);
//
//         Transcript *v_ts[transcript_csqs->n];
//         Variants *v_vs[transcript_csqs->n];
//         Variants *v_vs2[transcript_csqs->n];
//         for (int k = 0; k < transcript_csqs->n; k++) {
//             v_ts[k] = NULL;
//             v_vs[k] = NULL;
//             v_vs2[k] = NULL;
//         }
//
//         for (int k = 0; k < variant_transcripts->n; k++) {
//             Transcript t = kv_A(*variant_transcripts, k);
//             const char *tid = gff_id2string(gff, 0, t.definition->id);
//
//             for (int l = 0; l < transcript_csqs->n; l++) {
//                 if (strcmp(tid, kv_A(*transcript_csqs, l).tr_id) == 0) {
//                     v_ts[l] = &variant_transcripts->a[k];
//
//                     // Find region for matched transcript
//                     IRegion tr_region = transcript_to_iregion(&t);
//
//                     // Create subset of variants overlapping the current transcript, excluding the current variant
//                     Variants *vs_j = malloc(sizeof(Variants)); kv_init(*vs_j);
//                     Variants *vs_j2 = malloc(sizeof(Variants)); kv_init(*vs_j2);
//                     for (int m = 0; m < vs->n; m++) {
//
//                         bcf1_t *t_v = kv_A(*vs, m)->v;
//                         if (iregion_intersects(tr_region, variant_to_iregion(t_v))) {
//                             kv_push(bcf1_t *, *vs_j2, t_v);
//                             if (m == j) continue;
//                             kv_push(bcf1_t *, *vs_j, t_v);
//                         }
//                     }
//
//                     /* if (vs_j->n > 0)*/ v_vs[l] = vs_j;
//                     /*if (vs_j->n > 0)*/ v_vs2[l] = vs_j2;
//                     break;
//                 }
//             }
//         }
//
//         // ConsequencePredictions *csq_preds; init_consequence_predictions(csq_preds, transcript_csqs->n);
//         for (int k = 0; k < transcript_csqs->n; k++) {
//             if (v_ts[k] == NULL) continue; // There is no coding transcript overlapping the current variant;
//             if (v_vs[k] == NULL) continue; // There is no difference between the haplotype and reference 
//
//             Transcript *t = v_ts[k];
//             Variants *vv_vs = v_vs[k];
//             Variants *vv_vs2 = v_vs2[k];
//
//             // Calculate offsets between ref, hap_ref and hap_alt
//             int hap_ref_offset = 0;
//             for (int l = 0; l < vv_vs->n; l++) hap_ref_offset += get_allele_offset_variant(kv_A(*vv_vs, l), t->definition->beg, t->definition->end);
//
//             kstring_t hap_ref; get_haplotype_sequence(t->seq, *vv_vs, hap_ref_offset, t->definition->beg, &hap_ref);
//
//             int alt_allele_offset = get_allele_offset_variant(vv, t->definition->beg, t->definition->end);
//             int hap_alt_offset = hap_ref_offset + alt_allele_offset;
//             kstring_t hap_alt; get_haplotype_sequence(t->seq, *vv_vs2, hap_alt_offset, t->definition->beg, &hap_alt);
//
//
//             // Check vep consequences for triggers for rechecking consequences in haplotype context
//             int trigger = 0;
//             for (int l = 0; l < transcript_csqs->a[k].num_consequences; l++) {
//                 const char *consequence = transcript_csqs->a[k].consequences[l];
//                 if ((strcmp(consequence, SPLICE_DONOR_VARIANT) == 0) ||
//                     (strcmp(consequence, SPLICE_ACCEPTOR_VARIANT) == 0)) {
//                     trigger = SPLICE_TRIGGER;
//                     // log_warn("Splice trigger %s:%li = %s", bcf_hdr_id2name(hdr, vv->rid), vv->pos, consequence);
//                 } else if (strcmp(consequence, FRAMESHIFT_VARIANT) == 0) {
//                     trigger = FRAMESHIFT_TRIGGER;
//                     // log_warn("Frameshift trigger %s:%li = %s", bcf_hdr_id2name(hdr, vv->rid), vv->pos, consequence);
//                 } else if (strcmp(consequence, STOP_GAINED_VARIANT) == 0) {
//                     trigger = STOP_TRIGGER;
//                 }
//             }
//
//             int is_rescued = 0;
//             if (trigger == SPLICE_TRIGGER) {
//
//                 // Run reference
//                 int32_t start = VARIANT_POS + (t->definition->beg - vv->pos);
//                 int32_t end = VARIANT_POS - (t->definition->end - (vv->pos+1)); // Open ended, so + 1
//                 char *ref_seq = t->seq.s - start;
//
//                 if (start < 0) start = 0;
//                 if (end < 0) end = 0;
//
//                 Range ref_boundary = { start, end };
//                 char *ref = pad_sequence(ref_seq, ref_boundary, WIDTH);
//                 float *ref_predictions; int num_ref_predictions;
//                 int ref_len = strlen(ref);
//                 predict_for_seq(models, ref, ref_len, t->definition->strand, &ref_predictions, &num_ref_predictions);
//                 free(ref);
//
//                 // Run haplotype reference
//                 int pred_offset = prediction_offset(vv, vv_vs, t->definition->beg, t->definition->end);
//                 int64_t transcript_start = get_position_in_haplotype(t->definition->beg, vv_vs, t->definition->beg, t->definition->end);
//                 int64_t transcript_end = get_position_in_haplotype(t->definition->end, vv_vs, t->definition->beg, t->definition->end);
//                 int64_t vposition = get_position_in_haplotype(vv->pos, vv_vs, t->definition->beg, t->definition->end);
//                 int32_t hap_ref_start = VARIANT_POS + (transcript_start - vposition);
//                 int32_t hap_ref_end = VARIANT_POS - (transcript_end - vposition) + pred_offset;
//                 char *hap_ref_seq = hap_ref.s - hap_ref_start;
//                 int hap_ref_seq_len = strlen(hap_ref_seq);
//                 if (hap_ref_start < 0) hap_ref_start = 0;
//                 if (hap_ref_end < 0) hap_ref_end = 0;
//                 Range hap_ref_boundary = {hap_ref_start, hap_ref_end};
//                 char *hap_ref_padded = pad_sequence(hap_ref_seq, hap_ref_boundary, WIDTH + pred_offset);
//                 float *hap_ref_predictions; int num_hap_ref_predictions;
//                 int hap_ref_len = strlen(hap_ref_padded);
//                 predict_for_seq(models, hap_ref_padded, hap_ref_len, t->definition->strand, &hap_ref_predictions, &num_hap_ref_predictions);
//                 float *corrected_hap_ref_predictions = match_predictions_to_reference(hap_ref_predictions, vv, vv_vs, t, num_ref_predictions);
//                 Score score = calculate_delta_scores(vv->d.allele[1], t->definition->gene->name, ref_predictions, corrected_hap_ref_predictions, num_ref_predictions, CAI_DISTANCE);
//
//
//                 // // Run haplotype alt
//                 transcript_end += get_allele_offset_variant(vv, t->definition->beg, t->definition->end);
//                 int hap_alt_pred_offset = prediction_offset(vv, vv_vs2, t->definition->beg, t->definition->end);
//                 int32_t hap_alt_start = VARIANT_POS + (transcript_start - vposition);
//                 int32_t hap_alt_end = VARIANT_POS - (transcript_end - vposition) + hap_alt_pred_offset;
//                 char *hap_alt_seq = hap_alt.s - hap_alt_start;
//                 int hap_alt_seq_len = strlen(hap_alt_seq);
//                 if (hap_alt_start < 0) hap_alt_start = 0;
//                 if (hap_alt_end < 0) hap_alt_end = 0;
//                 Range hap_alt_boundary = {hap_alt_start, hap_alt_end};
//                 char *hap_alt_padded = pad_sequence(hap_alt_seq, hap_alt_boundary, WIDTH + hap_alt_pred_offset); //+ get_allele_offset(vv, t->definition->beg, t->definition->end));
//                 float *hap_alt_predictions; int num_hap_alt_predictions;
//                 int hap_alt_len = strlen(hap_alt_padded);
//                 predict_for_seq(models, hap_alt_padded, hap_alt_len, t->definition->strand, &hap_alt_predictions, &num_hap_alt_predictions);
//                 float *corrected_hap_alt_predictions = match_predictions_to_reference(hap_alt_predictions, vv, vv_vs2, t, num_ref_predictions);
//                 Score score_alt = calculate_delta_scores(vv->d.allele[1], t->definition->gene->name, corrected_hap_ref_predictions, corrected_hap_alt_predictions, num_ref_predictions, CAI_DISTANCE);
//
//                 // Calculate scores between ref and hap_alt
//                 Score score_hap = calculate_delta_scores(vv->d.allele[1], t->definition->gene->name, ref_predictions, corrected_hap_alt_predictions, num_ref_predictions, CAI_DISTANCE);
//
//
//                 if (r->c->hap_score[k]->l > 0) {
//                     kputc('&', r->c->hap_score[k]);
//                 }
//
//                 if (score_hap.dl > SPLICEAI_THRESHOLD && score_hap.dg > SPLICEAI_THRESHOLD && (score_hap.dg_idx - score_hap.dl_idx) % 3 == 0) is_rescued = 1;
//
//                 if (score_hap.al > SPLICEAI_THRESHOLD && score_hap.ag > SPLICEAI_THRESHOLD && (score_hap.ag_idx - score_hap.al_idx) % 3 == 0) is_rescued = 1;
//
//                 if (score_hap.al < SPLICEAI_RESCUE_THRESHOLD && score_hap.dl < SPLICEAI_RESCUE_THRESHOLD) is_rescued = 1;
//
//                 if (is_rescued) {
//                     log_info("%s:%li is rescued! Original consequence: %s", bcf_hdr_id2name(hdr, vv->rid), vv->pos, transcript_csqs->a[k].consequences[0]);
//                 }
//
//                 ksprintf(r->c->hap_score[k], "%s|%s|%.2f|%.2f|%.2f|%.2f|%d|%d|%d|%d", score_hap.alt, score_hap.gene, score_hap.ag, score_hap.al, score_hap.dg, score_hap.dl, score_hap.ag_idx, score_hap.al_idx, score_hap.dg_idx, score_hap.dl_idx);
//                 log_info(r->c->hap_score[k]->s);
//
//             free(ref_predictions);
//             free(hap_ref_predictions);
//             free(corrected_hap_ref_predictions);
//             free(hap_alt_predictions);
//             free(corrected_hap_alt_predictions);
//
//             }
//
//             if (trigger) {
//
//                 kstring_t spliced_transcript = {0};
//                 // Splice transcript
//                 int length = 0;
//                 int hap_length = 0;
//                 for (int l = 0; l < t->definition->ncds; l++) {
//                     int64_t start = t->definition->cds[l]->beg;
//                     int64_t end = start + t->definition->cds[l]->len;
//                     int64_t start_in_hap = get_position_in_haplotype(start, vv_vs2, t->definition->beg, t->definition->end);
//                     int64_t end_in_hap = get_position_in_haplotype(end, vv_vs2, t->definition->beg, t->definition->end);
//
//                     length += (end - start);
//                     hap_length += (end_in_hap - start_in_hap);
//                     kputsn(hap_alt.s + (start_in_hap - t->definition->beg), end_in_hap - start_in_hap, &spliced_transcript);
//                 }
//
//                 int v_pos = get_position_in_spliced_haplotype(vv, vv_vs2, t);
//                 if (trigger == FRAMESHIFT_TRIGGER) {
//                     int offset = get_allele_offset_variant(vv, t->definition->beg, t->definition->end);
//                     if (v_pos == -1) continue;
//
//                     kvec_t(int64_t) frameshifts_used; kv_init(frameshifts_used);
//                     kv_push(int64_t, frameshifts_used, vv->pos);
//
//                     if (t->definition->strand == STRAND_FWD) {
//                         for (int l = 0; l < vv_vs->n; l++) {
//                             bcf1_t *vvv = kv_A(*vv_vs, l);
//
//                             int already_rescued = 0;
//                             for (int m = 0; m < frameshifts_rescued.n; m++) {
//                                 if (kv_A(frameshifts_rescued, m) == vvv->pos) {
//                                     already_rescued = 1;
//                                     break;
//                                 }
//                             }
//                             if (already_rescued) continue;
//
//                             int v2_pos = get_position_in_spliced_haplotype(vvv, vv_vs2, t);
//                             int variant_offset = get_allele_offset_variant(vvv, t->definition->beg, t->definition->end);
//                             if (variant_offset != 0 && v2_pos != -1 && v2_pos > v_pos && v2_pos <= v_pos + MAX_FRAMESHIFT_RADIUS) {
//                                 offset += variant_offset;
//                                 kv_push(int64_t, frameshifts_used, vvv->pos);
//                             }
//                             if (v2_pos > v_pos + MAX_FRAMESHIFT_RADIUS) break;
//                         }
//                     }
//                     else if (t->definition->strand == STRAND_REV) {
//                         for (int l = vv_vs->n-1; l >= 0; l--) {
//                             bcf1_t *vvv = kv_A(*vv_vs, l);
//
//                             int already_rescued = 0;
//                             for (int m = 0; m < frameshifts_rescued.n; m++) {
//                                 if (kv_A(frameshifts_rescued, m) == vvv->pos) {
//                                     already_rescued = 1;
//                                     break;
//                                 }
//                             }
//                             if (already_rescued) continue;
//
//                             int v2_pos = get_position_in_spliced_haplotype(vvv, vv_vs2, t);
//                             int variant_offset = get_allele_offset_variant(vvv, t->definition->beg, t->definition->end);
//                             if (variant_offset != 0 && v2_pos != -1 && v2_pos < v_pos && v2_pos >= v_pos - MAX_FRAMESHIFT_RADIUS) {
//                                 offset += variant_offset;
//                                 kv_push(int64_t, frameshifts_used, vvv->pos);
//                             }
//                             if (v2_pos < v_pos - MAX_FRAMESHIFT_RADIUS) break;
//                         }
//                     }
//
//                     if (offset % 3 == 0) {
//                         log_info("%s:%li is rescued! Original consequence: %s", bcf_hdr_id2name(hdr, vv->rid), vv->pos, transcript_csqs->a[k].consequences[0]);
//
//                         is_rescued = 1;
//
//                         for (int l = 0; l < frameshifts_used.n; l++) {
//                             kv_push(int64_t, frameshifts_rescued, kv_A(frameshifts_used, l));
//                         }
//                     }
//
//                     kv_destroy(frameshifts_used);
//                 }
//
//                 if (trigger == STOP_TRIGGER) {
//                     int offset = get_allele_offset_variant(vv, t->definition->beg, t->definition->end);
//                     int start = v_pos / 3;
//                     int end = offset > 0 ? (v_pos + offset) / 3 + 1 : v_pos / 3 + 1;
//                     int stop_hit = 0;
//                     if (t->definition->strand == STRAND_FWD) {
//                         v_pos = get_position_in_spliced_haplotype(vv, vv_vs2, t);
//                         const char *stop_codons[3] = STOP_CODONS;
//                         for (int l = start; l < end; l++) {
//                             for (int m = 0; m < 3; m++) {
//                                 const char *stop = stop_codons[m];
//                                 char *frame = spliced_transcript.s + l * 3;
//                                 if (stop[0] == frame[0] && stop[1] == frame[1] && stop[2] == frame[2]) {
//                                     stop_hit = 1;
//                                     break;
//                                 }
//                             }
//                             if (stop_hit) break;
//                         }
//                     } else {
//                         stop_hit = 1;
//                     }
//                     // else if (t->definition->strand == STRAND_REV) {
//                     //     const char *stop_codons[3] = COMPLEMENTARY_STOP_CODONS;
//                     //     for (int l = end-1; l >= start; l--) {
//                     //         for (int m = 0; m < 3; m++) {
//                     //             const char *stop = stop_codons[m];
//                     //             char *frame = spliced_transcript.s + l * 3;
//                     //             if (stop[0] == frame[2] && stop[1] == frame[1] && stop[2] == frame[0]) {
//                     //                 stop_hit = 1;
//                     //                 break;
//                     //             }
//                     //         }
//                     //
//                     //         if (stop_hit) break;
//                     //     }
//                     // }
//
//                     if (!stop_hit) {
//                         log_info("%s:%li is rescued! Original consequence: %s", bcf_hdr_id2name(hdr, vv->rid), vv->pos, transcript_csqs->a[k].consequences[0]);
//                         is_rescued = 1;
//                     }
//                 }
//
//
//                 free(spliced_transcript.s);
//             }
//
//             if (r->c->is_rescue[k]->l > 0) {
//                 kputc('&', r->c->is_rescue[k]);
//             }
//             kputw(is_rescued, r->c->is_rescue[k]);
//
//             if (r->c->rescue_type[k]->l > 0) {
//                 kputc('&', r->c->rescue_type[k]);
//             }
//             if (is_rescued) {
//                 if (trigger == SPLICE_TRIGGER) {
//                     kputs("rescued_splice", r->c->rescue_type[k]);
//                 } else if (trigger == FRAMESHIFT_TRIGGER) {
//                     kputs("rescued_frameshift", r->c->rescue_type[k]);
//                 } else {
//                     kputs("rescued_stop", r->c->rescue_type[k]);
//                 }
//             } else {
//                 kputc('.', r->c->rescue_type[k]);
//             }
//
//             if (trigger != SPLICE_TRIGGER) {
//                 if (r->c->hap_score[k]->l > 0) kputc('&', r->c->hap_score[k]);
//                 kputc('.', r->c->hap_score[k]);
//             }
//
//             // Free per transcript predicted for
//             free(hap_ref.s);
//             free(hap_alt.s);
//             kv_destroy(*vv_vs);
//             kv_destroy(*vv_vs2);
//             free(vv_vs);
//             free(vv_vs2);
//         }
//
//         // Free per variant predicted for
//         for (int i = 0; i < variant_transcripts->n; i++) { free(kv_A(*variant_transcripts, i).seq.s); }
//         kv_destroy(*variant_transcripts); free(variant_transcripts);
//
//         for (int i = 0; i < transcript_csqs->n; i++) { 
//             VepTranscript tcs = kv_A(*transcript_csqs, i);
//             free(tcs.tr_id);
//             for (int j = 0; j < transcript_csqs->a[i].num_consequences; j++) {
//                 free(transcript_csqs->a[i].consequences[j]);
//             }
//             free(tcs.consequences);
//         }
//         kv_destroy(*transcript_csqs);
//         free(transcript_csqs);
//
//     }
//
//     kv_destroy(frameshifts_rescued);
// }

int add_rescueranger_tags(bcf_hdr_t *hdr) {
    const char *tags[7] = INFO_TAGS;
    const char *descriptions[7] = { SPLICEAI_HAP_DESC, CS_DESC, RESCUE_DESC, RESCUE_TYPE_DESC, PTV_DESC, PTV_TYPE_DESC, PREDICTION_FLAGS_DESC };

    for (int i = 0; i < 7; i++) {
        if (bcf_hdr_append(hdr, descriptions[i]) != 0) {
            log_error("Failed to append description for tag %s to vcf header.", tags[i]);
            return -1;
        }
    }

    return 0;
}

int main(int argc, char *argv[]) {
    log_set_level(2);
    setenv("TF_CPP_MIN_LOG_LEVEL", "1", 1);

    // log_debug("Loading SpliceAI models...");
    // Model *models = load_models(argv[1]);

    // Parse gff file
    const char *gff_path = argv[2];
    log_debug("Loading GFF from: %s", gff_path);
    gff_t *gff = gff_init(gff_path);
    gff_parse(gff);

    // Open fasta file
    const char *fasta_path = argv[3];
    log_debug("Loading FASTA from: %s", fasta_path);
    faidx_t *fai = fai_load(fasta_path);

    // Open VCF
    const char *input_vcf_path = argv[4];
    log_debug("Loading VCF from: %s", input_vcf_path);
    htsFile *fp_in = hts_open(input_vcf_path, "r");
    bcf_hdr_t *hdr = bcf_hdr_read(fp_in);
    HdrFormat *vep_fmt = parse_csq_format(hdr);

    const char *output_vcf_path = argv[5];
    htsFile *fp_out = hts_open(output_vcf_path, "w");
    add_rescueranger_tags(hdr);

    bcf_hdr_write(fp_out, hdr);

    bcf1_t *v = bcf_init();

    kvec_t(RecordBlock) batches; kv_init(batches);

    Records *batch = malloc(sizeof(Records));
    kv_init(*batch);

    IRegion batch_region = {-1, -1, -1};

    while (bcf_read(fp_in, hdr, v) >= 0) {
        int end_reached = 0;

        IRegion var_region = {0};
        int is_non_coding = find_gene_region_for_variant(v, hdr, gff, &var_region);
        if (batch_region.rid == -1) batch_region = var_region;
        int overlaps_batch = is_non_coding ? 0 : iregion_intersects(batch_region, var_region);

        if (batch->n > 0 && !overlaps_batch) { // Variant is either non-coding or falls outside of batch
            RecordBlock block = {kv_A(*batch, 0)->vep == NULL ? 0 : 1, batch};
            kv_push(RecordBlock, batches, block);

            batch_region = var_region;
            batch = malloc(sizeof(Records));
            kv_init(*batch);
        }

        VepCsqs *vep_csqs = is_non_coding ? NULL : parse_vep_transcript_consequences(v, vep_fmt, hdr);
        if (overlaps_batch && var_region.end > batch_region.end) batch_region.end = var_region.end;
        kv_push(Record *, *batch, init_record(v, vep_csqs));
    }

    if (batch->n > 0) {
        RecordBlock block = {kv_A(*batch, 0)->vep == NULL ? 0 : 1, batch};
        kv_push(RecordBlock, batches, block);
    }

    Genes genes; kv_init(genes); // TODO: Turn this into a hash for faster accessing.
    for (int b = 0; b < batches.n; b++) {
        RecordBlock block = kv_A(batches, b);
        if (!block.is_coding) continue;

        // Find matching gff transcripts/genes for records in batch
        for (int i = 0; i < block.records->n; i++) {
            Record *r = kv_A(*block.records, i);
            find_gf_tscripts_for_record(r, hdr, gff, &r->tscripts);
            if (r->tscripts == NULL) continue; // A variant in a coding gene can still have no overlapping coding transcripts

            r->genes = malloc(sizeof(Genes));
            kv_init(*r->genes);
            for (int j = 0; j < r->tscripts->n; j++) {
                gf_tscript_t *t = kv_A(*r->tscripts, j);
                if (t == NULL) { kv_push(Gene *, *r->genes, NULL); continue; }

                uint32_t gene_id = t->gene->id;
                int has_matched = 0;

                for (int k = 0; k < genes.n; k++) {
                    Gene *gene = kv_A(genes, k);
                    if (gene_id != gene->gff_id) continue;

                    has_matched = 1;
                    kv_push(Gene *, *r->genes, gene);
                }

                if (has_matched) continue;

                // No match found, create and store a new gene
                Range gene_range = {t->gene->beg, t->gene->end};

                int slen = 0;
                char *seq = faidx_fetch_seq(fai, bcf_hdr_id2name(hdr, r->v->rid), gene_range.start, gene_range.end, &slen);
                seq[slen] = '\0';
                kstring_t s = { slen, slen+1, seq };

                Gene *gene = malloc(sizeof(Gene));
                gene->gff_id = gene_id;
                gene->coords = gene_range;
                gene->seq = s;
                gene->strand = t->gene->strand;
                kv_push(Gene *, genes, gene);
            }
        }
    }

    bcf_destroy(v);
    destroy_hdr_fmt(vep_fmt);
    bcf_close(fp_in);
    fai_destroy(fai);

    for (int b = 0; b < batches.n; b++) {
        RecordBlock block = kv_A(batches, b);
        if (!block.is_coding) continue;

        for (int i = 0; i < batch->n; i++) {
            Record *r = kv_A(*block.records, i);

            // Check VEP consequences for rescue triggers
            for (int j = 0; j < r->vep->n; j++) {
                TscriptCsqs transcript_consequence = kv_A(*r->vep, j);
                for (int k = 0; k < transcript_consequence.num_consequences; k++) {
                    const char *consequence = transcript_consequence.consequences[k];

                    if (strcmp(consequence, SPLICE_DONOR_VARIANT) == 0 ||
                        strcmp(consequence, SPLICE_ACCEPTOR_VARIANT) == 0) {
                        r->triggers[SPLICE_TRIGGER] = 1;
                    } else if (strcmp(consequence, FRAMESHIFT_VARIANT) == 0) {
                        r->triggers[FRAMESHIFT_TRIGGER] = 1;
                    } else if (strcmp(consequence, STOP_GAINED_VARIANT) == 0) {
                        r->triggers[STOP_TRIGGER] = 1;
                    }
                }
            }

            // TODO: Check SpliceAI output for rescue triggers
        }
    }

    // Split records by phase
    for (int b = 0; b < batches.n; b++) {
        RecordBlock block = kv_A(batches, b);
        if (!block.is_coding) continue;

        PhasedRecords *phased_records = split_records_on_phase(block.records, hdr);
        for (int p = 0; p < NUM_PHASES; p++) {
            Records *batch = kv_A(*phased_records, p);

            for (int i = 0; i < batch->n; i++) {
                Record *r = kv_A(*batch, i);
                if (!r->triggers[SPLICE_TRIGGER] && !r->triggers[FRAMESHIFT_TRIGGER] && !r->triggers[STOP_TRIGGER]) continue;

                for (int g = 0; g < r->genes->n; g++) {
                    Gene *gene = kv_A(*r->genes, g);

                    // Find records for the current allele
                    Records hap_in; kv_init(hap_in); // Variants for current allele including current variant
                    Records hap_ex; kv_init(hap_ex); // Variant for current allele excluding current variant

                    for (int j = 0; j < batch->n; j++) {
                        Record *hap_r = kv_A(*batch, j);

                        for (int g2 = 0; g2 < hap_r->genes->n; g2++) {
                            Gene *gene2 = kv_A(*hap_r->genes, g2);
                            if (gene->gff_id == gene2->gff_id) {
                                kv_push(Record *, hap_in, hap_r);

                                if (hap_r == r) continue;
                                kv_push(Record *, hap_ex, hap_r);
                            }
                        }
                    }

                    // Predict splicing



                    kv_destroy(hap_in);
                    kv_destroy(hap_ex);
                }
            }

            kv_destroy(*batch);
            free(batch);
        }

        kv_destroy(*phased_records);
        free(phased_records);
    }


    // Free genes
    for (int i = 0; i < genes.n; i++) {
        free(kv_A(genes, i)->seq.s);
        free(kv_A(genes, i));
    }
    kv_destroy(genes);

    // Free batches
    for (int i = 0; i < batches.n; i++) {
        RecordBlock block = kv_A(batches, i);

        for (int j = 0; j < block.records->n; j++) {
            Record *r = kv_A(*block.records, j);
            destroy_record(r);
        }
        kv_destroy(*block.records);
        free(block.records);
    }
    kv_destroy(batches);

    // while (bcf_read(fp_in, hdr, v) >= 0) {
    //     int is_non_coding = find_variant_transcript_region(v, hdr, gff, &var_region);
    //     if (is_non_coding) {
    //         if (batch.n > 0) { // There exists a batch
    //             // Process batch
    //             // Clear batch
    //         }
    //
    //         bcf_write(fp_out, hdr, v);
    //         continue; // Current variant is non-coding, go to next iteration
    //     }
    //
    //     // Check for overlap between current variant and batch regions
    //     if (iregion_intersects(batch_region, var_region) == 1) { // Variant overlaps with batch region, append to current batch
    //         if (var_region.end > batch_region.end) batch_region.end = var_region.end; // Update overlap region if necessary
    //
    //         // Convert variant to record, add to batch
    //         VepCsqs *csqs = parse_vep_transcript_consequences(v, vep_fmt, hdr);
    //         kv_push(Record *, batch, init_record(v, csqs));
    //
    //         continue; // Saved variant to batch, go to next iteration
    //     }
    //
    //     /* Batch completed, process it */
    //
    //     // Find matching gff transcripts/genes for records in batch
    //     kvec_t(Gene) genes; kv_init(genes);
    //     for (int i = 0; i < batch.n; i++) {
    //         Record *r = kv_A(batch, i);
    //         find_gf_tscripts_for_record(r, hdr, gff, &r->tscripts);
    //
    //         int slen = 0;
    //         for (int j = 0; j < r->tscripts->n; j++) {
    //             gf_tscript_t t = kv_A(*r->tscripts, j);
    //             uint32_t gene_id = t.gene->id;
    //
    //             int has_matched = 0;
    //             for (int k = 0; k < genes.n; k++) {
    //                 Gene gene = kv_A(genes, k);
    //                 if (gene_id != gene.gff_id) continue;
    //
    //                 has_matched = 1;
    //                 r->gene = &gene;
    //             }
    //
    //             if (has_matched) continue;
    //
    //             // No match found, create and store a new gene
    //             Range gene_range = {t.gene->beg, t.gene->end};
    //
    //             char *seq = faidx_fetch_seq(fai, bcf_hdr_id2name(hdr, r->v->rid), gene_range.start, gene_range.end, &slen);
    //             seq[slen] = '\0';
    //             kstring_t s = { slen, slen+1, seq };
    //
    //             Gene gene = { gene_id, gene_range, s, t.gene->strand };
    //             kv_push(Gene, genes, gene);
    //         }
    //     }
    //
    //     // Phase records
    //     PhasedRecords *phased_records = group_variants_by_phase(&batch, hdr);
    //     for (int i = 0; i < NUM_PHASES; i++) {
    //         Records *records = kv_A(*phased_records, i);
    //
    //         // Preprocess records
    //         if (records->n < 1) continue;
    //         remove_overlapping_variants(records, hdr);
    //
    //         if (records->n < 1) continue;
    //
    //         // Create inclusive and exclusive haplotype record sets
    //         // Calculate haplotype sequence
    //
    //         // process_phased_variants(ppvs, vep_fmt, gff, fai, hdr, models);
    //     }
    //
    //     for (int i = 0; i < genes.n; i++) {
    //         Gene gene = kv_A(genes, i);
    //         free(gene.seq.s);
    //     }
    //     kv_destroy(genes);
    //
    //     for (int i = 0; i < phased_records->n; i++) {
    //         kv_destroy(*kv_A(*phased_records, i));
    //         free(kv_A(*phased_records, i));
    //     };
    //     kv_destroy(*phased_records); free(phased_records);
    //
    //     // Write out and clear current batch
    //     // for (int i = 0; i < batch.n; i++) {
    //     //     kstring_t hap_score = {0}, is_rescue = {0}, rescue_type = {0};
    //     //     Record *r = kv_A(batch, i);
    //     //     for (int j = 0; j < r->c->num_consequences; j++) {
    //     //         if (j != 0) {
    //     //             kputc(',', &hap_score);
    //     //             kputc(',', &is_rescue);
    //     //             kputc(',', &rescue_type);
    //     //         }
    //     //
    //     //         if (r->c->hap_score[j]->l > 0) kputs(r->c->hap_score[j]->s, &hap_score);
    //     //         else kputc('.', &hap_score);
    //     //         if (r->c->is_rescue[j]->l > 0) kputs(r->c->is_rescue[j]->s, &is_rescue);
    //     //         else kputc('.', &is_rescue);
    //     //         if (r->c->rescue_type[j]->l > 0) kputs(r->c->rescue_type[j]->s, &rescue_type);
    //     //         else kputc('.', &rescue_type);
    //     //     }
    //     //
    //     //     bcf_update_info_string(hdr, r->v, SPLICEAI_HAP_TAG, hap_score.s);
    //     //     bcf_update_info_string(hdr, r->v, RESCUE_TAG, is_rescue.s);
    //     //     bcf_update_info_string(hdr, r->v, RESCUE_TYPE_TAG, rescue_type.s);
    //     //
    //     //
    //     //     bcf_write(fp_out, hdr, r->v);
    //     //     bcf_destroy(r->v);
    //     //
    //     //     // destroy_consequence_predictions(r->c);
    //     //     // free(r);
    //     //     // free(hap_score.s);
    //     // }
    //     for (int i = 0; i < batch.n; i++) {
    //         Record *r = kv_A(batch, i);
    //         destroy_record(r);
    //     }
    //     kv_destroy(batch);
    //
    //     // Start new batch, with current variant as the first batch entry
    //     batch_region = var_region;
    //     kv_init(batch);
    //     Record *r = init_record(bcf_dup(v), parse_vep_transcript_consequences(v, vep_fmt, hdr));
    //     kv_push(Record *, batch, r);
    // }

    // Process final batch
    // unpack_variants(&batch);
    // PhasedRecords *pvs = group_variants_by_phase(&batch, hdr);
    // for (int i = 0; i < NUM_PHASES; i++) {
    //     Records *ppvs = kv_A(*pvs, i);
    //     if (ppvs->n == 0) continue;
    //     remove_overlapping_variants(ppvs, hdr);
    //
    //     if (ppvs->n == 9 || ppvs->n == -1) continue;
    //     // process_phased_variants(ppvs, vep_fmt, gff, fai, hdr, models);
    // }

    // for (int i = 0; i < pvs->n; i++) { kv_destroy(*kv_A(*pvs, i)); free(kv_A(*pvs, i)); };
    // kv_destroy(*pvs); free(pvs);

    // Clear out final batch
    // for (int i = 0; i < batch.n; i++) {
    //     bcf_write(fp_out, hdr, kv_A(batch, i)->v);
    //     bcf_destroy(kv_A(batch, i)->v);
    // }
    // kv_destroy(batch);

    // Free remaining allocs

    gff_destroy(gff);
    bcf_hdr_destroy(hdr);
    bcf_close(fp_out);


    // Destroy SpliceAI models
    // destroy_models(models);
}

