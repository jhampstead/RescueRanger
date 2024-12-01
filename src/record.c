#include "record.h"
#include <stdlib.h>

const int n_csq_fields = 29;
const char vcf_delim = ',';
const char *csq_field = "CSQ";
const char *hdr_desc = "Description";
const char vep_delim = '|';
const char *vep_desc_prefix = "Format: ";
const int vep_desc_prefix_len = 8;
const char consequence_delim = '&';

const char *vep_consequence = "Consequence";
const char *vep_feature_type = "Feature_type";
const char *vep_feature = "Feature";
const char *vep_feature_type_transcript = "Transcript";


void init_splice_site(SpliceSite **site, size_t *n, size_t *m) {
    *m = 1;
    *n = 0;
    *site = malloc(sizeof(SpliceSite));
}

void resize_splice_site(SpliceSite **site, size_t *m) {
    *m *= 2;
    *site = realloc(*site, *m * sizeof(SpliceSite));
}

void destroy_splice_site(SpliceSite *site, size_t n) {
    free(site);
}


void init_consequence_predictions(ConsequencePredictions **dst, int num_consequences) {
    ConsequencePredictions *csq_preds = malloc(sizeof(ConsequencePredictions));
    csq_preds->num_consequences = num_consequences;
    // csq_preds->hap_score = malloc(num_consequences * sizeof(kstring_t *));
    csq_preds->splice_change = malloc(num_consequences * sizeof(kstring_t *));
    csq_preds->is_rescue = malloc(num_consequences * sizeof(kstring_t *));
    csq_preds->rescue_type = malloc(num_consequences * sizeof(kstring_t *));
    csq_preds->is_ptv = malloc(num_consequences * sizeof(kstring_t *));
    csq_preds->ptv_type = malloc(num_consequences * sizeof(kstring_t *));
    csq_preds->flags = malloc(num_consequences * sizeof(kstring_t *));

    for (int i = 0; i < num_consequences; i++) {
        // csq_preds->hap_score[i] = calloc(1, sizeof(kstring_t));
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
        // if (csq_preds->hap_score[i]->m > 0) free(csq_preds->hap_score[i]->s);
        // free(csq_preds->hap_score[i]);

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

    // free(csq_preds->hap_score);
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
HdrFormat *parse_csq_hdr_fmt(bcf_hdr_t *hdr) {
    char *hdr_info;
    if (bcf_hdr_id2int(hdr, BCF_DT_ID, csq_field) < 0) return NULL;

    // Get CSQ Description
    bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, BCF_HL_INFO, NULL, csq_field, NULL);
    int key = bcf_hrec_find_key(hrec, hdr_desc);
    const char *description = strstr(hrec->vals[key], vep_desc_prefix) + vep_desc_prefix_len; // Skip "Format: "
    int len = strlen(description)-1;

    // Split CSQ description
    kstring_t s = {0}; kputsn(description, len, &s);
    int *field_idxs, n; field_idxs = ksplit(&s, vep_delim, &n);

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
    bcf_get_info_string(hdr, v, csq_field, &vep_str, &n_vep_values);

    VepCsqs *transcripts = malloc(sizeof(VepCsqs)); kv_init(*transcripts);

    kstring_t s = {0}; kputs(vep_str, &s);
    int *field_idxs, n; field_idxs = str_split(&s, vcf_delim, &n);
    for (int i = 0; i < n; i++) {
        int flen = strlen(s.s + field_idxs[i]);
        kstring_t transcript_csq = { flen, flen+1, s.s + field_idxs[i] }; // WARN: NOT MALLOC'ED, SHOULD NOT BE FREED
        int *vep_field_idxs, vep_n; vep_field_idxs = str_split(&transcript_csq, vep_delim, &vep_n);

        int feature_type_id; khash_str2int_get(vep_fmt->str2id, vep_feature_type, &feature_type_id);
        if (strcmp(transcript_csq.s + vep_field_idxs[feature_type_id], vep_feature_type_transcript) != 0) { free(vep_field_idxs); continue; }

        int feature_id; khash_str2int_get(vep_fmt->str2id, vep_feature, &feature_id);
        int consequences_id; khash_str2int_get(vep_fmt->str2id, vep_consequence, &consequences_id);

        // Parse consequences for transcript, vep allows multiple with '&' separators
        int clen = strlen(transcript_csq.s + vep_field_idxs[consequences_id]);
        kstring_t consequences = { clen, clen+1, transcript_csq.s + vep_field_idxs[consequences_id] };
        int *consequence_field_idxs, consequence_n; consequence_field_idxs = str_split(&consequences, consequence_delim, &consequence_n);
        char **cons = malloc(consequence_n * sizeof(char *));
        for (int j = 0; j < consequence_n; j++) {
            cons[j] = strdup(consequences.s + consequence_field_idxs[j]);
        }

        // Store in hashmap
        char *tr_id = transcript_csq.s + vep_field_idxs[feature_id];
        TscriptCsqs *t = malloc(sizeof(TscriptCsqs));
        t->tid = strdup(tr_id);
        t->num_consequences = consequence_n;
        t->consequences = cons;
        for (int j = 0; j < 3; j++) {
            t->triggers[j] = 0;
        }
        kv_push(TscriptCsqs *, *transcripts, t);

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
        TscriptCsqs *vep_csq_field = kv_A(*vep_csqs, i);

        for (int j = 0; j < vep_csq_field->num_consequences; j++) {
            free(vep_csq_field->consequences[j]);
        }
        free(vep_csq_field->consequences);
        free(vep_csq_field->tid);

        free(vep_csq_field);
    }

    kv_destroy(*vep_csqs);
    free(vep_csqs);
}

Record *init_record(bcf1_t *v, VepCsqs *csqs) {
    Record *r = malloc(sizeof(Record));

    bcf1_t *v_dup = bcf_dup(v);
    bcf_unpack(v_dup, BCF_UN_ALL);
    r->v = v_dup;

    r->c = NULL;
    r->genes = NULL;
    r->tscripts = NULL;
    r->transcript_genes = NULL;
    r->vep = NULL;
    r->isoforms = NULL;
    r->complete_scores = NULL;

    r->skipped = 0;
    r->unphased = 0;

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

