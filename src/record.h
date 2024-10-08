#ifndef RECORD_HEADER
#define RECORD_HEADER

#include <htslib/khash.h>
#include <htslib/vcf.h>
#include <htslib/khash_str2int.h>
#include <klib/khash.h>
#include <stdint.h>

#include "klib/kvec.h"
#include "bcftools/gff.h"

#include "cpliceai/range.h"
#include "cpliceai/utils.h"


typedef struct Exon {
    Range *cds;
    Range *bounds;
} Exon;

typedef struct Transcript
{
    uint32_t id;        // transcript id
    uint32_t start,end;   // transcript's beg and end coordinate (ref strand, 0-based, inclusive)
    uint32_t strand:2,  // STRAND_REV,FWD,UNK
             ncds:30,   // number of exon coding regions
             mcds;
    Range **cds;     // ordered list of exon coding regions
    uint32_t nexons, mexons;
    Exon **exons;
} Transcript;

typedef struct Gene {
    uint32_t gff_id;
    char *name;
    Range coords;
    kstring_t seq;
    uint32_t strand;
} Gene;

typedef kvec_t(Gene *) Genes;

KHASH_MAP_INIT_INT64(RangeHash, Range *);
KHASH_MAP_INIT_INT(GeneHash, Gene *);

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

typedef struct TscriptCsqs {
    char *tid;
    int num_consequences;
    char **consequences;
} TscriptCsqs;

typedef kvec_t(TscriptCsqs) VepCsqs;

typedef kvec_t(gf_tscript_t *) gf_tscripts;

typedef kvec_t(Transcript *) Transcripts;

typedef struct Record {
    bcf1_t *v;
    VepCsqs *vep;
    int triggers[3];
    ConsequencePredictions *c;

    Score **hap_scores;
    Score **complete_scores;
    Transcripts *tscripts;
    Genes *transcript_genes;
    Genes *genes;
} Record;

typedef kvec_t(Record *) Records;

typedef struct HdrFormat {
    char *data_; // Only used for storing data, access this data should be done through id2str
    char **id2str;
    kh_str2int_t *str2id;
} HdrFormat;

void init_consequence_predictions(ConsequencePredictions **dst, int num_consequences);

void destroy_consequence_predictions(ConsequencePredictions *csq_preds);

HdrFormat *parse_csq_format(bcf_hdr_t *hdr);

void destroy_hdr_fmt(HdrFormat *hdr_fmt);

VepCsqs *parse_vep_transcript_consequences(bcf1_t *v, const HdrFormat *vep_fmt, bcf_hdr_t *hdr);

void destroy_vep_csqs(VepCsqs *vep_csqs);

Record *init_record(bcf1_t *v, VepCsqs *csqs);

void destroy_record(Record *r);

#endif
