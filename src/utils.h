#include <stdint.h>
#include <htslib/vcf.h>

#include <klib/kstring.h>
#include <klib/kvec.h>


// typedef struct Transcript {
//     int rid;
//     gf_tscript_t *definition;
//     kstring_t seq;
// } Transcript;
//
// typedef kvec_t(Transcript) Transcripts;
//
// typedef struct VepTranscript {
//     char *tr_id;
//     char **consequences;
// } VepTranscript;
//
// typedef kvec_t(VepTranscript) VepTranscripts;

enum Trigger {
    SPLICE,
    FRAMESHIFT,
    STOP
};

typedef kvec_t(bcf1_t *) Variants;

typedef struct Gene {
    Range coords;

} Gene;

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

typedef struct Record {
    bcf1_t *v;
    Variants *hap_incl;
    Variants *hap_excl;
    enum Trigger triggers[3];
    ConsequencePredictions *c;
} Record;

typedef kvec_t(Record *) Records;

int get_allele_offset(int64_t pos, int rlen, int alen, int boundary_start, int boundary_end);

int get_allele_offset_variant(bcf1_t *v, int boundary_start, int boundary_end);

int64_t get_position_in_haplotype(int64_t pos, Variants *vs, int64_t start, int64_t end);

int get_haplotype_sequence(kstring_t ref, Variants vs, int seq_offset, int pos_offset, kstring_t *out_hap);


