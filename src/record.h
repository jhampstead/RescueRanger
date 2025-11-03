#ifndef RECORD_HEADER
#define RECORD_HEADER

#include <stdbool.h>
#include <htslib/khash.h>
#include <htslib/vcf.h>
#include <htslib/khash_str2int.h>
#include <klib/khash.h>
#include <stdint.h>

#include "klib/kvec.h"
#include "bcftools/gff.h"

#include "cpliceai/range.h"
#include "cpliceai/utils.h"

enum UTRType {
    UTR5 = 1,
    UTR3 = 2
};

enum RescueType {
    NO_RESCUE = 0,
    COMPLETE_RESCUE,
    POSSIBLE_REINITIATION_SITE,
    FRAME_RESTORING_INDEL,
    PARTIAL_EXON_SKIPPING,
    PARTIAL_INTRON_RETENTION,
    EXON_SKIPPING,
    INTRON_RETENTION,
};

enum PTVType {
    NO_PTV = 0,
    STOP_GAINED,
    FRAMESHIFT,
    LONG_INFRAME_INDEL
};

enum PredictionFlags {
    UNPHASED_VARIANT_IN_PREDICTION = 1 << 0,
    SKIPPED_VARIANT_IN_PREDICTION = 1 << 1,
    VARIANT_LARGER_THAN_SPLICE_WINDOW = 1 << 2,
    RESCUED_BY_SPLICE_FRAMESHIFT = 1 << 3,
    POSSIBLE_REINITIATION = 1 << 4,
    MANUFACTURED_SPLICE_RESCUE = 1 << 5,
    HAPLOTYPE_SPECIFIC_DONOR = 1 << 6,
    HAPLOTYPE_SPECIFIC_ACCEPTOR = 1 << 7,
    CRYPTIC_DONOR = 1 << 8,
    CRYPTIC_ACCEPTOR = 1 << 9,
    WEAK_CANONICAL_DONOR = 1 << 10,
    WEAK_CANONICAL_ACCEPTOR = 1 << 11,
    ABLATED_CANONICAL_DONOR = 1 << 12,
    ABLATED_CANONICAL_ACCEPTOR = 1 << 13,
    ABLATED_CRYPTIC_DONOR = 1 << 14,
    ABLATED_CRYPTIC_ACCEPTOR = 1 << 15,
    UPSTREAM_PTV = 1 << 16,
    UTR5_VARIANT = 1 << 17,
    UTR3_VARIANT = 1 << 18,
    LAST_EXON = 1 << 19,
    SINGLE_EXON = 1 << 20,
    NON_CONSTITUTIVE = 1 << 21,
    INFRAME = 1 << 22
};

typedef struct Isoform {
    // NonCodingRegion *ncrs; // non-coding regions defining the spliced boundaries of this isoform
    float prior_probability; // prior probability of the isoform
    float posterior_probability; // posterior probability of the isoform
    uint8_t rescue; // 0 = not rescued, 1 = rescued
    enum RescueType rescue_type; // representing the type of rescue event predicted (e.g. partial_exon_skipping)
    uint32_t affected_bases; // number of bases offset/excluded/included within the isoform sequence
    uint8_t ptv; // 0 = no non-canonical PTV on isoform, 1 = non-canonical PTV on isoform
    enum PTVType ptv_type; // representing the type of non-canonical PTV predicted (stop_gained or frameshift_variant)
} Isoform;

typedef struct RescueProbabilities {
    uint32_t *n_isoforms; // counts of isoforms supporting each rescue type
    uint8_t *rescues;
    enum RescueType *rescue_types; // rescue types
    float *rescue_probabilities; // cumulative probability for each rescue type
    uint32_t *affected_bases; // counts of affected bases
} RescueProbabilities;

typedef struct Isoforms {
    uint8_t valid;
    uint32_t n_isoforms; // number of isoforms
    uint32_t m_isoforms; // capacity for isoforms
    uint32_t n_rescue_types; // count of unique rescue types
    uint32_t m_rescue_types; // capacity for unique rescue types
    uint8_t rescue; // overall rescue prediction (0 = not rescued, 1 = rescued)
    float rescue_probability; // total rescue probability
    uint8_t ptv; // prediction of non-canonical PTV
    float ptv_probability; // cumulative probability of non-canonical PTV
    enum PTVType ptv_type; // type of new PTV, one of stop_gained or frameshift_variant TODO: add support for variants with multiple PTV types
    Isoform **isoforms; // array of pointers to Isoform structs
    RescueProbabilities *rescue_type_probabilities; // pointer to RescueProbabilities struct
    char **flags; // array of prediction flags
} Isoforms;

typedef struct SpliceSite {
    void *records_ptr; // Circular dependency?
    uint32_t position;
    uint32_t canonical_position;
    uint32_t predicted; // whether haplotype prediction has been run for the splice site
    uint32_t gff_canonical; // whether the splice site is canonical in the GFF
    uint32_t lost;
    float ref_prior_probability;
    float prior_probability; // raw SpliceAI probability
    float posterior_probability; // scaled probability across window
} SpliceSite;

typedef struct SpliceSites {
    SpliceSite *canonical_donors;
    SpliceSite *canonical_acceptors;
    SpliceSite *cryptic_donors;
    SpliceSite *cryptic_acceptors;
    size_t n_canonical_donors, m_canonical_donors;
    size_t n_canonical_acceptors, m_canonical_acceptors;
    size_t n_cryptic_donors, m_cryptic_donors;
    size_t n_cryptic_acceptors, m_cryptic_acceptors;
} SpliceSites;

typedef struct Exon {
    Range *cds;
    Range *ref_bounds;
    Range *hap_bounds;
    SpliceSites *splice_sites;
    bool constitutive;
    bool in_frame;
} Exon;

typedef struct Transcript {
    uint32_t id;        // transcript id
    uint32_t start, end;   // transcript's start and end coordinate (ref strand, 0-based, inclusive)
    uint32_t strand:2,  // STRAND_REV,FWD,UNK
             ncds:30,   // number of exon coding regions
             mcds;
    Range *cds;     // ordered list of exon coding regions
    size_t nexons, mexons;
    Exon **exons;
    int triggers[3];
} Transcript;

typedef struct Gene {
    uint32_t gff_id;
    char *name;
    Range coords;
    kstring_t seq;
    uint32_t strand;
} Gene;

typedef kvec_t(Range) Ranges;
typedef kvec_t(Gene *) Genes;

KHASH_MAP_INIT_INT64(RangeHash, Range *);
KHASH_MAP_INIT_INT(GeneHash, Gene *);
KHASH_MAP_INIT_INT(TranscriptHash, Transcript *);

typedef struct ConsequencePredictions {
    // kstring_t **hap_score;
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
    int triggers[3];
} TscriptCsqs;

typedef kvec_t(TscriptCsqs *) VepCsqs;

typedef kvec_t(gf_tscript_t *) gf_tscripts;

typedef kvec_t(Transcript *) Transcripts;

typedef struct SpliceChange {
    float AG;
    float AL;
    float DG;
    float DL;
    float sum;
} SpliceChange;

typedef struct Record {
    bcf1_t *v;
    VepCsqs *vep;
    ConsequencePredictions *c;
    Isoforms *isoforms;

    Score **native_scores;
    Score **diff_scores;
    Score **haplotype_scores;
    Transcripts *tscripts;
    Genes *transcript_genes;
    Genes *genes;
    uint32_t *flags;
    uint8_t skipped, unphased;
} Record;

typedef kvec_t(Record *) Records;

typedef struct NonCodingRegion {
    int is_used;
    Range coords;
    Range canonical_coords;
    SpliceSite *donor, *acceptor;
    enum UTRType is_utr;
    int uses_competitor;
    Records records;
} NonCodingRegion;

typedef kvec_t(NonCodingRegion) NonCodingRegions;

typedef struct HdrFormat {
    char *data_; // Only used for storing data, access this data should be done through id2str
    char **id2str;
    kh_str2int_t *str2id;
} HdrFormat;

void init_splice_site(SpliceSite **site, size_t *n, size_t *m);
void resize_splice_site(SpliceSite **site, size_t *m);
void destroy_splice_site(SpliceSite *site, size_t n);

void init_consequence_predictions(ConsequencePredictions **dst, int num_consequences);
void destroy_consequence_predictions(ConsequencePredictions *csq_preds);

HdrFormat *parse_csq_hdr_fmt(bcf_hdr_t *hdr);

void destroy_hdr_fmt(HdrFormat *hdr_fmt);

VepCsqs *parse_vep_transcript_consequences(bcf1_t *v, const HdrFormat *vep_fmt, bcf_hdr_t *hdr);

void destroy_vep_csqs(VepCsqs *vep_csqs);

Record *init_record(bcf1_t *v, VepCsqs *csqs);

void destroy_record(Record *r);

#endif
