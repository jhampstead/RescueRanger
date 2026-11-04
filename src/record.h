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

// Number of haplotype phases; the per-phase arrays on Record are indexed as [index * NUM_PHASES + phase]
#ifndef NUM_PHASES
#define NUM_PHASES 2
#endif

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
    INFRAME = 1 << 22,
    // The NMD escape region the new stop of an unrescued or induced PTV is in, on the most probable isoform with such a PTV. At
    // most one is set per transcript and phase; none when the stop is outside every escape region or no stop was found
    NMD_START_PROXIMAL = 1 << 23,
    NMD_LONG_EXON = 1 << 24,
    NMD_DISTAL = 1 << 25
};

typedef struct Isoform {
    float prior_probability; // prior probability of the isoform
    float posterior_probability; // posterior probability of the isoform
    uint8_t rescue; // 0 = not rescued, 1 = rescued
    enum RescueType rescue_type; // representing the type of rescue event predicted (e.g. partial_exon_skipping)
    uint32_t affected_bases; // number of bases offset/excluded/included within the isoform sequence
    uint8_t ptv; // 0 = no non-canonical PTV on isoform, 1 = non-canonical PTV on isoform
    enum PTVType ptv_type; // representing the type of non-canonical PTV predicted (stop_gained or frameshift_variant)
    uint8_t nmd_ptv; // 1 = the isoform keeps or induces a PTV, which the NMD flag is taken from (also an unrescued VEP stop or frameshift)
    uint32_t nmd_flag; // NMD_START_PROXIMAL, NMD_LONG_EXON or NMD_DISTAL for the escape region of the new stop of that PTV, 0 for none
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
    kvec_t(Isoform) isoforms; // the isoforms, one per combination of non-coding regions
    uint32_t n_rescue_types; // count of unique rescue types
    uint32_t m_rescue_types; // capacity for unique rescue types
    uint8_t rescue; // overall rescue prediction (0 = not rescued, 1 = rescued)
    float rescue_probability; // total rescue probability
    uint8_t ptv; // prediction of non-canonical PTV
    float ptv_probability; // cumulative probability of non-canonical PTV
    enum PTVType ptv_type; // type of new PTV, one of stop_gained or frameshift_variant TODO: add support for variants with multiple PTV types
    RescueProbabilities *rescue_type_probabilities; // pointer to RescueProbabilities struct
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
    SpliceSites *splice_sites;
    bool constitutive;
    bool in_frame;
} Exon;

typedef kvec_t(Exon *) Exons;

// The gap between two consecutive exons, in ascending reference order. All values are half-open.
typedef struct Intron {
    int64_t start, end;     // [start, end), i.e. the first base after the upstream exon up to the next exon's start
    int64_t len;            // end - start
    int64_t coding_offset;  // total length of all introns up to and including this one
} Intron;

typedef kvec_t(Intron) Introns;

typedef struct Transcript {
    uint32_t id;        // transcript id
    uint32_t start, end;   // transcript's start and end coordinate (ref strand, 0-based, end exclusive)
    uint32_t strand:2,  // STRAND_REV,FWD,UNK
             ncds:30,   // number of exon coding regions
             mcds;
    Range *cds;     // ordered list of exon coding regions
    Exons exons;
    Introns introns;   // filled by Transcript_build
} Transcript;

typedef struct Gene {
    uint32_t gff_id;
    char *name;
    Range coords;       // half-open; seq holds the bases of [coords.start, coords.end)
    kstring_t seq;
    uint32_t strand;
} Gene;

typedef kvec_t(Range) Ranges;
typedef kvec_t(Gene *) Genes;

// Derives the introns from exon ranges that are half-open and sorted ascending.
// Returns exactly n_exons-1 introns (none for 0 or 1 exons), so intron i lies between exon i and exon i+1.
// The caller frees the result with kv_destroy.
Introns Introns_from_exons(const Range *exons, size_t n_exons);

// Splice sites are single bases, never exclusive ends, so these hold the +1/-1 arithmetic in one place.
// A donor is the last exon base on the forward strand and the first on the reverse strand, and an acceptor is the
// opposite. An unknown strand is treated as forward.
int64_t Exon_donor_base(const Exon *exon, uint32_t strand);
int64_t Exon_acceptor_base(const Exon *exon, uint32_t strand);

// The ascending half-open intron between a donor base and an acceptor base of the same transcript.
Range Intron_between_sites(uint32_t strand, int64_t donor_base, int64_t acceptor_base);

// The inverse of Intron_between_sites, for an intron range.
int64_t Intron_donor_base(Range intron, uint32_t strand);
int64_t Intron_acceptor_base(Range intron, uint32_t strand);

// The ascending half-open range a UTR covers. boundary_cds is the CDS of the coding exon the UTR touches: for the
// 5' UTR the first coding exon in transcript order, for the 3' UTR the last. gene is the half-open gene range.
Range Utr_range(enum UTRType type, uint32_t strand, Range gene, Range boundary_cds);

// The range of reference positions in which a predicted splice site counts as belonging to this exon, taken from the
// ref_bounds of the exon and its neighbours (predicted sites carry reference coordinates, so the variants of a haplotype
// do not move the window). Sites on the upper border of an exon (forward donor, reverse acceptor) look towards the exon
// above, sites on the lower border (forward acceptor, reverse donor) towards the exon below. Returns false if there is no
// such neighbour, or the strand is unknown.
bool Exon_splice_window(const Transcript *transcript, size_t exon_index, bool donor, Range *window);

// One transcript as read from the GFF. Every range is 0-based inclusive, exactly as the GFF reader returns it.
typedef struct TranscriptData {
    uint32_t id;
    uint32_t start, end;
    uint32_t strand;
    const Range *cds;           // coding regions, ascending
    const int *constitutive;    // one flag per coding region
    size_t ncds;
    const Range *exons;         // ascending
    size_t nexons;
} TranscriptData;

// Builds a Transcript from GFF data. This is the one place the inclusive GFF coordinates are converted:
// the returned transcript, its CDS ranges, exons' ref_bounds and introns are all half-open.
// Exons get a CDS the same way set_transcript_exons does: walking the exons in order, an exon takes the next
// unused coding region if the two intersect, otherwise it has none (cds NULL) and the region stays unused.
// Exon.cds points into Transcript.cds. splice_sites is NULL.
Transcript *Transcript_build(const TranscriptData *data);

// Frees a transcript made by Transcript_build, including its exons and introns.
void Transcript_destroy(Transcript *transcript);

// The NMD escape regions of Regional nonsense constraint offers biological and clinical insights into genetic disease (Nature
// Communications, https://doi.org/10.1038/s41467-026-69983-z), on the reference transcript:
//   start-proximal: the first 150 coding bases (the A of the ATG is base 0)
//   long exon: the bases of an exon that are more than 400 bases before its last base in transcript order (so never in the
//              last exon, which has no donor)
//   distal: the last exon, the last 50 bases of the penultimate exon, and the intron between them (a retained intron there is
//           distal too). A transcript of one exon is distal as a whole.
// A stop codon is placed by its first base. Every range is half-open and ascending, whatever the strand, so one test serves
// both strands. An unknown strand is read as forward.
#define NMD_START_PROXIMAL_LENGTH 150
#define NMD_LONG_EXON_LENGTH 400
#define NMD_PENULTIMATE_EXON_LENGTH 50

typedef struct NmdRegions {
    Ranges start_proximal;   // One range per coding region that holds some of the first 150 coding bases, ascending
    Ranges long_exon;        // One range per exon that has bases more than 400 before its last base, ascending
    Ranges distal;           // One range, or none for a transcript without exons
} NmdRegions;

NmdRegions Transcript_nmd_regions(const Transcript *transcript);

// Frees the three ranges of the regions. Nothing is freed for the struct itself.
void NmdRegions_destroy(NmdRegions *regions);

// Frees the four site arrays of a SpliceSites, the Records vector every site points at, and the struct itself.
void SpliceSites_destroy(SpliceSites *sites);

// A copy of a shared transcript for one record in one phase, which is the only thing that is mutated while predicting.
// The copy has its own Transcript and its own Exon structs, but still points at the shared, read-only ref_bounds, Exon.cds,
// cds and introns, so those must outlive the copy and must not be freed through it. splice_sites starts as NULL on every
// exon. Free it with Transcript_phase_destroy, never Transcript_destroy.
Transcript *Transcript_phase_copy(const Transcript *shared);

// Frees a copy made by Transcript_phase_copy: its exons together with their splice_sites. The shared transcript the copy
// was made from is not touched.
void Transcript_phase_destroy(Transcript *copy);

KHASH_MAP_INIT_INT(GeneHash, Gene *);
KHASH_MAP_INIT_INT(TranscriptHash, Transcript *);

typedef struct TscriptCsqs {
    char *tid;
    int num_consequences;
    char **consequences;
    int triggers[3]; // From the consequences and the allele frequency: the same in every phase
    uint8_t reading_frame[NUM_PHASES]; // Per phase: the record has a reading frame conflict with another record of the haplotype
} TscriptCsqs;

typedef kvec_t(TscriptCsqs *) VepCsqs;

typedef kvec_t(gf_tscript_t *) gf_tscripts;

typedef kvec_t(Transcript *) Transcripts;

// The alleles of a record without the prefix REF and ALT share, which for a normalised indel is its VCF anchor base. An
// insertion has ref_len 0 and a deletion alt_len 0. Only the first ALT is read; a record without ALT reads as an empty one.
typedef struct Alleles {
    int64_t pos; // 0-based position of the first base after the shared prefix
    int prefix_len; // Number of shared leading bases left out
    int ref_len; // REF bases after the prefix
    int alt_len; // ALT bases after the prefix
    const char *alt; // ALT bases after the prefix; points into the bcf1_t the view was made from
} Alleles;

typedef struct Record {
    bcf1_t *v;
    Alleles alleles; // Anchor-free view of v's alleles
    VepCsqs *vep;
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

// Which program wrote the consequence annotation of the input
typedef enum AnnotationKind {
    ANNOTATION_VEP,      // Ensembl VEP: INFO/CSQ, with a Feature_type and a Feature column
    ANNOTATION_BCFTOOLS, // bcftools csq: INFO/BCSQ, with a transcript column; the names of the consequences differ from VEP's
} AnnotationKind;

typedef struct HdrFormat {
    char *data_; // Only used for storing data, access this data should be done through id2str
    char **id2str;
    kh_str2int_t *str2id;
    AnnotationKind kind;
    char *tag; // The INFO field the annotation is in
    int consequence_id, feature_type_id, feature_id; // Columns of an entry; feature_type_id is -1 for bcftools, feature_id is its transcript column
} HdrFormat;

void init_splice_site(SpliceSite **site, size_t *n, size_t *m);
void resize_splice_site(SpliceSite **site, size_t *m);

// The tag the annotation is in when none is given: CSQ for VEP, BCSQ for bcftools csq.
const char *Annotation_default_tag(AnnotationKind kind);

// Reads the Format of the INFO field `tag` (NULL or empty for the default of the kind) and checks that it has the columns the kind
// needs: Consequence, Feature_type and Feature for VEP; Consequence and transcript for bcftools. Returns NULL after logging the
// reason if the field is not in the header or has other columns, which is what a file of the other kind looks like.
HdrFormat *parse_csq_hdr_fmt(bcf_hdr_t *hdr, AnnotationKind kind, const char *tag);

void destroy_hdr_fmt(HdrFormat *hdr_fmt);

// The consequences of a record per transcript, from the annotation of the kind the format is for, as VEP names them: the consequences
// of bcftools csq are mapped to those, a leading '*' is dropped, entries that point to another record ('@position') and entries without
// a transcript are left out, and entries of one transcript are joined.
VepCsqs *parse_vep_transcript_consequences(bcf1_t *v, const HdrFormat *fmt, bcf_hdr_t *hdr);

void destroy_vep_csqs(VepCsqs *vep_csqs);

// The anchor-free view of an unpacked record's REF and first ALT. Valid as long as v is.
Alleles Alleles_from_variant(const bcf1_t *v);

// Whether the variant of v (its REF and first ALT, as Alleles_from_variant reads them) belongs to the region, such as a gene:
// variant_in_region on the view. A base the variant changes has to lie in the region; an insertion changes none, so both bases
// around it have to. The load-time gene lookup uses this, so it agrees with the later stages. v must be unpacked up to BCF_UN_STR.
bool bcf1_variant_in_region(const bcf1_t *v, const Range region);

Record *Record_init(bcf1_t *v, VepCsqs *csqs);

// Sets a flag (an enum PredictionFlags value) on a transcript of the record for a phase. The record needs its flag slots, one
// per transcript and phase.
void set_record_flag(Record *r, int transcript_index, int phase_index, uint32_t flag_value);

// Frees what a built (valid) Isoforms owns: the isoform vector and the rescue probabilities. Not the struct.
void Isoforms_destroy(Isoforms *isoforms);

// Frees the record with everything it owns, including the flags and every built isoform slot. The Gene and Transcript objects
// its vectors point at belong to the annotation of the block (or the phase copies) and are not freed.
void Record_destroy(Record *r);

#endif
