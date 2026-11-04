#ifndef ANNOTATION_HEADER
#define ANNOTATION_HEADER

#include "record.h"
#include "region.h"
#include "session.h"

// What the annotation holds for one record. The entries of tscripts and transcript_genes line up with the record's
// VEP consequences: an entry is NULL where the consequence has no coding transcript in the GFF.
// All pointers point at objects owned by the Annotation, nothing here is a copy.
typedef struct RecordAnnotation {
    Transcripts tscripts;
    Genes transcript_genes;     // The gene of each transcript, NULL wherever tscripts is NULL
    Genes genes;                // The distinct genes of the record, in order of first appearance
} RecordAnnotation;

// The genes and transcripts of a block of records, built once and shared by all of its records. Everything in it is
// read-only after the build. Transcripts and exons are half-open (see Transcript_build), gene ranges are half-open too.
typedef struct Annotation {
    kh_TranscriptHash_t *transcripts;   // Owns the transcripts, keyed by GFF transcript id
    kh_GeneHash_t *genes;               // Owns the genes, keyed by GFF gene id
    RecordAnnotation *records;          // One per record of the block, in the order of the block
    size_t n_records;
} Annotation;

// Builds the annotation for the records of a block, which must come from a coding gene region (so that they have VEP
// consequences). A record without matching coding transcripts gets empty vectors. Returns NULL, after logging, if a gene
// sequence cannot be fetched.
Annotation *Annotation_build(const Records *records, const Session *s);

void Annotation_destroy(Annotation *annotation);

// Finds the GFF transcripts of a record: one entry per VEP consequence of the record, NULL where the consequence does not
// match a protein-coding transcript. The caller frees the result with kv_destroy and free.
int find_gf_tscripts_by_iregion(IRegion iregion, bcf_hdr_t *hdr, gff_t *gff, gf_tscripts **out_gf_tscripts);
void match_gf_tscripts_to_vep_tscripts(gf_tscripts *tscripts, Record *r, gff_t *gff, gf_tscripts **out_gf_scripts);
int find_gf_tscripts_for_record(Record *r, bcf_hdr_t *hdr, gff_t *gff, gf_tscripts **out_gf_tscripts);

#endif
