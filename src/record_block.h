#ifndef RECORD_BLOCK_HEADER
#define RECORD_BLOCK_HEADER

#include <stdbool.h>
#include <htslib/vcf.h>

#include "annotation.h"
#include "record.h"
#include "session.h"

typedef struct RecordBuffer {
    bcf1_t *v;
    bool unused; // v holds a record that has not been taken into a block yet
    bool failed; // Reading failed. It is kept, so that the failure is reported after the records that were read before it
} RecordBuffer;

enum ReadStatus {
    READ_OK,    // A record was read
    READ_END,   // The end of the input
    READ_FAILED // The input could not be read
};

// What the result of bcf_read and the error code of the record it read into mean: a result of 0 is a record, -1 without an error
// code is the end of the input, and every other negative result, or an error code on the record, is a failure. A failure used to
// be taken for the end of the input, so a damaged file gave a shorter output and an exit status of success.
enum ReadStatus classify_read(int bcf_read_result, int errcode);

// Reads the next record into the buffer, unless it holds one that has not been taken yet (then it is READ_OK and nothing is read).
// A failure is logged and kept in the buffer: after it nothing is read again and the result is READ_FAILED.
enum ReadStatus RecordBuffer_fill(htsFile *fp, bcf_hdr_t *hdr, RecordBuffer *buffer);

// A batch of consecutive records whose coding gene regions overlap, or a single record outside every coding gene.
typedef struct RecordBlock {
    int is_coding;
    Records records;
    Annotation *annotation; // NULL unless is_coding. Owned by the block, and shared by all of its records and phases
} RecordBlock;

// Reads the next block. For a coding block this also builds the annotation, and gives every record that has coding
// transcripts its gene vectors (the Gene objects are shared, owned by the annotation) and its score, isoform and flag
// arrays, once. A record's tscripts vector has one NULL entry per VEP consequence until RecordBlock_begin_phase.
// Returns 0 if a block was read, 1 at the end of the input, and 2 if the annotation of the block could not be built or the input
// could not be read. The records that were read before a read failure are returned first, as a block, and the failure with the
// next call.
int RecordBlock_next(const Session *s, RecordBuffer *buffer, RecordBlock *block_out);

void RecordBlock_destroy(RecordBlock block);

// Gives every record with coding transcripts its own copies of the transcripts for one phase, in record->tscripts.
// The copies are what the predictions mutate (the splice_sites), the annotation itself stays read-only.
void RecordBlock_begin_phase(RecordBlock *block);

// Destroys the copies made by RecordBlock_begin_phase and sets the entries back to NULL. Call it once the phase is done,
// because the non-coding regions of a phase point into the splice sites of its copies.
void RecordBlock_end_phase(RecordBlock *block);

#endif
