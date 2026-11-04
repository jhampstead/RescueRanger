#include "reference_scores.h"

#include "logging/log.h"

#ifdef CPLICEAI_HAS_REFERENCE_SCORES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/khash_str2int.h>
#include <cpliceai/gene_reference.h>
#include <cpliceai/reference.h>
#include <cpliceai/utils.h>

struct ReferenceScores {
    Reference reference;
    faidx_t *fai;

    // The gene that is loaded, which is the one the last call asked for
    GeneReference gene;
    int has_gene;
    char gene_chrom[FIELD_MAX_LEN];

    // The genes that could not be used, so that the warning is given once
    void *refused;
    long windows;
};

int ReferenceScores_supported(void) {
    return 1;
}

// reference_read maps the file, and reference_free would free the mapping as if it were memory from malloc. The mapping begins
// at the file header, which is where the block of the reference ends up after: off_block is the size of the header
static void unmap_reference(Reference *reference) {
    munmap(reference->block - sizeof(FileHeader), reference->block_size + sizeof(FileHeader));
}

ReferenceScores *ReferenceScores_open(const char *reference_path, faidx_t *fai) {
    Reference reference = {0};
    if (reference_read(reference_path, &reference) != EXIT_SUCCESS) {
        log_error("Could not read the reference scores: %s", reference_path);
        return NULL;
    }

    // Only the FASTA is checked here. The file also holds a fingerprint of the gene regions it was made from, which is not
    // checked: it cannot be made from the annotation, and every gene is checked against the annotation when it is used. The value
    // the file has is given as the one to compare with, so that only the FASTA can fail
    if (reference_check_inputs(&reference, fasta_digest(fai), reference.regions_digest, reference_path) != EXIT_SUCCESS) {
        unmap_reference(&reference);
        return NULL;
    }

    ReferenceScores *scores = calloc(1, sizeof(ReferenceScores));
    scores->reference = reference;
    scores->fai = fai;
    scores->refused = khash_str2int_init();
    gene_reference_init(&scores->gene);

    log_info("Reference scores: %zu genes on %zu contigs in %s", reference.n_genes, reference.n_contigs, reference_path);

    return scores;
}

void ReferenceScores_counts(const ReferenceScores *scores, long *windows, long *refused_genes) {
    *windows = scores == NULL ? 0 : scores->windows;
    *refused_genes = scores == NULL ? 0 : (long) khash_str2int_size(scores->refused);
}

void ReferenceScores_close(ReferenceScores *scores) {
    if (scores == NULL) return;

    log_info("Reference scores: %ld windows read from them, %d genes were not in them like the annotation's and were predicted from the sequence",
             scores->windows, (int) khash_str2int_size(scores->refused));

    gene_reference_destroy(&scores->gene);
    khash_str2int_destroy_free(scores->refused);
    unmap_reference(&scores->reference);
    free(scores);
}

// The index of the first gene of that name on the contig, or -1. With start and end it is the first that also has those coordinates
// and that strand
static long find_gene(const Reference *reference, const char *chrom, const char *name, int match_coordinates, int64_t start, int64_t end, char strand) {
    for (size_t c = 0; c < reference->n_contigs; c++) {
        const Contig *contig = &reference->contigs[c];
        if (strncmp(reference->contig_names + contig->name_start, chrom, FIELD_MAX_LEN) != 0) continue;

        for (size_t g = contig->region_start; g < contig->region_start + contig->n_regions; g++) {
            const Region *region = &reference->genes[g];
            if (strncmp(reference->region_names + region->name_start, name, FIELD_MAX_LEN) != 0) continue;
            if (match_coordinates && ((int64_t) reference->gene_starts[g] != start || (int64_t) reference->gene_ends[g] != end || region->strand != strand)) continue;

            return (long) g;
        }

        return -1;
    }

    return -1;
}

static float *refuse(ReferenceScores *scores, const char *chrom, const char *name, int64_t start, int64_t end, char strand) {
    char key[3 * FIELD_MAX_LEN + 64];
    snprintf(key, sizeof(key), "%s:%s:%ld-%ld:%c", chrom, name, (long) start, (long) end, strand);

    if (!khash_str2int_has_key(scores->refused, key)) {
        const long same_name = find_gene(&scores->reference, chrom, name, 0, 0, 0, 0);
        if (same_name < 0) {
            log_warn("The reference scores have no gene %s on %s, so its reference predictions are made from the sequence (the annotation has it at %ld-%ld, %c)",
                     name, chrom, (long) start, (long) end, strand);
        } else {
            log_warn("The reference scores have gene %s on %s at %lu-%lu (%c), the annotation at %ld-%ld (%c): the reference predictions of this gene are made from the sequence",
                     name, chrom, (unsigned long) scores->reference.gene_starts[same_name], (unsigned long) scores->reference.gene_ends[same_name],
                     scores->reference.genes[same_name].strand, (long) start, (long) end, strand);
        }
        khash_str2int_set(scores->refused, strdup(key), 1);
    }

    scores->has_gene = 0;
    return NULL;
}

float *ReferenceScores_window(ReferenceScores *scores, const char *chrom, const char *gene_name, int64_t gene_start, int64_t gene_end,
                              char strand, int64_t pos, int radius) {
    if (scores == NULL) return NULL;

    const int loaded = scores->has_gene
        && strncmp(scores->gene_chrom, chrom, FIELD_MAX_LEN) == 0
        && strncmp(scores->gene.name, gene_name, FIELD_MAX_LEN) == 0
        && (int64_t) scores->gene.start == gene_start
        && (int64_t) scores->gene.end == gene_end
        && scores->gene.strand == strand;

    if (!loaded) {
        scores->has_gene = 0;

        if (find_gene(&scores->reference, chrom, gene_name, 1, gene_start, gene_end, strand) < 0) return refuse(scores, chrom, gene_name, gene_start, gene_end, strand);

        // This finds the first gene of that name on the contig, which is not necessarily the one that matched above
        if (gene_reference_update(chrom, gene_name, scores->fai, &scores->reference, &scores->gene) != EXIT_SUCCESS) {
            return refuse(scores, chrom, gene_name, gene_start, gene_end, strand);
        }

        if ((int64_t) scores->gene.start != gene_start || (int64_t) scores->gene.end != gene_end || scores->gene.strand != strand) {
            return refuse(scores, chrom, gene_name, gene_start, gene_end, strand);
        }

        strncpy(scores->gene_chrom, chrom, FIELD_MAX_LEN - 1);
        scores->gene_chrom[FIELD_MAX_LEN - 1] = '\0';
        scores->has_gene = 1;
    }

    float *window = NULL;
    if (gene_reference_get_score_window(pos, radius, &scores->gene, &window) != EXIT_SUCCESS) return NULL;

    scores->windows++;
    return window;
}

#else

int ReferenceScores_supported(void) {
    return 0;
}

ReferenceScores *ReferenceScores_open(const char *reference_path, faidx_t *fai) {
    (void) fai;
    log_error("Cannot use the reference scores %s: this build uses a cpliceai without them. Build against a newer cpliceai "
              "(-DCPLICEAI_ROOT=...), or leave the option out and the reference predictions are made from the sequence", reference_path);

    return NULL;
}

void ReferenceScores_close(ReferenceScores *scores) {
    (void) scores;
}

float *ReferenceScores_window(ReferenceScores *scores, const char *chrom, const char *gene_name, int64_t gene_start, int64_t gene_end,
                              char strand, int64_t pos, int radius) {
    (void) scores; (void) chrom; (void) gene_name; (void) gene_start; (void) gene_end; (void) strand; (void) pos; (void) radius;
    return NULL;
}

void ReferenceScores_counts(const ReferenceScores *scores, long *windows, long *refused_genes) {
    (void) scores;
    *windows = 0;
    *refused_genes = 0;
}

#endif
