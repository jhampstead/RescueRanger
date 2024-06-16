#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <zlib.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>

#include "../bcftools/gff.h"
#include "../bcftools/regidx.h"

#include "logging/log.h"


#define BUFFER_SIZE 2024
#define ID 0
#define TYPE 2
#define BEGIN 3
#define END 4
#define DIRECTION 6
#define PHASE 7
#define ATTRIBUTE 8

#define N_CSQ_FIELDS 29

enum SeqDir {
    FORWARD,
    BACKWARD
};

typedef struct Row {
    char *name;
    char *type;
    int64_t start, end;
    enum SeqDir dir;
    int16_t phase;
    char *id;
    char *transcript_id;
} Row;

typedef struct Table {
    int nrows;
    int nallocs;
    Row *rows;
} Table;

typedef struct CDS {
    int64_t start, end;
    int16_t phase;
} CDS;

typedef struct Transcript {
    char *id;
    char *region;
    int64_t start_codon, end_codon;
    int direction;
    int ncds;
    CDS *cds;
} Transcript;


Row parse_row(char *s) {
    int col = 0, start = 0, len;
    Row row = (Row) {NULL, NULL, 0, 0, 0, 0, NULL, NULL};
    char *as;

    while (*s && *s!='\n') {
        len = 0;
        // Determine length of fragment to delimiter
        while (s[len] && s[len]!='\t' && s[len]!='\n') len++;

        s[len] = '\0';

        switch (col) {
            case ID:
                row.name = strdup(s);
                break;
            case TYPE:
                row.type = strdup(s);
                break;
            case BEGIN:
                row.start = atol(s);
                break; 
            case END:
                row.end = atol(s);
                break;
            case DIRECTION:
                if (*s == '+') row.dir = FORWARD;
                if (*s == '-') row.dir = BACKWARD;
                break;
            case PHASE:
                if (*s == '.') {
                    row.phase = -1;
                    break;
                }
                row.phase = atol(s);
                break;
            case ATTRIBUTE:
                as = s;

                int alen;
                char *delim;
                while (*as && *as!='\n' && *as!='\0') {
                    // Separate single instance of attribute from rest of string
                    alen = 0;
                    while (as[alen] && as[alen]!=';' && as[alen]!='\n') alen++;
                    as[alen] = '\0';

                    // Split key=value pair
                    delim = strchr(as, '=');
                    *delim = '\0';

                    // If transcript_id, store value
                    if (strcmp(as, "transcript_id") == 0) {
                        int id_len = alen - strlen(as);
                        row.id = (char *) malloc(id_len + 1);
                        memcpy(row.id, delim+1, id_len + 1);
                    }

                    as = as + alen + 1;
                }
                break;
        }

        s = s + len + 1;
        col++;
    }

    return row;
}

int read_gff(char *path, size_t *nrows, Row **rows) {
    char buffer[BUFFER_SIZE];
    gzFile gff = gzopen(path, "rb");
    if (!gff) {
        log_error("Failed to open file: %s", path); 
        exit(-1);
    }

    (*rows) = (Row *) malloc(sizeof(Row));
    (*nrows) = 0;
    size_t size = 1;
    Row row;
    while (gzgets(gff, buffer, sizeof(buffer)) > 0) {
        char *s = buffer;
        if (*s == '#') continue;

        row = parse_row(s);

        (*nrows)++;
        if ((*nrows) > size) {
            size *= 2;
            (*rows) = (Row *) realloc(*rows, size * sizeof(Row));
        }
        (*rows)[(*nrows)-1] = row;
    }
    gzclose(gff);

    (*rows) = (Row *) realloc(*rows, (*nrows) * sizeof(Row));

    return 0;
}

void destroy_row(Row *row) {
    free(row->name);
    free(row->type);
    free(row->id);
}

void destroy_table(Table *table) {
    for (int i = 0; i < table->nrows; i++) {
        destroy_row(&table->rows[i]);
    }
    free(table->rows);
}

int get_transcripts(size_t nrows, Row *rows, size_t *argc, Transcript *argv[]) {
    Transcript *t = NULL;
    size_t ntranscripts = 0;
    size_t size = 1;
    Transcript *transcripts = (Transcript *) malloc(sizeof(Transcript));

    for (int i = 0; i < nrows; i++) {
        Row *row = &rows[i];

        // Only interested in start and stop codons and coding sequences
        if (strcmp(row->type, "CDS") != 0 && strcmp(row->type, "start_codon") != 0 && strcmp(row->type, "stop_codon") != 0) continue;

        // Find transcript if it exists
        for (int j = ntranscripts-1; j >= 0; j--) {
            if (strcmp(transcripts[j].id, row->id) == 0) {
                t = &(transcripts[j]);
                break;
            }
        }

        // Create transcript if it doesnt exist
        if (t == NULL) {
            ntranscripts++;
            if (ntranscripts > size) {
                size *= 2;
                transcripts = (Transcript *) realloc(transcripts, size * sizeof(Transcript));
            }

            transcripts[ntranscripts-1] = (Transcript) {strdup(row->id), strdup(row->name), 0, 0, row->dir, 0, NULL};
            t = &transcripts[ntranscripts-1];
            t->cds = (CDS *) malloc(sizeof(CDS));
        }

        if (strcmp(row->type, "CDS") == 0) {
            t->ncds++;
            t->cds = (CDS *) realloc(t->cds, t->ncds * sizeof(CDS));
            CDS cds = {row->start, row->end, row->phase};
            t->cds[t->ncds-1] = cds;
        } else if (strcmp(row->type, "start_codon") == 0) {
            t->start_codon = row->dir == FORWARD ? row->start : row->end;
        } else { // stop_codon
            t->end_codon = row->dir == FORWARD ? row->end : row->start;
        }

        t = NULL;
    }

    transcripts = (Transcript *) realloc(transcripts, ntranscripts * sizeof(Transcript));
    *argv = transcripts;
    *argc = ntranscripts;

    return 0;
}

void destroy_transcript(Transcript *t) {
    free(t->cds);
    free(t->region);
    free(t->id);
}

int clean_transcripts(size_t *ntranscripts, Transcript **transcripts) {
    size_t read, write = 0;

    Transcript *t;
    for (read = 0; read < *ntranscripts; read++) {
        t = &((*transcripts)[read]);
        if (t->start_codon == 0 || t->end_codon == 0) {
            destroy_transcript(t);
            continue; 
        }

        // Skip CDS that lie outside transcript regions
        int cds_w = 0;
        int total_cds_length = 0;
        for (int cds_r = 0; cds_r < t->ncds; cds_r++) { 
            if (t->direction == FORWARD && !(t->start_codon <= t->cds[cds_r].start && t->end_codon >= t->cds[cds_r].end)) continue;
            if (t->direction == BACKWARD && !(t->end_codon <= t->cds[cds_r].start && t->start_codon >= t->cds[cds_r].end)) continue;

            total_cds_length += t->cds[cds_r].end - (t->cds[cds_r].start-1);
            t->cds[cds_w++] = t->cds[cds_r];
        }
        t->ncds = cds_w;
        t->cds = (CDS *) realloc(t->cds, cds_w * sizeof(CDS));

        if (t->ncds > 0 && total_cds_length % 3 == 0) (*transcripts)[write++] = *t;
        else destroy_transcript(t);
    }

    *transcripts = (Transcript *) realloc(*transcripts, write * sizeof(Transcript));
    *ntranscripts = write;

    return 0;
}

int parse_gff(char *path, size_t *ntranscripts, Transcript **transcripts) {
    size_t nrows = 0;
    Row *rows;

    read_gff(path, &nrows, &rows);
    get_transcripts(nrows, rows, ntranscripts, transcripts);

    // Free rows
    for (int i = 0; i < nrows; i++) {
        destroy_row(&rows[i]);
    }
    free(rows);

    clean_transcripts(ntranscripts, transcripts);

    return 0;
}

void parse_consequences(char *csq, size_t *argc, char ***argv[]) {
    char *s = strtok(csq, ",");

    char ***consequences = (char ***) malloc(sizeof(char **));
    size_t nrows = 0, size = 1;
    while (s != NULL) {
        char **fields = (char **) malloc(N_CSQ_FIELDS * sizeof(char *));

        int i, len = 0;
        for(i = 0; s[len] && s[len] != '\0'; i++) {
            len = 0;
            while (s[len] && s[len] != '|' && s[len] != '\0') len++; // Find delimiter
            s[len] = '\0';

            fields[i] = strdup(s);


            s += len + 1;
        }

        if (i == N_CSQ_FIELDS-2) fields[++i] = ""; // If string ended after delimiter, set last value to empty string

        nrows++;
        if (nrows > size) {
            size *= 2;
            consequences = realloc(consequences, size * sizeof(char **));
        }
        consequences[nrows-1] = fields;

        s = strtok(NULL, ",");
    }
    free(s);

    consequences = realloc(consequences, nrows * sizeof(char **));
    *argc = nrows;
    *argv = consequences;
}

int main(int argc, char *argv[]) {

    // Parse gff file
    char *gff_path = argv[1];
    Transcript *transcripts = NULL;
    size_t ntranscripts = 0;
    // parse_gff(gff_path, &ntranscripts, &transcripts);


    gff_t *gff = gff_init(gff_path);
    gff_parse(gff);

    printf("Parsed succesfully\n");
    regidx_t *tscripts = (regidx_t *) gff_get(gff, idx_tscript);
//    regidx_t *idx = regidx_init(gff_path,NULL,NULL,0,NULL);
    regitr_t *itr = regitr_init(tscripts);
    
    if (regidx_overlap(tscripts, "chr1", 11868, 11868, itr)) printf("Overlap found\n");
    for (int i = 0; i < 25 && regitr_overlap(itr); i++) {
        printf("Transcript: chr=%s  beg=%d  end=%d\n", itr->seq, itr->beg+1, itr->end+1);
    }
    
    regitr_destroy(itr);

    // regidx_t *cds = (regidx_t *) gff_get(gff, idx_cds);
    // itr = regitr_init(cds);
    //
    // for (int i = 0; i < 25 && regitr_loop(itr); i++) {
    //     printf("CDS: chr=%s  beg=%d  end=%d\n", itr->seq, itr->beg+1, itr->end+1);
    // }
    //
    // regitr_destroy(itr);

    gff_destroy(gff);
    // for (int i = 0; i < ntranscripts; i++) {
    //     Transcript *t = &transcripts[i];
    //
    //     int total_size = 0;
    //     for (int j = 0; j < t->ncds; j++) {
    //         total_size += t->cds[j].end - (t->cds[j].start-1);
    //     }
    //
    //     printf("%s: %d\t%d\n", t->id, total_size, total_size % 3);
    // }

    // char *vcf_path = argv[2];
    // htsFile *fp_in = hts_open(vcf_path, "r");
    // if (!fp_in) {
    //     log_error("Could not open vcf file: %s", vcf_path);
    //     return 1;
    // }
    //
    // bcf_hdr_t *hdr = bcf_hdr_read(fp_in);
    // int csq_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "CSQ");
    // if (csq_id < 0) {
    //     log_error("VCF contains no INFO field CSQ, is it VEP annotated?");
    //     return 1;
    // }
    //
    //
    // bcf1_t *r = bcf_init();
    //
    // while (bcf_read(fp_in, hdr, r) >= 0) {
    //     if (!r) continue;
    //     
    //     bcf_unpack(r, BCF_UN_ALL);
    //
    //     char *csq_values = NULL;
    //     int n_csq_values = 0;
    //     
    //     char ***consequences;
    //     size_t nconsequences;
    //
    //     if (bcf_get_info_string(hdr, r, "CSQ", &csq_values, &n_csq_values) > 0) {
    //         parse_consequences(csq_values, &nconsequences, &consequences);
    //     }
    //
    //     for (int i = 0; i < nconsequences; i++) {
    //         for (int j = 0; j < N_CSQ_FIELDS; j++) {
    //             free(consequences[i][j]);
    //         }
    //         free(consequences[i]);
    //     }
    //     // Free the allocated memory
    //     // if (csq_values) {
    //     //     for (int i = 0; i < n_csq_values; i++) {
    //     //         free(csq_values[i]);
    //     //     }
    //     //     free(csq_values);
    //     // }
    // }


    // Free transcripts from gff file
    // for (int i = 0; i < ntranscripts; i++) {
    //     destroy_transcript(&(transcripts[i]));
    // }
    // free(transcripts);




    //char *fasta = argv[1];
    // faidx_t *fai = fai_load(fasta);


}

