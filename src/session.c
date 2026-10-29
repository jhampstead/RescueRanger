#include "session.h"

void Session_close(Session *s) {
    if (s->scores_tsv) fclose(s->scores_tsv);
    if (s->fp_out)     hts_close(s->fp_out);
    if (s->fp_in)      hts_close(s->fp_in);
    if (s->vep_fmt)    destroy_hdr_fmt(s->vep_fmt);
    if (s->hdr)        bcf_hdr_destroy(s->hdr);
    if (s->reference)  ReferenceScores_close(s->reference); /* Before the index it reads the sequence with */
    if (s->fai)        fai_destroy(s->fai);
    if (s->gff)        gff_destroy(s->gff);
    if (s->models)     destroy_models(s->models);
    *s = (Session){0};
}
