/*  error() and error_errno() of bcftools (bcftools.c), which the copy of the GFF reader calls. Without them the calls went to
    error() of the C library, which takes other arguments. */

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bcftools.h"

void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

void error_errno(const char *format, ...)
{
    va_list ap;
    int eno = errno;
    if ( format )
    {
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
    }
    if ( eno )
        fprintf(stderr, "%s%s\n", format ? ": " : "", strerror(eno));
    else
        fprintf(stderr, "\n");
    exit(-1);
}
