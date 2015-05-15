#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#define SINGLE_LINE_MAX 500000   /* maximum length of transcript */
#define SINGLE_HEADER_MAX 1000   /* maximum length of header */
#define ISOFORM_MAX 500000   /* maximum number of isoforms */
#define FILENAMEMAX 1000  /* maximum length of output file path */
#define INIT_FASTA_SIZE 1000000 /* length of concatenated sequence */
#define FASTA_SIZE_ADD 1000000 /* length of concatenated sequence to be added each time more memory is needed */
#define MAX_HEADER_PREFIX_LEN 100

void generate_reads(char*, int, int, int,char, char*, char*, char*, char);
char revcomp(char c);
char checkvalid(char*,int);
