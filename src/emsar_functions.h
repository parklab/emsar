#include "sam.h"

static char NUMBERS_STR[10]="0123456789";

void read_BAM_SE (char*, char, char);
alignment* convert_bam_alignment_2_alignment (bam1_t*, char*, char);
int parse_SAM_mmstr(char*);

