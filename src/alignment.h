#ifndef ALIGNMENT_H
#define ALIGNMENT_H


typedef struct alignment_node
{
     struct alignment_node* next;
     int tid;
     int mm;
     int fraglen;
     int pos;
} alignment;

typedef struct alignment_node_list
{
     alignment* first;
     alignment* last;
     int size;
} alignment_list;



alignment* newAlignment(int, int, int, int);
alignment_list* newAlignmentList(void);
char add_alignment_to_list (alignment_list*, alignment*, int*);
void delete_alignment_list(alignment_list*);
void delete_alignment(alignment*);
char check_fraglen_discrepancy(alignment_list*);

// these functions are related but not necessarily connected to the alignment or alignment_list structures.
char check_mate_readid_matching(char*, char*);
int parse_mmstr(char*);


#endif

