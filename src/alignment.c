#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alignment.h"

alignment* newAlignment (int tid, int mm, int fraglen, int pos){
    alignment* node = (alignment*)malloc(sizeof(alignment));
    node->tid=tid;
    node->mm=mm;
    node->fraglen = fraglen;
    node->pos = pos; //start position on the tid, zero-based.
    node->next=NULL;
    return(node);
}

alignment_list* newAlignmentList(void){
   alignment_list* list = (alignment_list*)malloc(sizeof(alignment_list));
   list->first=NULL;
   list->last=NULL;
   list->size=0;
   return(list);
}


// if alignment is added, return 1, otherwise, return 0.
// keep the current_min_mm so it filters the best alignments only, in terms of mm.
// note: size_limit filtering cannot be done here, because the size_limit should be applied to best alignments in terms of mm.
// same for mmstr discrepancy filtering.
char add_alignment_to_list (alignment_list* list, alignment* new_alignment, int* current_min_mm)
{
   alignment* p;
   if(new_alignment==NULL) return 0;

   //fprintf(stderr,"new alignment: tid=%d, mm=%d, mmstr=%s, fraglen=%d, pos=%d\n",new_alignment->tid, new_alignment->mm, new_alignment->mmstr, new_alignment->fraglen, new_alignment->pos); //DEBUGGING

   // the following five lines is just for the exceptional case in which the original fastq file contains duplicate read id's that happens to occur in a row in the bowtie output file. In this case, the mapped tid and position must be identical and the duplicate is removed. Otherwise, it will be treated as an internal repeat, which causes an inaccurate estimation of fragment lenght effect because the corresponding EUMA value is zero.
   p=list->first;
   while(p!=NULL) {
     if(new_alignment->tid==p->tid && new_alignment->pos==p->pos && new_alignment->fraglen==p->fraglen) { delete_alignment(new_alignment); return 0; }
     p=p->next;
   }

   if(new_alignment->mm > (*current_min_mm)) { delete_alignment(new_alignment); return 0; } // keep only the best alignment in terms of mm.
   else if(new_alignment->mm < (*current_min_mm)) {
      delete_alignment_list(list);
      (*current_min_mm) = new_alignment->mm;
   }

   if(list->first==NULL) { 
      list->first = new_alignment;
      list->last = new_alignment;
      list->size = 1;
   }
   else {
      list->last->next = new_alignment;
      list->last = new_alignment;
      list->size++;
   }
   return 1;
}

void delete_alignment(alignment *aln){
   free(aln);
}

void delete_alignment_list(alignment_list *list){
   alignment *p,*q;
   if(list->first!=NULL){
     p=list->first;
     list->first=NULL;
     list->last=NULL;
     while(p!=NULL){
        q=p->next;
        delete_alignment(p);
        p=q;
     }
     list->size=0;
   }
}



// Returns 1 if discrepant fraglen exists. Otherwise returns 0.
// If the alignment list contains no element, also returns 0.
char check_fraglen_discrepancy(alignment_list* list){
  if(list->first==NULL) return 0;
  alignment* p=list->first;
  int fraglen1=p->fraglen;
  p=p->next;
  while(p!=NULL) {
    if(p->fraglen!=fraglen1) return 1; // discrepant fraglen
    p=p->next;
  }
  return 0;
}



// related functions

int parse_mmstr(char* mmstr){
     int mm=0,i;
     if(strlen(mmstr)>0) mm++;
     for(i=0;i<=strlen(mmstr);i++){
        if(mmstr[i]==',') mm++;
     }
     return(mm);
}


//if matching, returns 1, if not matching, returns 0.
// The two id's must be in format of xxx/1 and xxx/2.
char check_mate_readid_matching(char* read_id1, char* read_id2)
{
   int i,slen;
   if(strlen(read_id1)!=strlen(read_id2)) return 0;
   slen = strlen(read_id1);
   if(read_id1[slen-2]=='/' && read_id2[slen-2]=='/' && ((read_id1[slen-1]=='1' && read_id2[slen-1]=='2') || (read_id1[slen-1]=='2' && read_id2[slen-1]=='1') && strncmp(read_id1,read_id2,slen-2)==0)) return(slen-2);
   else {
     for(i=0;i<slen;i++){
       if(read_id1[i]==' ' && read_id2[i]==' ') return i;  // for parsing  newer Casava format (must be same up to the first blank)
       if(read_id1[i] != read_id2[i]) return 0;
     }
     return(slen);  // just the smae read name without '/1' or '/2'
   }
}

