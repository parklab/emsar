#include "emsar.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "bool.h"
#include "alignment.h"
#include "sam.h"

pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;


/*------- global parameter setting -------*/

char set_library_strand_type(char* strand_type_str,char pe){
  if(strcmp(strand_type_str,"ns")==0) { library_strand_type=0; return 0; }
  else if(strcmp(strand_type_str,"ssf")==0 & pe==0) { library_strand_type='+'; return 0; }
  else if(strcmp(strand_type_str,"ssr")==0 & pe==0) { library_strand_type='-'; return 0; }
  else if(strcmp(strand_type_str,"ssfr")==0 & pe==1) { library_strand_type='+'; return 0; }
  else if(strcmp(strand_type_str,"ssrf")==0 & pe==1) { library_strand_type='-'; return 0; }
}





/*------- reading reference sequences & determining suffix array size -------*/

/* read raw fasta file (can be multi-line) and store it as a concatenated seq_array and creates an index_table file. */
void read_raw_fasta(char* fastafilename, char option)
/* option : 'E' : Ensembl format, 'R' : Refseq format. */
{
  FILE* fastafile=0;
  int i,j,jr,jf,tid;
  int c,prev_c;
  char mode;
  char *header_str;
  char *parsed_header_str;
  int cumlI_size;
  int th;
  int nm_index,transcript_len;


  fastafile = fopen((const char*)fastafilename,"r");
  /*Check for validity of the file.*/
  if(fastafile == 0)
  {
     printf("can't open fasta file.\n");
     exit(1);
  }

  max_tid=0; // initialize
  j=0; // position on seq_array.
  mode='h'; // start with a header mode.
  prev_c=0;

  while((c = fgetc(fastafile)) != EOF) {
     if(prev_c==0 && c!='>') { fprintf(stderr,"ERROR: wrong fasta file format (at tid=%d)\n",max_tid); exit(1); } //firstline not starting with '>'
     if(c=='>') {
        mode='h'; // header
        if(j!=0)  j++; 
     }
     else if(prev_c=='\n' && mode=='h') mode='s'; //sequence

     if(mode=='h'){
       if(c=='\n'){
         max_tid++;
       } 
     }
     else {   // mode=='s'
       if(c!='\n' && c!=' ' && c!='\t')  j++;
     }
     prev_c=c;
  }
  max_tid--;
  fclose(fastafile);
  fasta_size= (j+1)*2+1;  // (total sequence length + number of tid's (delimiters) ) *2 + null character

  fastafile = fopen((const char*)fastafilename,"r");
  /*Check for validity of the file.*/
  if(fastafile == 0)
  {
     printf("can't open fasta file.\n");
     exit(1);
  }

  tid=0; 
  i=0; // letter position in header string.
  j=0; // position on seq_array.
  mode='h'; // start with a header mode.
  prev_c=0;

  seq_array=(char*)malloc(fasta_size*sizeof(char));
  if(seq_array==NULL) { fprintf(stderr,"Failed to allocate memory for seq_array.\n"); exit(1); }
  if(verbose_flag>0) fprintf(stdout,"initialized seq_array..\n");  //DEBUGGING

  initialize_tname_indextable();

  header_str = (char*)malloc(SINGLE_LINE_MAX*sizeof(char));
  parsed_header_str = (char*)malloc(SINGLE_LINE_MAX*sizeof(char));


  while((c = fgetc(fastafile)) != EOF) {
     if(prev_c==0 && c!='>') { fprintf(stderr,"ERROR: wrong fasta file format.\n"); exit(1); } //firstline not starting with '>'
     if(c=='>') {
        mode='h'; // header
        if(j!=0) { 
           seq_array[j]='@'; // delimiter (getting ready for new sequence)
           j++; 
        }
     }
     else if(prev_c=='\n' && mode=='h') mode='s'; //sequence

     if(mode=='h'){
       if(c=='\n'){
         header_str[i]='\0';
         if(option=='E') parse_ensembl_header(header_str,&parsed_header_str);
         else parse_refseq_header(header_str,&parsed_header_str);

         // store the index table internally
         insert_key(parsed_header_str,tid,tname_tree);  // tname -> tid hash
         IndexTable[tid] = (char*)malloc((strlen(parsed_header_str)+1)*sizeof(char)); // tid -> tname array
         strcpy(IndexTable[tid],parsed_header_str);
         //fprintf(stderr,"IndexTable[tid=%d]=%s\n",tid,IndexTable[tid]);
         tid++;
         if(tid>max_tid+1){fprintf(stderr,"tid larger than max_tid.\n"); exit(1); }
         i=0;
       } else if(c!='>') { header_str[i]=c; i++; }  
     }
     else {   // mode=='s'
       if(c!='\n' && c!=' ' && c!='\t') { 
          seq_array[j]=uc(c); 
          j++;
       }
     }
     prev_c=c;
  }
  fclose(fastafile);

  free(header_str); header_str=NULL;
  free(parsed_header_str); parsed_header_str=NULL;

  // rc
  seq_array[j]='$'; // delimiter between fw and rc.
  borderpos=j; /* separator between fw and rc */
  
  jr=j+1; 
  for(jf=j-1;jf>=0;jf--){
    seq_array[jr]=revcomp(seq_array[jf]); jr++;
    if(jr+1>=fasta_size) { fprintf(stderr,"seq_array written beyond max at jr=%d\n",jr); exit(1);}
  }
  seq_array[jr]='$';  // the concatenated sequence is in the format: f0@f1@f2@f3$r3@r2@r1@r0$
  seq_array[jr+1]='\0';
  seqlength=jr; // last '$' position.
  //seq_array is done
  //fprintf(stdout,"seqarray=\n%s\n",seq_array);  //DEBUGGING

  // building cuml & computing sfa_size
  cuml= (int*)malloc((max_tid+2)*sizeof(int));  /*max index of cuml is max_tid+1*/
  if(cuml==NULL) { fprintf(stderr,"Failed to allocate memory to cuml.\n"); exit(1);}
  cuml[0]=0;
  nm_index=1;
  for(i=0;i<borderpos;i++){
     if(seq_array[i]=='@') {
       cuml[nm_index]=i+1;  /* cuml[x] is the first position of mRNA x. */
       nm_index++;
     }
  }
  cuml[nm_index]=borderpos+1;

  /* cuml index */
  cumlI_size=(int)(borderpos/bin)+1;  // double-check if this is okay.
  cumlI=(INTPAIR*)malloc(cumlI_size*sizeof(INTPAIR));
  if(cumlI==NULL) { fprintf(stderr,"Failed to allocate memory to cumlI.\n"); exit(1);}
  i=0;
  for(j=0;j<cumlI_size;j++){
    th=bin*(j+1);
    if(j==0) cumlI[j].start=0;
    else cumlI[j].start=i-1;
    while(i<=max_tid+1 && cuml[i]<th) i++;
    cumlI[j].end=i-1;
  }

  if(verbose_flag>0) fprintf(stdout,"cumlI_size=%d, max_tid=%d, fasta_size=%d\n",cumlI_size,max_tid,fasta_size); //DEBUGGING
  //for(j=0;j<cumlI_size;j++){ fprintf(stdout," cumlI[%d].start=%d, cumlI[%d].end=%d\n",j,cumlI[j].start,j,cumlI[j].end); } // DEBUGGING
  //for(i=0;i<=max_tid+1;i++) { fprintf(stdout, " cuml[%d]=%d\n",i,cuml[i]); } // DEBUGGING
}


void initialize_tname_indextable(void){
  tname_tree=new_tree();
  IndexTable=(char**)malloc((max_tid+1)*sizeof(char*));
  if(IndexTable==NULL) { fprintf(stderr,"Failed to allocate memory for IndexTable.\n"); exit(1); }
}


void determine_sfa_size(void){
   int i, transcript_len;
   max_sfa_size=0;
   for(i=0;i<=max_tid;i++){
      /* determine sfa size here, so that mRNAs shorter than read length is not counted. */
      transcript_len = cuml[i+1]-cuml[i]; // this length includes '@' so it has to be larger than readlength, not identical to readlength.
      if(transcript_len > readlength) max_sfa_size+=transcript_len-readlength;
   }
   if(verbose_flag>0) fprintf(stdout,"max_sfa_size=%d\n",max_sfa_size); //DEBUGGING
}


void read_bowtie_get_readlength (char* bowtiefilename){
    FILE* bowtie_file=0;
    char line[SINGLE_LINE_MAX],line2[SINGLE_LINE_MAX];

    if(strlen(bowtiefilename)==0){ bowtie_file = stdin; }
    else {
      bowtie_file = fopen((const char*)bowtiefilename,"r");
      /*Check for validity of the file.*/
      if(bowtie_file == 0)
      {
         fprintf(stderr,"can't open bowtie file.\n");
         exit(1);
      }
    }

    readlength=-1;
    while (fgets(line, SINGLE_LINE_MAX, bowtie_file) != NULL && readlength==-1) {
       if(line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0'; /*chomp*/
       readlength = parse_bowtieline_get_readlength (line);
    }
    fclose(bowtie_file);
}

void read_BAM_get_readlength (char* bamfilename, char sam_or_bam){
    samfile_t* bamfile=0;
    bam1_t* bam_alignment;

    if(strlen(bamfilename)==0) strcpy(bamfilename,"-"); 
    bamfile = samopen(bamfilename,sam_or_bam=='s'?"r":"rb",0);
    /*Check for validity of the file.*/
    if(bamfile == 0)
    {
       fprintf(stderr,"can't open %s file.\n",sam_or_bam=='s'?"SAM":"BAM");
       exit(1);
    }

    readlength=-1;
    while(samread(bamfile,bam_alignment=bam_init1())>0 && readlength==-1){
       if(bam_alignment->core.tid==-1) {bam_destroy1(bam_alignment); continue; } // skip unaligned reads.
       readlength = bam_alignment->core.l_qseq;

       //free bam_alignment
       bam_destroy1(bam_alignment);
    }
    samclose(bamfile);
}




void read_bowtie_get_readlengths_se (char* bowtiefilename){
    FILE* bowtie_file=0;
    char line[SINGLE_LINE_MAX],line2[SINGLE_LINE_MAX];

    if(strlen(bowtiefilename)==0){ bowtie_file = stdin; }
    else {
      bowtie_file = fopen((const char*)bowtiefilename,"r");
      /*Check for validity of the file.*/
      if(bowtie_file == 0)
      {
         fprintf(stderr,"can't open bowtie file.\n");
         exit(1);
      }
    }

    readlength=-1;
    Fraglengths.max=0;  // some arbitrary small initial number
    Fraglengths.min=30000; // some arbitrary big initial number
    while (fgets(line, SINGLE_LINE_MAX, bowtie_file) != NULL) {
       if(line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0'; /*chomp*/
       readlength = parse_bowtieline_get_readlength (line);
       if(readlength<Fraglengths.min) Fraglengths.min=readlength;
       if(readlength>Fraglengths.max) Fraglengths.max=readlength;
    }
    fclose(bowtie_file);
    nFraglen = Fraglengths.max - Fraglengths.min +1;
    Readlengths.max=Fraglengths.max;
    Readlengths.min=Fraglengths.min;
}

void read_BAM_get_readlengths_se (char* bamfilename, char sam_or_bam){
    samfile_t* bamfile=0;
    bam1_t* bam_alignment;

    if(strlen(bamfilename)==0) strcpy(bamfilename,"-"); 
    bamfile = samopen(bamfilename,sam_or_bam=='s'?"r":"rb",0);
    /*Check for validity of the file.*/
    if(bamfile == 0)
    {
       fprintf(stderr,"can't open %s file.\n",sam_or_bam=='s'?"SAM":"BAM");
       exit(1);
    }

    readlength=-1;
    Fraglengths.max=0;  // some arbitrary small initial number
    Fraglengths.min=30000; // some arbitrary big initial number
    while(samread(bamfile,bam_alignment=bam_init1())>0){
       if(bam_alignment->core.tid==-1) {bam_destroy1(bam_alignment); continue; } // skip unaligned reads.
       readlength = bam_alignment->core.l_qseq;
       if(readlength<Fraglengths.min) Fraglengths.min=readlength;
       if(readlength>Fraglengths.max) Fraglengths.max=readlength;

       //free bam_alignment
       bam_destroy1(bam_alignment);
    }
    samclose(bamfile);
    nFraglen = Fraglengths.max - Fraglengths.min +1;
    Readlengths.max=Fraglengths.max;
    Readlengths.min=Fraglengths.min;
}

/*------- parsing an alignment file -------*/

void read_BAM_SE (char* bamfilename, char sam_or_bam, char mode){
    samfile_t* bamfile=0;
    bam1_t* bam_alignment;
    char *read_id,*prev_read_id;
    int i,j,readlen;
    alignment* new_alignment;
    alignment_list *alignmentlist = newAlignmentList();
    int current_min_mm=10000; /* some really big number */

    // Initialize
    if(mode==0){
      TotalReadCount=0;
      readlength=-1;
    }
    else if(mode==1){;}
    else if(mode==2){
      TotalReadCount=0;
    }

    if(strlen(bamfilename)==0) strcpy(bamfilename,"-"); 
    bamfile = samopen(bamfilename,sam_or_bam=='s'?"r":"rb",0);
    /*Check for validity of the file.*/
    if(bamfile == 0)
    {
       fprintf(stderr,"can't open %s file.\n",sam_or_bam=='s'?"SAM":"BAM");
       exit(1);
    }

    read_id=(char*)malloc(SINGLE_LINE_MAX*sizeof(char));
    if(read_id==NULL) { fprintf(stderr,"Failed to allocate memory to read_id\n"); exit(1); }
    prev_read_id=(char*)malloc(SINGLE_LINE_MAX*sizeof(char));
    if(prev_read_id==NULL) { fprintf(stderr,"Failed to allocate memory to prev_read_id\n"); exit(1); }

    strcpy(prev_read_id," ");

    while(samread(bamfile,bam_alignment=bam_init1())>0){
       if(bam_alignment->core.tid==-1) { bam_destroy1(bam_alignment); continue; } // skip unaligned reads.

       strcpy(read_id,bam1_qname(bam_alignment));
       new_alignment = convert_bam_alignment_2_alignment (bam_alignment, (char*)(bamfile->header->target_name[bam_alignment->core.tid]));

       //free bam_alignment
       bam_destroy1(bam_alignment);

       //PE and SE are treated the same below, except for the fraglen discrepancy filtering for PE.
       if(new_alignment!=NULL){
          if( strcmp(prev_read_id,read_id)==0 ) add_alignment_to_list(alignmentlist, new_alignment, &current_min_mm);
          else {
             if(strcmp(prev_read_id," ")!=0) {
                if(alignmentlist->size <= MAX_REPEAT) update_ReadCounts(alignmentlist);
                delete_alignment_list(alignmentlist);
             }
             current_min_mm=10000;
             add_alignment_to_list(alignmentlist, new_alignment, &current_min_mm);
          }
          strcpy(prev_read_id,read_id);
       }
    }
    if(alignmentlist->size <= MAX_REPEAT) update_ReadCounts(alignmentlist);
    delete_alignment_list(alignmentlist);

    samclose(bamfile);
    free(read_id); read_id=NULL;
    free(prev_read_id); prev_read_id=NULL;
    free(alignmentlist);
}


alignment* convert_bam_alignment_2_alignment (bam1_t* bam_alignment, char* tname){
   int tid,mm,fraglen,pos;
   char* SAM_mmstr;
   char strand;
  
   //tid
   tid = search_treehash(tname,tname_tree);
   if(tid==-1 || tid>max_tid){fprintf(stderr,"error: unexisting tid in the bowtie output file. Check bowtieout file.\n"); exit(1); }

   //strand
   strand = (bam_alignment->core.flag&0x10?'-':'+');
   if(library_strand_type!=0 && library_strand_type!=strand) return NULL; // filter strand here

   //mm
   SAM_mmstr = bam_aux2Z( bam_aux_get(bam_alignment, "MD") );
   mm = parse_SAM_mmstr(SAM_mmstr);

   //fraglen
   fraglen = bam_alignment->core.l_qseq;

   //pos
   pos = bam_alignment->core.pos;

   return ( newAlignment(tid,mm,fraglen,pos) );
}


int parse_SAM_mmstr(char* SAM_mmstr){
  int mm=0,i;
  for(i=0;i<strlen(SAM_mmstr);i++){
    if(strchr(NUMBERS_STR,SAM_mmstr[i])==NULL) mm++;  // any non-numeric character adds to mm.
  } 
  return mm;
}

alignment* convert_bam_alignment_2_alignment_PE (bam1_t* bam_alignment, bam1_t* bam_alignment2, char* tname){
   bam1_t *b1,*b2;
   int tid,mm,fraglen,pos1,pos2,pos;
   char *SAM_mmstr1, *SAM_mmstr2;
   char SAM_mmstr[SINGLE_LINE_MAX];
   char strand1,strand2;

   //tid
   tid = search_treehash(tname,tname_tree);
   if(tid==-1 || tid>max_tid){fprintf(stderr,"error: unexisting tid in the bowtie output file. Check bam/sam file.\n"); exit(1); }

   //readlen (assuming mate1 and mate2 and all reads have the same readlen)
   if(readlength==-1) readlength = bam_alignment->core.l_qseq;
   if(readlength!=bam_alignment->core.l_qseq || readlength!=bam_alignment2->core.l_qseq) { fprintf(stderr,"Error: Paired-end data with variable read length is not supported. Check your bam/sam file.\n"); exit(1); }  

   // mate info (mate1=b1, mate2=b2)
   if((bam_alignment->core.flag & 0x40) && (bam_alignment2->core.flag & 0x80)) { b1=bam_alignment; b2=bam_alignment2; }
   else if((bam_alignment2->core.flag & 0x40) && (bam_alignment->core.flag & 0x80)) { b1=bam_alignment2; b2=bam_alignment; }
   else { fprintf(stderr, "error: mates are not grouped in the BAM/SAM file.\n"); exit(1); } 

   //mm
   SAM_mmstr1 = bam_aux2Z( bam_aux_get(b1, "MD") );
   SAM_mmstr2 = bam_aux2Z( bam_aux_get(b2, "MD") );
   mm = parse_SAM_mmstr(SAM_mmstr1) + parse_SAM_mmstr(SAM_mmstr2);

   //fraglen, pos  & strand
   pos1 = b1->core.pos;
   pos2 = b2->core.pos;
   strand1 = (b1->core.flag&0x10?'-':'+');
   strand2 = (b2->core.flag&0x10?'-':'+');
   if(pos2>pos1) { // mate1(f)...mate2(r) 
      fraglen = pos2-pos1+readlength; pos=pos1;
      if(library_strand_type=='-') return NULL;
      if(!(strand1=='+' && strand2=='-')) return NULL;
    }
    else { //mate2(f)...mate1(r)
      fraglen = pos1-pos2+readlength; pos=pos2;
      if(library_strand_type=='+') return NULL;
      if(!(strand1=='-' && strand2=='+')) return NULL;
    }

   return ( newAlignment(tid,mm,fraglen,pos) );

}


/* desired_strand is NULL (\0) if nonstrand-specific */
// sam_or_bam is 's' for sam and 'b' for bam.
void read_BAM_PE (char* bamfilename, char sam_or_bam, char mode){  
// mode : 0 (first bam), 1 (another bam, continued), 2 (another bam, new)
// frag_or_rsh : 0 : frag, 1: rsh
    samfile_t* bamfile=0;
    bam1_t *bam_alignment, *bam_alignment2;

    char *read_id,*prev_read_id;
    int i,j,readlen;
    alignment* new_alignment;
    alignment_list *alignmentlist = newAlignmentList();
    int current_min_mm=10000; /* some really big number */

    // Initialize
    if(mode==0){
      TotalReadCount=0; 
      readlength=-1;  
    } 
    else if(mode==1){;}
    else if(mode==2){
      TotalReadCount=0;
    }


    if(strlen(bamfilename)==0) strcpy(bamfilename,"-"); 
    bamfile = samopen(bamfilename,sam_or_bam=='s'?"r":"rb",0);
    /*Check for validity of the file.*/
    if(bamfile == 0)
    {
       fprintf(stderr,"can't open %s file.\n",sam_or_bam=='s'?"SAM":"BAM");
       exit(1);
    }

    read_id=(char*)malloc(SINGLE_LINE_MAX*sizeof(char));
    if(read_id==NULL) { fprintf(stderr,"Failed to allocate memory to read_id\n"); exit(1); }
    prev_read_id=(char*)malloc(SINGLE_LINE_MAX*sizeof(char));
    if(prev_read_id==NULL) { fprintf(stderr,"Failed to allocate memory to prev_read_id\n"); exit(1); }

    strcpy(prev_read_id," ");

    fprintf(stderr,"reading bam..\n"); //DEBUGGING
    while (samread(bamfile,bam_alignment=bam_init1())>0){
       if(bam_alignment->core.tid==-1) { bam_destroy1(bam_alignment); continue; } // skip unaligned reads.
      
       if(samread(bamfile,bam_alignment2=bam_init1())>0) { //mate in the second line
         strcpy(read_id,bam1_qname(bam_alignment));  // We assume that the two mates occur in the bam file one after another (stricter condition than simply sorted by read ID. We use the qname of the first occuring mate as the read name of the mate pair. So, if the two mates happen to have different qnames and their orders change, it could cause a problem.
         new_alignment = convert_bam_alignment_2_alignment_PE (bam_alignment, bam_alignment2, (char*)(bamfile->header->target_name[bam_alignment->core.tid]));
       }

       //free bam_alignment
       bam_destroy1(bam_alignment);
       bam_destroy1(bam_alignment2);

       //PE and SE are treated the same below, except for the fraglen discrepancy filtering for PE.
       
       if(new_alignment!=NULL){
          if( strcmp(prev_read_id,read_id)==0 ) add_alignment_to_list(alignmentlist, new_alignment, &current_min_mm);
          else {
             if(strcmp(prev_read_id," ")!=0) {
                if(alignmentlist->size <= MAX_REPEAT && check_fraglen_discrepancy(alignmentlist)==0) update_ReadCounts(alignmentlist); 
                delete_alignment_list(alignmentlist);
             }
             current_min_mm=10000;
             add_alignment_to_list(alignmentlist, new_alignment, &current_min_mm);
          }
          strcpy(prev_read_id,read_id);
       }
    }
    if(alignmentlist->size <= MAX_REPEAT && check_fraglen_discrepancy(alignmentlist)==0) update_ReadCounts(alignmentlist);
    delete_alignment_list(alignmentlist);

    samclose(bamfile);
    free(read_id); read_id=NULL;
    free(prev_read_id); prev_read_id=NULL;
    free(alignmentlist); alignmentlist=NULL;
}



alignment* parse_bowtieline (char* line, char** read_id){
       char tmpstr[SINGLE_LINE_MAX],tname[SINGLE_LINE_MAX],mmstr[SINGLE_LINE_MAX];
       int numtab,j,i,pos,fraglen,tid,mm;
       char strand;

       mmstr[0]=0; // some initialization

       /*The following parsing is a replacement of strtok which may cause memory leak. */
       numtab=0;j=0;
       for(i=0;i<=strlen(line);i++){
          if(line[i]=='\t'||line[i]=='\0') {
            numtab++;
            tmpstr[j]='\0';
            if(numtab==1) strcpy((*read_id),tmpstr);  // instead of returning the read_id, store it to the address passed as an argument
            else if(numtab==2) { 
               strand=tmpstr[0];
               if(library_strand_type!=0 && library_strand_type!=strand) return NULL;  //filter strand here.
            }
            else if(numtab==3) strcpy(tname,tmpstr);
            else if(numtab==4) pos = atoi(tmpstr);  
            else if(numtab==5) fraglen = strlen(tmpstr);
            else if(numtab==8) if(strlen(tmpstr)>0) strcpy(mmstr,tmpstr); else mmstr[0]='\0';
            j=0;
          }
          else { tmpstr[j]=line[i]; j++; }
       }
       if(numtab<7) { fprintf(stderr,"Error: input alignment file doesn't look like bowtieout file.\n"); exit(1); } 

       tid = search_treehash(tname,tname_tree);
       if(tid==-1 || tid>max_tid){ fprintf(stderr,"error: unexisting tid in the bowtie output file. Check bowtieout file.\n"); exit(1); }

       //mm
       mm=parse_mmstr(mmstr);

       return ( newAlignment(tid,mm,fraglen,pos) );
}



int parse_bowtieline_get_readlength (char* line){
       char tmpstr[SINGLE_LINE_MAX];
       int numtab,j,i,readlen;
       char strand;

       /*The following parsing is a replacement of strtok which may cause memory leak. */
       numtab=0;j=0;
       for(i=0;i<=strlen(line);i++){
          if(line[i]=='\t'||line[i]=='\0') {
            numtab++;
            tmpstr[j]='\0';
            if(numtab==5) readlen = strlen(tmpstr);
            j=0;
          }
          else { tmpstr[j]=line[i]; j++; }
       }
       if(numtab<7) { fprintf(stderr,"Error: input alignment file doesn't look like bowtieout file.\n"); exit(1); } 

       return ( readlen ); 
}

alignment* parse_bowtieline_PE (char* line1, char* line2, char** read_id){
       char tmpstr[SINGLE_LINE_MAX],tname1[SINGLE_LINE_MAX],tname2[SINGLE_LINE_MAX],*mmstr1,*mmstr2,read_id2[SINGLE_LINE_MAX],*tmpmmstr;
       int numtab,j,i,pos1,pos2,pos,tmppos,tid,mm,readlen1,readlen2,fraglen;
       char strand1,strand2,tmpstrand;
       int adjusted_read_id_len;
       char order_reversed;
       alignment* new_aln;

       mmstr1=(char*)malloc((strlen(line1)+strlen(line2))*sizeof(char));
       if(mmstr1==NULL) { fprintf(stderr,"Failed to allocate memory to mmstr1\n"); exit(1); }
       mmstr2=(char*)malloc(strlen(line2)*sizeof(char));
       if(mmstr2==NULL) { fprintf(stderr,"Failed to allocate memory to mmstr2\n"); exit(1); }

       /*The following parsing is a replacement of strtok which may cause memory leak. */
       numtab=0;j=0;
       for(i=0;i<=strlen(line1);i++){
          if(line1[i]=='\t'||line1[i]=='\0') {
            numtab++;
            tmpstr[j]='\0';
            if(numtab==1) strcpy((*read_id),tmpstr);  // instead of returning the read_id, store it to the address passed as an argument
            else if(numtab==2) strand1 = tmpstr[0];
            else if(numtab==3) strcpy(tname1,tmpstr);
            else if(numtab==4) pos1 = atoi(tmpstr);  
            else if(numtab==5) readlen1 = strlen(tmpstr);
            else if(numtab==8) if(strlen(tmpstr)>0) strcpy(mmstr1,tmpstr); else mmstr1[0]='\0';
            j=0;
          }
          else { tmpstr[j]=line1[i]; j++; }
       }
       if(numtab<7) { fprintf(stderr,"Error: input alignment file doesn't look like bowtieout file.\n"); exit(1); } 

       /*The following parsing is a replacement of strtok which may cause memory leak. */
       numtab=0;j=0;
       for(i=0;i<=strlen(line2);i++){
          if(line2[i]=='\t'||line2[i]=='\0') {
            numtab++;
            tmpstr[j]='\0';
            if(numtab==1) { 
               strcpy(read_id2,tmpstr); 
               if((adjusted_read_id_len = check_mate_readid_matching(*read_id,read_id2))==0) {fprintf(stderr,"Error: mate read ID's don't match. Check bowtie out format.\n"); exit(1); }
               if((*read_id)[strlen(*read_id)-1]==1) order_reversed=0; else order_reversed=1;
               (*read_id)[adjusted_read_id_len]='\0';  // remove '/1' or '/2' at the end or ' 1:N:0:GCCAAT/2' / ' 2:N:0:GCCAAT/2'.
            }
            else if(numtab==2) strand2 = tmpstr[0];
            else if(numtab==3) strcpy(tname2,tmpstr);
            else if(numtab==4) pos2 = atoi(tmpstr);  
            else if(numtab==5) readlen2 = strlen(tmpstr);
            else if(numtab==8) if(strlen(tmpstr)>0) strcpy(mmstr2,tmpstr); else mmstr2[0]='\0';
            j=0;
          }
          else { tmpstr[j]=line2[i]; j++; }
       }
       if(numtab<7) { fprintf(stderr,"Error: input alignment file doesn't look like bowtieout file.\n"); exit(1); } 


       if(strcmp(tname1,tname2)!=0) return NULL;  // Two mates aren't mapped to the same RNA, so don't include this alignment.

       //check read length (PE version doesn't support variable read lengths)
       if(readlength==-1) readlength=readlen1;
       if(readlength!=readlen1 || readlength!=readlen2) { fprintf(stderr,"Error: Paired-end data with variable read length is not supported. Check your bowtieout file.\n"); exit(1); }

       //convert tid
       tid = search_treehash(tname1,tname_tree);
       //fprintf(stderr,"reading bowtie file, tid=%d\n",tid);  //DEBUGGING
       if(tid==-1 || tid>max_tid){ fprintf(stderr,"error: unexisting tid in the bowtie output file. Check bowtieout file.\n"); exit(1); }

       if(order_reversed) { 
         tmppos = pos1; pos1=pos2; pos2=tmppos;
         tmpstrand = strand1; strand1=strand2; strand2=tmpstrand;
         tmpmmstr = mmstr1; mmstr1= mmstr2; mmstr2=tmpmmstr;
       }

       //mm
       mm=parse_mmstr(mmstr1)+parse_mmstr(mmstr2);

       //fraglen, pos & strand checking
       if(pos2>pos1) { // mate1(f)...mate2(r) 
         fraglen = pos2-pos1+readlength; pos=pos1;
         if(library_strand_type=='-') return NULL;
         if(!(strand1=='+' && strand2=='-')) return NULL;
       }
       else { //mate2(f)...mate1(r)
         fraglen = pos1-pos2+readlength; pos=pos2;
         if(library_strand_type=='+') return NULL;
         if(!(strand1=='-' && strand2=='+')) return NULL;
       }

       new_aln = newAlignment(tid,mm,fraglen,pos);

       free(mmstr1); free(mmstr2);
       return ( new_aln );
}



void read_bowtie_SE (char* bowtiefilename, char mode){
    FILE* bowtie_file=0;
    char line[SINGLE_LINE_MAX],line2[SINGLE_LINE_MAX],*read_id,*prev_read_id;
    int i,j,readlen;
    alignment* new_alignment;
    alignment_list *alignmentlist = newAlignmentList();
    int current_min_mm=10000; /* some really big number */

    // Initialize
    if(mode==0){
      TotalReadCount=0;
      readlength=-1;
    }
    else if(mode==1){;}
    else if(mode==2){
      TotalReadCount=0;
    }

    if(strlen(bowtiefilename)==0){ bowtie_file = stdin; }
    else {
      bowtie_file = fopen((const char*)bowtiefilename,"r");
      /*Check for validity of the file.*/
      if(bowtie_file == 0)
      {
         fprintf(stderr,"can't open bowtie file.\n");
         exit(1);
      }
    }

    read_id=(char*)malloc(SINGLE_LINE_MAX*sizeof(char));
    if(read_id==NULL) { fprintf(stderr,"Failed to allocate memory to read_id"); exit(1); }
    prev_read_id=(char*)malloc(SINGLE_LINE_MAX*sizeof(char));
    if(prev_read_id==NULL) { fprintf(stderr,"Failed to allocate memory to prev_read_id"); exit(1); }

    strcpy(prev_read_id," ");
    while (fgets(line, SINGLE_LINE_MAX, bowtie_file) != NULL) {
       if(line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0'; /*chomp*/

       new_alignment = parse_bowtieline (line, &read_id);

       //PE and SE are treated the same below, except for the fraglen discrepancy filtering for PE.
       if(new_alignment!=NULL){
          if( strcmp(prev_read_id,read_id)==0 ) add_alignment_to_list(alignmentlist, new_alignment, &current_min_mm);
          else {
             if(strcmp(prev_read_id," ")!=0) {
                if(alignmentlist->size <= MAX_REPEAT) update_ReadCounts(alignmentlist);
                delete_alignment_list(alignmentlist);
             }
             current_min_mm=10000;
             add_alignment_to_list(alignmentlist, new_alignment, &current_min_mm);
          }
          strcpy(prev_read_id,read_id);
       }
    }
    if(alignmentlist->size <= MAX_REPEAT) update_ReadCounts(alignmentlist);
    delete_alignment_list(alignmentlist);

    fclose(bowtie_file);
    free(read_id); read_id=NULL;
    free(prev_read_id); prev_read_id=NULL;
    free(alignmentlist); alignmentlist=NULL;
}


void read_bowtie_PE (char* bowtiefilename, char mode){
    FILE* bowtie_file=0;
    char line[SINGLE_LINE_MAX],line2[SINGLE_LINE_MAX],*read_id,*prev_read_id;
    int i,j,readlen;
    alignment* new_alignment;
    alignment_list *alignmentlist = newAlignmentList();
    int current_min_mm=10000; /* some really big number */

    // Initialize
    if(mode==0){
      TotalReadCount=0;
      readlength=-1;
    }
    else if(mode==1){;}
    else if(mode==2){
      TotalReadCount=0;
    }

    if(strlen(bowtiefilename)==0){ bowtie_file = stdin; }
    else {
      bowtie_file = fopen((const char*)bowtiefilename,"r");
      /*Check for validity of the file.*/
      if(bowtie_file == 0)
      {
         fprintf(stderr,"can't open bowtie file.\n");
         exit(1);
      }
    }

    read_id=(char*)malloc(SINGLE_LINE_MAX*sizeof(char));
    if(read_id==NULL) { fprintf(stderr,"Failed to allocate memory to read_id\n"); exit(1); }
    prev_read_id=(char*)malloc(SINGLE_LINE_MAX*sizeof(char));
    if(prev_read_id==NULL) { fprintf(stderr,"Failed to allocate memory to prev_read_id\n"); exit(1); }

    strcpy(prev_read_id," ");
    while (fgets(line, SINGLE_LINE_MAX, bowtie_file) != NULL) {
       if(line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0'; /*chomp*/

       //read mate2 line
       fgets(line2, SINGLE_LINE_MAX, bowtie_file);
       if(line2[strlen(line2)-1]=='\n') line2[strlen(line2)-1]='\0'; /*chomp*/
       new_alignment = parse_bowtieline_PE (line, line2, &read_id);

       //PE and SE are treated the same below, except for the fraglen discrepancy filtering for PE.
       if(new_alignment!=NULL){
          if( strcmp(prev_read_id,read_id)==0 ) add_alignment_to_list(alignmentlist, new_alignment, &current_min_mm);
          else {
             if(strcmp(prev_read_id," ")!=0) {
                if(alignmentlist->size <= MAX_REPEAT && check_fraglen_discrepancy(alignmentlist)==0) update_ReadCounts(alignmentlist);
                delete_alignment_list(alignmentlist);
             }
             current_min_mm=10000;
             add_alignment_to_list(alignmentlist, new_alignment, &current_min_mm);

          }
          strcpy(prev_read_id,read_id);
       }
    }
    if(alignmentlist->size <= MAX_REPEAT && check_fraglen_discrepancy(alignmentlist)==0) update_ReadCounts(alignmentlist);
    delete_alignment_list(alignmentlist);

    fclose(bowtie_file);
    free(read_id); read_id=NULL;
    free(prev_read_id); prev_read_id=NULL;
    free(alignmentlist);
}

void update_ReadCounts(alignment_list* list){
   int fraglength;
   int i,j,ti2,tsize,nopush;
   int *tarray;
   alignment* p;
   int tlen,d3,dmid5,dmid3;
   
   if(list->first==NULL) {fprintf(stderr,"Error: NULL alignment list.\n"); exit(1); }

   fraglength = list->first->fraglen;  // here we assume that all alignments have the same fraglength. This may not be true for PE, but such cases are filtered out.

   if(fraglength<=Max_Fraglength && fraglength>=Min_Fraglength) {  /* filter by fraglength */
      if(list->size==1) { 
         p=list->first;
         if(posmodel==1){
            tlen = (cuml[p->tid+1]-cuml[p->tid]-1);
            
            if(p->pos<perpos_freq_len) perpos_freq_5[p->pos]++;
            d3 = tlen - (p->pos + p->fraglen - 1); //length(tid) - (3end of read)
            if(d3<perpos_freq_len) perpos_freq_3[d3]++; 
            if(tlen<perpos_freq_len) {  // for normalizing the counts by available length
               for(i=tlen;i<perpos_freq_len;i++) {
                 perpos_unavail_freq_5[i]++;
                 perpos_unavail_freq_3[i]++;
               }
            }
            

            /*
            dmid5 = p->pos + ( p->fraglen-1 ) / 2;  // relative position of center of read with respect to 5'end of transcript
            dmid3 = tlen - dmid5; // relative position of center of read with respect to 3'end of transcript
            if(dmid5<perpos_freq_len) perpos_freq_5[dmid5]++;
            if(dmid3<perpos_freq_len) perpos_freq_3[dmid3]++; 
            if(tlen<perpos_freq_len) {  // for normalizing the counts by available length
               for(i=tlen;i<perpos_freq_len;i++) {
                 perpos_unavail_freq_5[i]++;
                 perpos_unavail_freq_3[i]++;
               }
            }
            */
          }

          if((*update_rshbucket_single_PTR)(p->tid,'r',fraglength,NULL)) {
              rsh_size_single++; 
          }
      }
      else {
        tarray = (int*)malloc(list->size*sizeof(int));
        p=list->first; tsize=0; nopush=0;
        while(p!=NULL) { 
          for(j=0;j<tsize;j++) {
             if(tarray[j] >= p->tid) {          
                for(ti2=tsize-1; ti2>=j; ti2--) { tarray[ti2+1] = tarray[ti2]; }
                tarray[j]=p->tid;
                tsize++;
                nopush=1;
                break;
             }
          }
          if(nopush==0) {  /* push to the array */
             tarray[tsize]=p->tid;
             tsize++;
          } 
          p=p->next; nopush=0;
        }

        if(posmodel==1){
          p=list->first;
          while(p!=NULL) {
            tlen = cuml[p->tid+1]-cuml[p->tid]-1;

            
            if(p->pos<perpos_freq_len) perpos_freq_5[p->pos]+=1.0/(double)tsize;
            d3 = tlen - (p->pos + p->fraglen - 1); //length(tid) - (3end of read)
            if(d3<perpos_freq_len) perpos_freq_3[d3]+=1.0/(double)tsize;
            if(tlen<perpos_freq_len) {  // for normalizing the counts by available length
               for(i=tlen;i<perpos_freq_len;i++) {
                 perpos_unavail_freq_5[i]+=1.0/(double)tsize;
                 perpos_unavail_freq_3[i]+=1.0/(double)tsize;
               }
            }
            
            /*
            dmid5 = p->pos + ( p->fraglen-1 ) / 2;  // relative position of center of read with respect to 5'end of transcript
            dmid3 = tlen - dmid5; // relative position of center of read with respect to 3'end of transcript
            if(dmid5<perpos_freq_len) perpos_freq_5[dmid5]+=1.0/(double)tsize;
            if(dmid3<perpos_freq_len) perpos_freq_3[dmid3]+=1.0/(double)tsize;
            if(tlen<perpos_freq_len) {  // for normalizing the counts by available length
               for(i=tlen;i<perpos_freq_len;i++) {
                 perpos_unavail_freq_5[i]+=1.0/(double)tsize;
                 perpos_unavail_freq_3[i]+=1.0/(double)tsize;
               }
            }
            */
            p=p->next;
          }
        }
        

        if((*update_rshbucket_PTR)(tsize,tarray,'r',fraglength,NULL)) rsh_size++;
        free(tarray); tarray=NULL;
      }
      FraglengthCounts[fraglength]++;
      TotalReadCount++;
   }
}



/*------- suffix array initializing & sorting -------*/ 

void initialize_suffixarray_SS_4(char* tag)
{
           int i,t,k,length_t,s_index;
           int taglen=strlen(tag);
           int numtags;

           if(taglen==1) numtags=4;
           else if(taglen==2) numtags=16;
           else numtags=64;

           sfa_size = max_sfa_size/numtags *2;  // *2 to avoid reallocating more memory as much as possible
           sfa = (int*)malloc(sfa_size*sizeof(int));
           if(sfa==NULL) { fprintf(stderr,"Failed to allocate memory to suffix array.\n"); exit(1); }
           s_index=0;
           for(t=0;t<=max_tid;t++){
              length_t=cuml[t+1]-cuml[t];
              for(k=0;k<length_t-readlength;k++){
                if(strncmp(seq_array+cuml[t]+k,tag,taglen)==0){
                  if(is_noncanonical(seq_array+cuml[t]+k,readlength)) continue; // include only ACGT
                  sfa[s_index]=cuml[t]+k;
                  s_index++;
                  if(s_index>=sfa_size) { 
                    sfa_size+=max_sfa_size/numtags; 
                    sfa = realloc(sfa,sfa_size*sizeof(int));
                    if(sfa==NULL) { fprintf(stderr,"Failed to reallocate memory to suffix array at tag=%s, tid=%d, s_index=%d\n",tag,t,s_index); exit(1); }
                  }
                }
              }
           }
           if(s_index>0){
             sfa_size=s_index;
             sfa = realloc(sfa,sfa_size*sizeof(int));
             if(sfa==NULL) { fprintf(stderr,"Failed to reallocate memory to suffix array at tag=%s in the end of the cycle\n",tag); exit(1); }
             local_sfa_start=0; local_sfa_end = sfa_size-1;
           }
           else { sfa_size=0; local_sfa_start=0; local_sfa_end = -1; free(sfa); sfa=NULL; }
}

void initialize_suffixarray_NS_5(char* tag)
{
           int i,t,s_index,flippos;
           int taglen=strlen(tag);
           int numtags;

           if(taglen==1) numtags=4;
           else if(taglen==2) numtags=16;
           else numtags=64;

           sfa_size = max_sfa_size/numtags *2; // *2 to avoid reallocating more memory as much as possible
           sfa = (int*)malloc(sfa_size*sizeof(int));
           if(sfa==NULL) { fprintf(stderr,"Failed to allocate memory to suffix array.\n"); exit(1); }
           s_index=0;
           for(t=0;t<=max_tid;t++){
             if(cuml[t+1]-cuml[t]>readlength){
               for(i=cuml[t];i<cuml[t+1]-readlength;i++){
                 flippos=flip(i);
                 if(strncmp(seq_array+i,seq_array+flippos,readlength)>0) {
                   if(strncmp(seq_array+flippos,tag,taglen)==0) {
                      if(is_noncanonical(seq_array+flippos,readlength)) continue; // include only ACGT
                      sfa[s_index]=flippos; s_index++;
                      if(s_index>=sfa_size) {
                        sfa_size+=max_sfa_size/numtags;
                        sfa = realloc(sfa,sfa_size*sizeof(int));
                        if(sfa==NULL) { fprintf(stderr,"Failed to reallocate memory to suffix array at tag=%s, tid=%d, s_index=%d\n",tag,t,s_index); exit(1); }
                      }
                   }
                 }
                 else {
                   if(strncmp(seq_array+i,tag,taglen)==0) {
                      if(is_noncanonical(seq_array+i,readlength)) continue; // include only ACGT
                      sfa[s_index]=i; s_index++;
                      if(s_index>=sfa_size) {
                        sfa_size+=max_sfa_size/numtags;
                        sfa = realloc(sfa,sfa_size*sizeof(int));
                        if(sfa==NULL) { fprintf(stderr,"Failed to reallocate memory to suffix array at tag=%s, tid=%d, s_index=%d\n",tag,t,s_index); exit(1); }
                      }
                   }
                 }
               }
             }
           }

           if(s_index>0){
             sfa_size=s_index;
             sfa = realloc(sfa,sfa_size*sizeof(int));
             if(sfa==NULL) { fprintf(stderr,"Failed to reallocate memory to suffix array at tag=%s in the end of the cycle\n",tag); exit(1); }
             local_sfa_start=0; local_sfa_end = sfa_size-1;
           }
           else { sfa_size=0; local_sfa_start=0; local_sfa_end = -1; free(sfa); sfa=NULL; }
}

void initialize_suffixarray_NS_PE_2(char* tag)
{
           int i,t,s_index,flippos;
           int taglen=strlen(tag);
           int numtags;

           if(taglen==1) numtags=4;
           else if(taglen==2) numtags=16;
           else numtags=64;
           local_sfa_start = local_sfa_end+1;  // local_sfa_end must be initialized to -1.

           if(sfa==NULL) { sfa = (int*)malloc(max_sfa_size*2*sizeof(int)); sfa_size=max_sfa_size*2; }
           if(sfa==NULL) { fprintf(stderr,"Failed to allocate memory to suffix array.\n"); exit(1); }
           s_index=local_sfa_start;
           for(t=0;t<=max_tid;t++){
             if(cuml[t+1]-cuml[t]>readlength){
               for(i=cuml[t];i<cuml[t+1]-readlength;i++){
                   //all fw substrings
                   if(strncmp(seq_array+i,tag,taglen)==0) {
                      if(getb(noncanonarr,i)==0) { // include only ACGT
                        sfa[s_index]=i; s_index++;
                      }
                   }
                   //all rc substrings
                   flippos=flip(i);
                   if(strncmp(seq_array+flippos,tag,taglen)==0) {
                      if(getb(noncanonarr,flippos)==0){ // include only ACGT
                        sfa[s_index]=flippos; s_index++;
                      }
                   }
               }
             }
           }
           local_sfa_end = s_index-1;

}


void initialize_suffixarray_SS_PE_2(char* tag)
{
           fprintf(stderr,"In initialize_suffixarray_SS_PE_2\n"); //DEBUGGING
           int i,t,s_index,flippos;
           int taglen=strlen(tag);
           int numtags;

           if(taglen==1) numtags=4;
           else if(taglen==2) numtags=16;
           else numtags=64;
           local_sfa_start = local_sfa_end+1;  // local_sfa_end must be initialized to -1.

           if(sfa==NULL) { sfa = (int*)malloc(max_sfa_size*sizeof(int)); sfa_size=max_sfa_size; }
           if(sfa==NULL) { fprintf(stderr,"Failed to allocate memory to suffix array.\n"); exit(1); }
           s_index=local_sfa_start;
           for(t=0;t<=max_tid;t++){
             if(cuml[t+1]-cuml[t]>readlength){
               for(i=cuml[t];i<cuml[t+1]-readlength;i++){
                   //all fw substrings
                   if(strncmp(seq_array+i,tag,taglen)==0) {
                      if(getb(noncanonarr,i)==0) { // include only ACGT
                        sfa[s_index]=i; s_index++;
                      }
                   }
               }
             }
           }
           local_sfa_end = s_index-1;
}

void *quick_sort_suffixarray_4(void *arg)
{
   int pivot_index, pivot_new_index;
   SFA_RANGE a = *((SFA_RANGE*)arg);
   SFA_RANGE a1,a2;
   int left = a.left;
   int right = a.right;
   int serial = 0;
   int i,thread2;

   if(right>left){
     pivot_index = left + (int)((right-left)/2);
     pivot_new_index = partition_suffix_array_3(left, right, pivot_index);
     a1.left=left; a1.right=pivot_new_index-1; 
     a2.left=pivot_new_index+1; a2.right=right; 

     //making decisions about threading
     if(MAX_Thread>1 && a2.right>a2.left){  // the latter condition is just to save nearly wasteful creation of threads.
       pthread_mutex_lock(&mutex1); 
       if(nThread<MAX_Thread-1) {
         nThread++; 
         serial = 0;
         for(i=0;i<MAX_Thread-1;i++) if(working_thread[i]==0) { thread2=i; working_thread[i]=1; break; }
       }
       else { serial=1; }
       pthread_mutex_unlock(&mutex1);
     } else { serial=1; }

    if( serial==1) {
       quick_sort_suffixarray_4(&a1);
       quick_sort_suffixarray_4(&a2);
    }
    else {
       pthread_create(&pth[thread2],NULL,quick_sort_suffixarray_4,&a2);
       quick_sort_suffixarray_4(&a1);
       pthread_join(pth[thread2],NULL); 
       pthread_mutex_lock(&mutex1);
       nThread--; working_thread[thread2]=0;  
       pthread_mutex_unlock(&mutex1);
    }
   }
}



void *quick_sort_suffixarray_4_nothread(void *arg)
{
   int pivot_index, pivot_new_index;
   SFA_RANGE a = *((SFA_RANGE*)arg);
   SFA_RANGE a1,a2;
   int left = a.left;
   int right = a.right;
   int i,thread2;

   if(right>left){
     pivot_index = left + (int)((right-left)/2);
     pivot_new_index = partition_suffix_array_3(left, right, pivot_index);
     a1.left=left; a1.right=pivot_new_index-1; 
     a2.left=pivot_new_index+1; a2.right=right; 

     quick_sort_suffixarray_4_nothread(&a1);
     quick_sort_suffixarray_4_nothread(&a2);
   }
}



int partition_suffix_array_3(int left, int right, int pivot_index)
{
    int tmp;
    int store_index,i;

    if(pivot_index>=sfa_size) { fprintf(stderr,"Error in partition_suffix_array_3 : pivot_index (%d) larger than sfa_size (%d)\n",pivot_index,sfa_size); exit(1); }
    if(right>=sfa_size) { fprintf(stderr,"Error in partition_suffix_array_3 : index right (%d) larger than sfa_size (%d)\n",right,sfa_size); exit(1); }
    tmp = sfa[pivot_index]; sfa[pivot_index] = sfa[right];  sfa[right]=tmp;
    store_index = left;
    for(i=left;i<right;i++){
      if(sfa[i]+readlength>fasta_size) { fprintf(stderr,"Error in partition_suffix_array_3 : sfa[i=%d] (%d) larger than fasta size (%d)\n",i,sfa[i],fasta_size); exit(1); }
      if(sfa[right]+readlength>fasta_size) { fprintf(stderr,"Error in partition_suffix_array_3 : sfa[right=%d] (%d) larger than fasta size (%d)\n",right,sfa[right],fasta_size); exit(1); }
      if(strncmp(seq_array+sfa[i], seq_array+sfa[right], readlength)<0) {
        tmp = sfa[i]; sfa[i]=sfa[store_index]; sfa[store_index]=tmp;
        store_index++;
      }
    }
    tmp = sfa[right]; sfa[right]=sfa[store_index]; sfa[store_index]=tmp;
    return store_index;
}


void quick_sort_suffixarray_by_mate_2(int** pSfa_s, int** pSfa_d, int left, int right)
{
   int pivot_index, pivot_new_index, left1,right1, left2,right2;

   if(right>left){
     pivot_index = left + (int)((right-left)/2);
     pivot_new_index = partition_suffix_array_by_mate_2(pSfa_s, pSfa_d, left, right, pivot_index);
     left1=left; right1=pivot_new_index-1; 
     left2=pivot_new_index+1; right2=right; 

     quick_sort_suffixarray_by_mate_2(pSfa_s,pSfa_d,left1,right1);
     quick_sort_suffixarray_by_mate_2(pSfa_s,pSfa_d,left2,right2);
   }
}

int partition_suffix_array_by_mate_2(int** pSfa_s, int** pSfa_d, int left, int right, int pivot_index)
{
    int tmp;
    int store_index,i;

    tmp = (*pSfa_s)[pivot_index]; (*pSfa_s)[pivot_index] = (*pSfa_s)[right];  (*pSfa_s)[right]=tmp;
    tmp = (*pSfa_d)[pivot_index]; (*pSfa_d)[pivot_index] = (*pSfa_d)[right];  (*pSfa_d)[right]=tmp;
    store_index = left;
    for(i=left;i<right;i++){
      if(strncmp(seq_array+(*pSfa_s)[i], seq_array+(*pSfa_s)[right],readlength)<0) { // mate1 is identical for all, so compare only mate 2, also assuming all the entries are 'canonical' (every entry has at least readlength from the mate2 start pos). 
        tmp = (*pSfa_s)[i]; (*pSfa_s)[i]=(*pSfa_s)[store_index]; (*pSfa_s)[store_index]=tmp;
        tmp = (*pSfa_d)[i]; (*pSfa_d)[i]=(*pSfa_d)[store_index]; (*pSfa_d)[store_index]=tmp;
        store_index++;
      }
    }
    tmp = (*pSfa_s)[right]; (*pSfa_s)[right]=(*pSfa_s)[store_index]; (*pSfa_s)[store_index]=tmp;
    tmp = (*pSfa_d)[right]; (*pSfa_d)[right]=(*pSfa_d)[store_index]; (*pSfa_d)[store_index]=tmp;
    return store_index;
}


int generate_seqtag(char*** tag, int taglen){
   int i,j,j2,k;
   int totaltags; 

   if(taglen==1){
     totaltags=4;
     (*tag)=(char**)malloc(totaltags*sizeof(char*));
     k=0;
     for(i=0;i<4;i++){
       (*tag)[k]=(char*)malloc(2*sizeof(char));
       (*tag)[k][0]=base[i]; (*tag)[k][1]='\0'; k++;
     }
   } else if(taglen==2){
     totaltags=16;
     (*tag)=(char**)malloc(totaltags*sizeof(char*));
     k=0;
     for(i=0;i<4;i++) for(j=0;j<4;j++) {
       (*tag)[k]=(char*)malloc(3*sizeof(char));
       (*tag)[k][0]=base[i]; (*tag)[k][1]=base[j]; (*tag)[k][2]='\0'; k++;
     }
   } else if(taglen==3){
     totaltags=64;
     (*tag)=(char**)malloc(totaltags*sizeof(char*));
     k=0;
     for(i=0;i<4;i++) for(j=0;j<4;j++) for(j2=0;j2<4;j2++) {
       (*tag)[k]=(char*)malloc(4*sizeof(char));
       (*tag)[k][0]=base[i]; (*tag)[k][1]=base[j]; (*tag)[k][2]=base[j2]; (*tag)[k][3]='\0'; k++;
     }
   }
   if(k!=totaltags) { fprintf(stderr,"error: tag number doesn't match between computed totaltags(%d) and actual totaltags(%d)\n",totaltags,k); exit(1); }
   return k; /* number of differen tags */
}

void delete_seqtag(char*** tag, int numtags)
{
   int i;
   for(i=0;i<numtags;i++) { free((*tag)[i]); (*tag)[i]=NULL; }
   free(*tag);
   (*tag)=NULL;
}



/*------- suffix array printing -------*/
void print_sfa ( void ) {
  int i;
  FILE *outfile_sfa;

  outfile_sfa = fopen(sfafile_name,"w");
  if(outfile_sfa==NULL){
      fprintf(stderr,"Can't write to output sfa file %s\n",sfafile_name);
      exit(1);
  }

  if(sfa_size==0 || sfa==NULL) {
    fprintf(stderr,"Suffix array is null. Not printing.\n");
    exit(1);
  }
  for(i=0;i<sfa_size;i++) {
    fprintf(outfile_sfa,"%d\t%d\n",i,sfa[i]);
  }
  fclose(outfile_sfa);
}


/*------- suffix array  marking -------*/

void mark_sfa_se(void){
     int i, sfa_m_size = bitsize(sfa_size);
     sfa_m = (char*)calloc(sfa_m_size,sizeof(char)); /* if sfa[i] and sfa[i-1] is identical strings, sfa_m[i]=0, otherwise, sfa_m[i]=1. bitstring. */
     if(sfa_m==NULL) { fprintf(stderr,"Failed to allocate memory to sfa_m\n"); exit(1); }
     for(i=1;i<sfa_size;i++)
        if(strdiff_se(seq_array+sfa[i],seq_array+sfa[i-1])!=0) putb(sfa_m,i,1);
}


void mark_sfa_pe(int d){
   int i,prev_i, sfa_m_pe_size = bitsize(sfa_size);
   sfa_m_pe = (char*)calloc(sfa_m_pe_size,sizeof(char)); /* if sfa[i] and sfa[prev_i] is identical strings, sfa_m[i]=0, otherwise, sfa_m[i]=1. bitstring. prev_i is the last position before i where sfa_f value is 1. */
   if(sfa_m_pe==NULL) { fprintf(stderr,"Failed to allocate memory to sfa_m_pe\n"); exit(1); }
   i=0;
   while(i<sfa_size && getb(sfa_f,i)==0) i++; prev_i=i;
   putb(sfa_m_pe,prev_i,1);
   do {i++;} while(i<sfa_size && getb(sfa_f,i)==0);
   while(i<sfa_size){
      if(strdiff_pe(seq_array+sfa[i],seq_array+sfa[prev_i],d)!=0) putb(sfa_m_pe,i,1);
      prev_i=i;
      do {i++;} while(i<sfa_size && getb(sfa_f,i)==0);
   }
}









/*------- rshbucket -------*/

void initialize_rshbucket(int max_t_size){
  int i;
  if(max_t_size==-1) rshbucket_max_t_size = INIT_RSHBUCKET_MAX_T_SIZE;
  else rshbucket_max_t_size = max_t_size;

  rshbucket = calloc(rshbucket_max_t_size-1, sizeof(node1**));
  if(rshbucket==NULL) { fprintf(stderr,"Failed to allocate memory to rshbucket\n"); exit(1); }

  rshbucket_single = (node1**)calloc(max_tid+1,sizeof(node1*));
  if(rshbucket_single==NULL) { fprintf(stderr,"Failed to allocate memory to rshbucket_single\n"); exit(1); }

  rsh_size=0;
  rsh_size_single=0;
}


/* build rsh, indexTable and tname_tree, library_strand_type. */
void construct_rsh_from_rshfile(char* input_rshfile_name){
   FILE *input_rshfile;
   node1 *lastp;
   char line[RSH_SINGLE_LINE_MAX];

   input_rshfile = fopen((const char*)input_rshfile_name,"r");
   /*Check for validity of the file.*/
   if(input_rshfile == 0)
   {
     printf("can't open input rsh file.\n");
     exit(1);
   }

   lastp=NULL;
   while (fgets(line, RSH_SINGLE_LINE_MAX, input_rshfile) != NULL) {
     if(line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0'; /*chomp*/
 
     if(line[0]=='#') parse_rsh_headerline(line);
     else if(line[0]=='@') parse_rsh_indexline(line);
     else if(line[0]!='c') { // column headings
       lastp=parse_rsh_mainline(line,lastp);
     }
   }
   fclose(input_rshfile);
   max_cid = max_tid + rsh_size;
   fprintf(stderr,"done reading rsh. rshsize=%d\n",rsh_size); //debugging

}


void parse_rsh_indexline (char* line){
       // the line looks like this : "@tid\ttname"
       char tmpstr[SINGLE_LINE_MAX],tname[SINGLE_LINE_MAX];
       int numtab,j,i,tid;

       /*The following parsing is a replacement of strtok which may cause memory leak. */
       numtab=0;j=0;
       for(i=0;i<=strlen(line);i++){
          if(line[i]=='\t'||line[i]=='\0') {
            numtab++;
            tmpstr[j]='\0';
            if(numtab==1) tid = atoi(tmpstr+1);
            else if(numtab==2) strcpy(tname,tmpstr);
            j=0;
          }
          else { tmpstr[j]=line[i]; j++; }
       }

       insert_key(tname,tid,tname_tree);  // tname -> tid hash
       IndexTable[tid] = (char*)malloc((strlen(tname)+1)*sizeof(char)); // tid -> tname array
       strcpy(IndexTable[tid],tname);

}


void parse_rsh_headerline (char* line){
       // the line looks like this : "#max_tid,max_t_size,minfrag,maxfrag,readlength"
       char tmpstr[SINGLE_LINE_MAX],tname[SINGLE_LINE_MAX];
       int numtab,j,i,max_t_size;

       /*The following parsing is a replacement of strtok which may cause memory leak. */
       numtab=0;j=0;
       for(i=0;i<=strlen(line);i++){
          if(line[i]==','||line[i]=='\0') {
            numtab++;
            tmpstr[j]='\0';
            if(numtab==1) max_tid = atoi(tmpstr+1);
            else if(numtab==2) max_t_size = atoi(tmpstr);
            else if(numtab==3) Min_Fraglength = atoi(tmpstr);
            else if(numtab==4) Max_Fraglength = atoi(tmpstr);
            else if(numtab==5) readlength = atoi(tmpstr);
            j=0;
          }
          else { tmpstr[j]=line[i]; j++; }
       }
 
       initialize_tname_indextable(); // This requires initializing max_tid, which is done within the current function.
       determine_fraglength_range();
       initialize_rshbucket(max_t_size);
}

node1* parse_rsh_mainline (char* line, node1* lastp){
       // the line looks like this : "cid\ttsize\ttid0\tother.tid\teuma.values"
       char tmpstr[RSH_SINGLE_LINE_MAX],othertidstr[RSH_SINGLE_LINE_MAX],eumastr[RSH_SINGLE_LINE_MAX];
       int numtab,j,i,t_size,t_ind,tid0,*tarray,*EUMA;
       node1 *q;

       /*The following parsing is a replacement of strtok which may cause memory leak. */
       numtab=0;j=0;
       for(i=0;i<=strlen(line);i++){
          if(line[i]=='\t'||line[i]=='\0') {
            numtab++;
            tmpstr[j]='\0';
            if(numtab==2) t_size = atoi(tmpstr);
            else if(numtab==3) tid0 = atoi(tmpstr);
            else if(numtab==4) strcpy(othertidstr,tmpstr);
            else if(numtab==5) strcpy(eumastr,tmpstr);
            j=0;
          }
          else { tmpstr[j]=line[i]; j++; }
       }


       //tarray
       if(t_size>1){
        tarray=(int*)malloc((t_size)*sizeof(int));
        tarray[0]=tid0;
        numtab=0;j=0;
        for(i=0;i<=strlen(othertidstr)-1;i++){  //-1 to remove the ',' at the end.
          if(othertidstr[i]==','||othertidstr[i]=='\0') {
            numtab++;
            tmpstr[j]='\0';
            tarray[numtab] = atoi(tmpstr);
            j=0;
          }
          else { tmpstr[j]=othertidstr[i]; j++; }
        }
       } else tarray=NULL; 


       //EUMAarray
       if(strlen(eumastr)>0){
        EUMA=(int*)malloc(nFraglen*sizeof(int));
        numtab=0;j=0;
        for(i=0;i<=strlen(eumastr)-1;i++){
          if(eumastr[i]==','||eumastr[i]=='\0') {
            numtab++;
            tmpstr[j]='\0';
            EUMA[numtab-1] = atoi(tmpstr);
            j=0;
          }
          else { tmpstr[j]=eumastr[i]; j++; }
        }
       }else EUMA=NULL;

       if(EUMA!=NULL) {
         q=newNode2(tarray,t_size,NULL,EUMA,NULL);
         if(t_size==1) rshbucket_single[tid0]=q;
         else{
           t_ind = t_size-2;
           if(rshbucket[t_ind]==NULL){
             rshbucket[t_ind]=calloc(max_tid+1,sizeof(node1*));
             if(rshbucket[t_ind]==NULL) { fprintf(stderr,"Failed to allocate memory to rshbucket[%d].\n",t_ind); exit(1); }
           }
           if( rshbucket[t_ind][tid0] ==NULL )  {
             rshbucket[t_ind][tid0]=q;
             lastp=q;
           }
           else {
             lastp->next = q;
             lastp=q;
           }
           rsh_size++;
         }
       }
       if(EUMA!=NULL) free(EUMA);
       if(tarray!=NULL) free(tarray);

       return(lastp);
}



char update_rshbucket_single(int tid, char type, int fraglen_ind, char* poscat){
    node1 *p;
    int i;

    if(type=='e'){
      if(rshbucket_single[tid]==NULL){
         rshbucket_single[tid]=newNode1(NULL,1,NULL,type,fraglen_ind,poscat); return 1;
      }
      else {
         p=rshbucket_single[tid];
         if(p->EUMA==NULL) p->EUMA=(int*)calloc(nFraglen,sizeof(int)); p->EUMA[fraglen_ind]++; 
         return 0;
      }
    }
    else // type=='r'
    {
      if(rshbucket_single[tid]==NULL){ return 0; }
      else {
         p=rshbucket_single[tid];
         p->ReadCount++; 
       }
       return 0;
    }
}



/* return 1 if new node is made, 0 otherwise */
char update_rshbucket(int t_size, int* tarray, char type, int fraglen_ind, char* poscat){
  int t_ind;
  node1 *p;
  short int cmp;
  int i;

  t_ind = t_size-2;

  if(type=='e'){
    if(t_size > rshbucket_max_t_size) {
       rshbucket = realloc(rshbucket,(t_size-1)*sizeof(node1**));
       if(rshbucket==NULL) { fprintf(stderr,"Failed to reallocate memory to rshbucket.\n"); exit(1); }
       for(i=rshbucket_max_t_size-1;i<=t_ind;i++) rshbucket[i]=NULL;
       rshbucket_max_t_size = t_size;
    }
    
    if(rshbucket[t_ind]==NULL){
       rshbucket[t_ind]=calloc(max_tid+1,sizeof(node1*));
       if(rshbucket[t_ind]==NULL) { fprintf(stderr,"Failed to allocate memory to rshbucket[%d].\n",t_ind); exit(1); }
    }
  
    if( rshbucket[t_ind][tarray[0]] ==NULL )  {
       rshbucket[t_ind][tarray[0]]=newNode1(tarray,t_size,NULL,type,fraglen_ind,poscat); return 1;
    }
    else {
       p=rshbucket[t_ind][tarray[0]];
       cmp=cmptarr(tarray,p->tarr,t_size);
       if(cmp<0) {
          rshbucket[t_ind][tarray[0]] = newNode1(tarray,t_size,p,type,fraglen_ind,poscat);
          return 1;
       }    
       else if(cmp==0){
          if(p->EUMA==NULL) p->EUMA=(int*)calloc(nFraglen,sizeof(int)); p->EUMA[fraglen_ind]++; 
          return 0;
       }    
       
       while(1){
         if(p->next!=NULL){
             cmp=cmptarr(tarray,p->next->tarr,t_size);
             if(cmp<0) {
               p->next = newNode1(tarray,t_size,p->next,type,fraglen_ind,poscat);
               return 1;
             }    
             else if(cmp==0){
               if(p->next->EUMA==NULL) p->next->EUMA=(int*)calloc(nFraglen,sizeof(int)); p->next->EUMA[fraglen_ind]++; 
               return 0;
             }    
             else  p= p->next; 
         }
         else { 
             p->next = newNode1(tarray,t_size,NULL,type,fraglen_ind,poscat); return 1; 
         }
       }
    }
  }
  else //if(type=='r')
  {
    if(t_size > rshbucket_max_t_size) { return 0; }
    if(rshbucket[t_ind]==NULL){ return 0; }
    if( rshbucket[t_ind][tarray[0]] ==NULL )  { return 0; }
    else {
       p=rshbucket[t_ind][tarray[0]];
       cmp=cmptarr(tarray,p->tarr,t_size);
       if(cmp<0) { return 0; }    
       else if(cmp==0){
          p->ReadCount++; //type == 'r'
          return 0;
       }    
       
       while(1){
         if(p->next!=NULL){
             cmp=cmptarr(tarray,p->next->tarr,t_size);
             if(cmp<0) { return 0; }    
             else if(cmp==0){
               p->next->ReadCount++; 
               return 0;
             }    
             else  p= p->next; 
         }
         else { return 0; } 
       }
    }
  }
}




node1* newNode1 (int *t_array, int t_size, node1 *next, char type, int fraglen_ind, char* poscat){
   int i;
   node1 *p;

   if(t_array==NULL) { // for rshbucket_single
      p = malloc(sizeof(node1));
      p->tarr[0]=0; // waste of memory. (This space is not used)
   }
   else {
      p =malloc(sizeof(node1) + sizeof(int)*(t_size-2));
      for(i=1;i<t_size;i++) p->tarr[i-1]=t_array[i]; // tarray is copied.
   }

   p->next=next;
   p->EUMA=NULL;
   if(type=='e')  {
     p->EUMA=(int*)calloc( nFraglen, sizeof(int) );
     p->EUMA[fraglen_ind]=1;
     p->ReadCount=0;
   }
   else {
     p->ReadCount=1; // type 'r'
   }
   return p;
}


node1* newNode2 (int *t_array, int t_size, node1 *next, int *EUMA, char* poscat){
   int i;
   node1 *p;

   if(t_array==NULL) { // for rshbucket_single
      p = malloc(sizeof(node1));
      p->tarr[0]=0; // waste of memory. (This space is not used)
   }
   else {
      p =malloc(sizeof(node1) + sizeof(int)*(t_size-2));
      for(i=1;i<t_size;i++) p->tarr[i-1]=t_array[i]; // tarray is copied.
   }

   p->next=next;

   p->EUMA=(int*)calloc( nFraglen, sizeof(int) );
   for(i=0;i<nFraglen;i++) p->EUMA[i]=EUMA[i];
   return p;
}

short int cmptarr (int* t_array1, int* t_array2,int t_size){
   int i;
   for(i=1;i<t_size;i++) {
     if(t_array1[i]<t_array2[i-1]) return -1;
     else if(t_array1[i]>t_array2[i-1]) return 1;
   }
   return 0;
}


void delete_rshbucket(void){
   int i,j;
   node1 *p,*q;
   for(i=0;i<=rshbucket_max_t_size-2;i++)
       if(rshbucket[i]!=NULL){
          for(j=0;j<=max_tid;j++)
            if(rshbucket[i][j]!=NULL){
               q=rshbucket[i][j];
               do {
                 p=q->next;
                 free(q->EUMA);
                 //if(q->poscat!=NULL) free(q->poscat);
                 free(q);
                 q=p;
               }while(q!=NULL);
            }
          free(rshbucket[i]);
          rshbucket[i]=NULL;
       }
   free(rshbucket);
   rshbucket=NULL;

   /* deleting rshbucket_single */
   for(j=0;j<=max_tid;j++)
     if(rshbucket_single[j]!=NULL){
        q=rshbucket_single[j];
        do {
          p=q->next;
          free(q->EUMA);
          //if(q->poscat!=NULL) free(q->poscat);
          free(q);
          q=p;
        }while(q!=NULL);
     }
   free(rshbucket_single);
   rshbucket_single=NULL;
}


void clear_readcounts_in_rshbucket(void){
   int i,j;
   node1 *p,*q;
   for(i=0;i<=rshbucket_max_t_size-2;i++)
       if(rshbucket[i]!=NULL){
          for(j=0;j<=max_tid;j++)
            if(rshbucket[i][j]!=NULL){
               q=rshbucket[i][j];
               do {
                 q->ReadCount=0;
                 p=q->next;
                 q=p;
               }while(q!=NULL);
            }
       }

   /* deleting rshbucket_single */
   for(j=0;j<=max_tid;j++)
     if(rshbucket_single[j]!=NULL){
        q=rshbucket_single[j];
        do {
          q->ReadCount=0;
          p=q->next;
          q=p;
        }while(q!=NULL);
     }
}



/* version 2 uses MAX_REPEAT cut :i.e. if the same sequence appears MAX_REPEAT or more times, skip it (eg. polyA) */
/* constructing rshbucket and determining rsh size */
void construct_rshbucket_2(int fraglength) //for SE, readlength==fraglength.
{
     int start_i,i,j,j2;
     int NMcount;
     int ti,ti2,t_size,nopush;
     int* t;
     int sfj,sfstart_i;
     int fraglen_ind = fraglength-Fraglengths.min;

     if(sfa_size>0){
        start_i=sfa_size-1; i=sfa_size-1; 
      
        for(i=sfa_size-2;i>=-1;i--){   /* i== -1 means it is now out of the loop. */
          if(verbose_flag>1 && sfa_size>10 && i%(sfa_size/10)==0) fprintf(stdout,"\r%3d%% done...",100-i/(sfa_size/10)*10);  /* progress log */
   
          /* out of loop or different sequence */
          if(i==-1 || strdiff_se(seq_array+sfa[start_i], seq_array+sfa[i]) != 0 ) { 
              sfstart_i=sf_i(sfa[start_i]);
              if(i == start_i-1) {    /* isoform-wise informative read */
                (*update_rshbucket_single_PTR)(sfstart_i,'e',fraglen_ind,NULL);
              }
              else {

                 if(start_i-i<MAX_REPEAT){
                   t=(int*)malloc((start_i-i)*sizeof(int));
                   if(t==NULL) { fprintf(stderr,"Failed to allocate memory to t in construct_rshbucket_2.\n"); exit(1); }
                   t[0]=sfstart_i;
                   t_size=1;

                   /* constructing t-array in a sorted, non-redundant manner */
                   for(j=start_i-1;j>i;j--){
                      nopush=0; /* indicator of whether to push the j'th isoform into the t_array */
                      sfj=sf_i(sfa[j]);
                      for(ti=0; ti<t_size; ti++){
                         if(sfj<=t[ti]) {  /* insert into array */
                             for(ti2=t_size-1; ti2>=ti; ti2--) { t[ti2+1] = t[ti2]; }
                             t[ti]=sfj;
                             t_size++;
                             nopush=1;
                             break;
                         }
                      }
                      if(nopush==0) {  /* push to the array */
                         t[t_size]=sfj; 
                         t_size++;
                      }
                   }
                   if((*update_rshbucket_PTR)(t_size,t,'e',fraglen_ind,NULL)) rsh_size++;  //rsh_size will be used to compute max_cid.
                   free(t); t=NULL;
                 }
              }
   
              if(i>=0) {
                start_i=i;
                //sfa=realloc(sfa,(start_i+1)*sizeof(int));  /* delete the used part of the suffix array as we go up */
              }
          }
   
        }
     }
     max_cid = max_tid + rsh_size;
}

/* version 2 uses MAX_REPEAT cut :i.e. if the same sequence appears MAX_REPEAT or more times, skip it (eg. polyA) */
/* constructing rshbucket and determining rsh size */
void construct_rshbucket_2_posbias(int fraglength) //for SE, readlength==fraglength.
{
     int start_i,i,j,j2;
     int NMcount;
     int ti,ti2,t_size,nopush;
     int* t;
     int sfj,sfstart_i;
     int fraglen_ind = fraglength-Fraglengths.min;
     char* poscat;

     if(sfa_size>0){
        start_i=sfa_size-1; i=sfa_size-1; 
      
        for(i=sfa_size-2;i>=-1;i--){   /* i== -1 means it is now out of the loop. */
          if(verbose_flag>1 && sfa_size>10 && i%(sfa_size/10)==0) fprintf(stdout,"\r%3d%% done...",100-i/(sfa_size/10)*10);  /* progress log */
   
          /* out of loop or different sequence */
          if(i==-1 || strdiff_se(seq_array+sfa[start_i], seq_array+sfa[i]) != 0 ) { 
              sfstart_i=sf_i(sfa[start_i]);
              if(i == start_i-1) {    /* isoform-wise informative read */
                poscat = (char*)malloc(sizeof(char));
                poscat[0] = determine_poscategory(sfstart_i,sf_p(sfa[start_i],sfstart_i),fraglen_ind+Fraglengths.min);
                if((*update_rshbucket_single_PTR)(sfstart_i,'e',fraglen_ind,poscat)) {
                  poscat=NULL; // poscat is transferred to rshbucket_single and will be freed later as rshbucket_single is deleted.
                } else { free(poscat); poscat=NULL; }
              }
              else {

                 if(start_i-i<MAX_REPEAT){
                   t=(int*)malloc((start_i-i)*sizeof(int));
                   if(t==NULL) { fprintf(stderr,"Failed to allocate memory to t in construct_rshbucket_2_posbias.\n"); exit(1); }
                   t[0]=sfstart_i;
                   t_size=1;

                   poscat = (char*)malloc((start_i-i)*sizeof(char));
                   if(poscat==NULL) { fprintf(stderr,"Failed to allocate memory to poscat in construct_rshbucket_2_posbias.\n"); exit(1); }
                   poscat[0] = determine_poscategory(sfstart_i,sf_p(sfa[start_i],sfstart_i),fraglen_ind+Fraglengths.min);

                   /* constructing t-array in a sorted, non-redundant manner */
                   for(j=start_i-1;j>i;j--){
                      nopush=0; /* indicator of whether to push the j'th isoform into the t_array */
                      sfj=sf_i(sfa[j]);
                      for(ti=0; ti<t_size; ti++){
                         if(sfj<=t[ti]) {  /* insert into array */
                             for(ti2=t_size-1; ti2>=ti; ti2--) { t[ti2+1] = t[ti2]; poscat[ti2+1] = poscat[ti2]; }
                             t[ti]=sfj;
                             poscat[ti] = determine_poscategory(sfj,sf_p(sfa[j],sfj),fraglen_ind+Fraglengths.min);
                             t_size++;
                             nopush=1;
                             break;
                         }
                      }
                      if(nopush==0) {  /* push to the array */
                         t[t_size]=sfj; 
                         poscat[t_size] = determine_poscategory(sfj,sf_p(sfa[j],sfj),fraglen_ind+Fraglengths.min);
                         t_size++;
                      }
                   }
                   if((*update_rshbucket_PTR)(t_size,t,'e',fraglen_ind,poscat)) {
                     rsh_size++;  //rsh_size will be used to compute max_cid.
                     poscat=NULL; // poscat is transferred to rshbucket_single and will be freed later as rshbucket_single is deleted.
                   } else { free(poscat); poscat=NULL; }
                   free(t); t=NULL;
                 }
              }
   
              if(i>=0) {
                start_i=i;
                //sfa=realloc(sfa,(start_i+1)*sizeof(int));  /* delete the used part of the suffix array as we go up */
              }
          }
   
        }
     }
     max_cid = max_tid + rsh_size;
}


/* constructing rshbucket and determining rsh size */
void construct_rshbucket_PE_3(int** pSfa_s, int** pSfa_d, int size) 
{
     int start_i,i,j,j2;
     int ti,ti2,t_size,nopush;
     int* t;
     int sfj,sfstart_i;
     int fraglen_ind,multi_d; 

     if(size>0){
        start_i=0;
        i=start_i+1;

        while(i<=size){ //i==size means it is now out of the loop.
   
          /* out of loop or different sequence */
          if(i==size || strdiff_se(seq_array+(*pSfa_s)[start_i],seq_array+(*pSfa_s)[i])) { 
              sfstart_i=sf_i((*pSfa_s)[start_i]);
              if(i-start_i==1) {    /* isoform-wise informative read */
                fraglen_ind = (*pSfa_d)[start_i]+readlength-Fraglengths.min;
                (*update_rshbucket_single_PTR)(sfstart_i,'e',fraglen_ind,NULL);
              }
              else {

                 multi_d=0;
                 for(j=start_i+1;j<i;j++) if((*pSfa_d)[j]!=(*pSfa_d)[start_i]) { multi_d=1; break;}  // filtering out substrings mapped with multiple d values.

                 if(multi_d==0 && i-start_i<MAX_REPEAT){
                   pthread_mutex_lock(&mutex1);
                   t=(int*)malloc((i-start_i)*sizeof(int));
                   pthread_mutex_unlock(&mutex1);
                   if(t==NULL) { fprintf(stderr,"Failed to allocate memory to t in construct_rshbucket_PE_3.\n"); exit(1); }
                   t[0]=sfstart_i;
                   t_size=1;

                   /* constructing t-array in a sorted, non-redundant manner */
                   for(j=start_i+1;j<i;j++){
                      nopush=0; /* indicator of whether to push the j'th isoform into the t_array */
                      sfj=sf_i((*pSfa_s)[j]);
                      for(ti=0; ti<t_size; ti++){
                         if(sfj<=t[ti]) {  /* insert into array */
                             for(ti2=t_size-1; ti2>=ti; ti2--) { t[ti2+1] = t[ti2]; }
                             t[ti]=sfj;
                             t_size++;
                             nopush=1;
                             break;
                         }
                      }
                      if(nopush==0) {  /* push to the array */
                         t[t_size]=sfj; 
                         t_size++;
                      }
                   }

                   fraglen_ind = (*pSfa_d)[start_i]+readlength-Fraglengths.min;
                   if((*update_rshbucket_PTR)(t_size,t,'e',fraglen_ind,NULL)) {
                      pthread_mutex_lock(&mutex1);
                      rsh_size++; //rsh_size will be used to compute max_cid.
                      pthread_mutex_unlock(&mutex1);
                   }
                   free(t); t=NULL;
                 }
              }
   
              if(i<size) {
                start_i=i;
                //sfa=realloc(sfa,(start_i+1)*sizeof(int));  /* delete the used part of the suffix array as we go up */
              }
          }
          i++;
        }
     }
    // max_cid=max_tid+rsh_size;
}

/* constructing rshbucket and determining rsh size */
void construct_rshbucket_PE_3_posbias(int** pSfa_s, int** pSfa_d, int size) 
{
     int start_i,i,j,j2;
     int ti,ti2,t_size,nopush;
     int* t;
     int sfj,sfstart_i;
     int fraglen_ind,multi_d; 
     char* poscat;

     if(size>0){
        start_i=0;
        i=start_i+1;

        while(i<=size){ //i==size means it is now out of the loop.
   
          /* out of loop or different sequence */
          if(i==size || strdiff_se(seq_array+(*pSfa_s)[start_i],seq_array+(*pSfa_s)[i])) { 
              fprintf(stderr,"In rsh: i=%d out of size=%d\n",i,size); //DEBUGGING
              sfstart_i=sf_i((*pSfa_s)[start_i]);
              if(i-start_i==1) {    /* isoform-wise informative read */
                fraglen_ind = (*pSfa_d)[start_i]+readlength-Fraglengths.min;
                poscat = (char*)malloc(sizeof(char));
                poscat[0] = determine_poscategory(sfstart_i,sf_p(sfa[start_i],sfstart_i),fraglen_ind+Fraglengths.min);
                if((*update_rshbucket_single_PTR)(sfstart_i,'e',fraglen_ind,poscat)){
                  poscat=NULL; // poscat is transferred to rshbucket_single and will be freed later as rshbucket_single is deleted.
                } else { free(poscat); poscat=NULL; }
              }
              else {

                 multi_d=0;
                 for(j=start_i+1;j<i;j++) if((*pSfa_d)[j]!=(*pSfa_d)[start_i]) { multi_d=1; break;}  // filtering out substrings mapped with multiple d values.

                 if(multi_d==0 && i-start_i<MAX_REPEAT){
                   fraglen_ind = (*pSfa_d)[start_i]+readlength-Fraglengths.min;

                   pthread_mutex_lock(&mutex1);
                   t=(int*)malloc((i-start_i)*sizeof(int));
                   pthread_mutex_unlock(&mutex1);
                   if(t==NULL) { fprintf(stderr,"Failed to allocate memory to t in construct_rshbucket_PE_3_posbias.\n"); exit(1); }
                   t[0]=sfstart_i;
                   t_size=1;

                   pthread_mutex_lock(&mutex1);
                   poscat = (char*)malloc((i-start_i)*sizeof(char));
                   pthread_mutex_unlock(&mutex1);
                   if(poscat==NULL) { fprintf(stderr,"Failed to allocate memory to poscat in construct_rshbucket_PE_3_posbias.\n"); exit(1); }
                   poscat[0] = determine_poscategory(sfstart_i,sf_p(sfa[start_i],sfstart_i),fraglen_ind+Fraglengths.min);

                   /* constructing t-array in a sorted, non-redundant manner */
                   for(j=start_i+1;j<i;j++){
                      nopush=0; /* indicator of whether to push the j'th isoform into the t_array */
                      sfj=sf_i((*pSfa_s)[j]);
                      for(ti=0; ti<t_size; ti++){
                         if(sfj<=t[ti]) {  /* insert into array */
                             for(ti2=t_size-1; ti2>=ti; ti2--) { t[ti2+1] = t[ti2]; poscat[ti2+1] = poscat[ti2];}
                             t[ti]=sfj;
                             poscat[ti] = determine_poscategory(sfj,sf_p(sfa[j],sfj),fraglen_ind+Fraglengths.min);
                             t_size++;
                             nopush=1;
                             break;
                         }
                      }
                      if(nopush==0) {  /* push to the array */
                         t[t_size]=sfj; 
                         poscat[t_size] = determine_poscategory(sfj,sf_p(sfa[j],sfj),fraglen_ind+Fraglengths.min);
                         t_size++;
                      }
                   }
                   fprintf(stderr,"In rsh: after all poscat\n"); //DEBUGGING

                   if((*update_rshbucket_PTR)(t_size,t,'e',fraglen_ind,poscat)) {
                        pthread_mutex_lock(&mutex1);
                        rsh_size++; //rsh_size will be used to compute max_cid.
                        pthread_mutex_unlock(&mutex1);
                        poscat=NULL; // poscat is transferred to rshbucket_single and will be freed later as rshbucket_single is deleted.
                   } else { free(poscat); poscat=NULL; }
                   free(t); t=NULL;
                 }
              }
   
              if(i<size) {
                start_i=i;
                //sfa=realloc(sfa,(start_i+1)*sizeof(int));  /* delete the used part of the suffix array as we go up */
              }
          }
          i++;
        }
     }
    // max_cid=max_tid+rsh_size;
}



/* print rsh */
void print_rsh (char* rshfile_name){
   int cid; 
   int i,j,k,t_size;
   node1* p;
   FILE *outfile_rsh;

   outfile_rsh = fopen(rshfile_name,"w");
   if(outfile_rsh==NULL){
      fprintf(stderr,"Can't write to output rsh file %s\n",rshfile_name);
      exit(1);
   }


   //first header
   if(pe==1) fprintf(outfile_rsh,"#%d,%d,%d,%d,%d\n",max_tid,rshbucket_max_t_size,Fraglengths.min,Fraglengths.max,readlength);
   else fprintf(outfile_rsh,"#%d,%d,%d,%d,%d\n",max_tid,rshbucket_max_t_size,Fraglengths.min,Fraglengths.max,-1);

   //indexlines
   for(j=0;j<=max_tid;j++){
     fprintf(outfile_rsh,"@%d\t%s\n",j,IndexTable[j]);
   }

   //column headings
   fprintf(outfile_rsh,"cid\tno.tids\tfirst.tid\tother.tids\tsegment.length\n"); 

   //single-t
   cid=0;
   for(j=0;j<=max_tid;j++){
       p=rshbucket_single[j];
       if(p==NULL) { 
             fprintf(outfile_rsh,"%d\t%d\t%d\t\t\t\n",cid,1,j);
             cid++;
       }
       else while(p!=NULL){
             fprintf(outfile_rsh,"%d\t%d\t%d\t\t",cid,1,j);
             for(k=0;k<nFraglen;k++) fprintf(outfile_rsh,"%d,",p->EUMA[k]);
             fprintf(outfile_rsh,"\n");
             cid++;
             p=p->next;
       }
   }

   //multi-t
   for(t_size=2;t_size<=rshbucket_max_t_size;t_size++)  /* i = t_size */
     if(rshbucket[t_size-2]!=NULL)
        for(j=0;j<=max_tid;j++){
            p=rshbucket[t_size-2][j];
            while(p!=NULL){
                  fprintf(outfile_rsh,"%d\t%d\t%d\t",cid,t_size,j);
                  for(k=1;k<t_size;k++) fprintf(outfile_rsh,"%d,",p->tarr[k-1]); 
                  fprintf(outfile_rsh,"\t");
                  for(k=0;k<nFraglen;k++) fprintf(outfile_rsh,"%d,",p->EUMA[k]);
                  fprintf(outfile_rsh,"\n");
                  cid++;
                  p=p->next;
            }
        }

    fclose(outfile_rsh);
}



/* scan rshbucket and construct EUMA and CT array, SE & PE */
void scan_rshbucket (void){
   int cid; 
   int i,j,k,t_size;
   node1* p;

   if(verbose_flag>0) fprintf(stdout,"max_tid=%d, rshsize=%d, rshsize_single=%d, max_cid=%d\n",max_tid,rsh_size,rsh_size_single,max_cid);

   adjEUMA = (double*)calloc(max_cid+1,sizeof(double));
   if(adjEUMA==NULL) { fprintf(stderr,"Failed to allocate memory to adjEUMA.\n"); exit(1); }
   ReadCount = (int*)calloc(max_cid+1,sizeof(int));
   if(ReadCount==NULL) { fprintf(stderr,"Failed to allocate memory to ReadCount.\n"); exit(1); }
   CT=(inta2*)malloc((max_cid+1)*sizeof(inta2));
   if(CT==NULL) { fprintf(stderr,"Failed to allocate memory to CT.\n"); exit(1); }

   //single-t
   cid=0;
   for(j=0;j<=max_tid;j++){
       p=rshbucket_single[j];
       if(p==NULL) {
             adjEUMA[cid] = 0;
             ReadCount[cid] = 0;
             CT[cid].size = 1;
             CT[cid].array = (int*)malloc(sizeof(int));
             CT[cid].array[0]=j;
             CT[cid].array2=NULL;
             cid++;
       }
       else while(p!=NULL){
             adjEUMA[cid] = compute_adjEUMA(p->EUMA);
             ReadCount[cid] = p->ReadCount;
             CT[cid].size = 1;
             CT[cid].array = (int*)malloc(sizeof(int));
             CT[cid].array[0]=j;
             CT[cid].array2=NULL;
             cid++;
             p=p->next;
       }
   }


   //multi-t
   for(t_size=2;t_size<=rshbucket_max_t_size;t_size++)  /* i = t_size */
     if(rshbucket[t_size-2]!=NULL)
        for(j=0;j<=max_tid;j++){
            p=rshbucket[t_size-2][j];
            while(p!=NULL){
                  adjEUMA[cid] = compute_adjEUMA(p->EUMA);
                  ReadCount[cid] = p->ReadCount;
                  CT[cid].size = t_size;
                  CT[cid].array = (int*)malloc(t_size*sizeof(int));
                  CT[cid].array[0]=j;
                  for(k=1;k<t_size;k++) CT[cid].array[k]=p->tarr[k-1];
                  CT[cid].array2=NULL;
                  cid++;
                  p=p->next;
            }
        }
}




/*------- sets of segments, TC/CT/SC/ST -------*/

/* version 2 uses max_cid instead of rsh_size */
/* transcripts with no i-informative regions are included in both TC and CT arrays. */
void build_TC_from_CT_2(void){
  unsigned int cid;
  int i,tid;

  TC=(TCE*)malloc((max_tid+1)*sizeof(TCE));
  if(TC==NULL) { fprintf(stderr,"Failed to allocate memory to TC.\n"); exit(1); }
  for(i=0;i<=max_tid;i++) { TC[i].list=NULL; TC[i].last=NULL; }

  for(cid=0;cid<=max_cid;cid++){
    for(i=0;i<CT[cid].size;i++){
        tid=CT[cid].array[i];
        if(TC[tid].last == NULL) {  //initial
           TC[tid].list=(LIST*)malloc(sizeof(LIST));
           if(TC[tid].list==NULL) { fprintf(stderr,"Failed to allocate memory to TC[tid=%d].list.\n",tid); exit(1); }
           TC[tid].last = TC[tid].list;
        }
        else { 
           TC[tid].last->next=(LIST*)malloc(sizeof(LIST));
           if(TC[tid].list->next==NULL) { fprintf(stderr,"Failed to allocate memory to TC[tid=%d].list->next.\n",tid); exit(1); }
           TC[tid].last = TC[tid].last->next;
        }
        TC[tid].last->val=cid;
        TC[tid].last->next=NULL;
    }

  }
}


/* in v0.9.7c, EUMAcut is back. */
/* Now in v0.7.2, MIN_EUMA is no longer used. */
/* version 2 returns 0 if the cid's EUMA is less than MIN_EUMA and the cid is not a single-isoform cid. IN this case, it is not included in the joining. */
/* return 0 if the cid is already taken */
int propagate_2(unsigned int cid, unsigned int T){
  int tid,i;
  unsigned int cid2;
  LIST* j;
  int wasnew=0;

  if(CS[cid]==-1){  //CS was initialized to -1 for all values.
    wasnew=1;
    if(CT[cid].size>1 && adjEUMA[cid]<EUMAcut) return 0;  /* EUMA cut */
    CS[cid]=T;
    for(i=0;i<CT[cid].size;i++){
      tid=CT[cid].array[i];
      if(TS[tid]==-1){  //TS was also initialized to -1 for all values.
        TS[tid]=T;
        j=TC[tid].list;
        ntid_per_sid++;
        while(j!=NULL){
          cid2=j->val;
          propagate_2(cid2,T);
          j=j->next;
        }
      }
    }
  }
  return wasnew;
}


void print_aEUMA_3(char* outfilename_aEUMA, int max_sid){
  unsigned int i;
  int j,tid,round;
  double FPKMmean,expected_readcount;
  FILE *outfile_aEUMA;

  outfile_aEUMA = fopen(outfilename_aEUMA,"w");
  if(outfile_aEUMA==NULL){
      fprintf(stderr,"Can't write to output aEUMA file %s\n",outfilename_aEUMA);
      exit(1);
  }

  fprintf(outfile_aEUMA,"segment_id\tsequence_sharing_set_id\ttranscript_id\ttranscript_names\teff.length\tReadcount\texpected_Readcount\n");
  for(i=0;i<=max_cid;i++) {
     fprintf(outfile_aEUMA,"c%d\ts%d\t",i,CS[i]);
     for(j=0;j<CT[i].size;j++) {
       if(j>0) fprintf(outfile_aEUMA,",");
       fprintf(outfile_aEUMA,"t%d",CT[i].array[j]);
     }
     fprintf(outfile_aEUMA,"\t");
     for(j=0;j<CT[i].size;j++) {
       if(j>0) fprintf(outfile_aEUMA,"+");
       fprintf(outfile_aEUMA,"%s",IndexTable[CT[i].array[j]]);
     }
     fprintf(outfile_aEUMA,"\t%lf",adjEUMA[i]);

     // read count and expected read count
     expected_readcount=0;
     for(j=0;j<CT[i].size;j++) {
       tid=CT[i].array[j];
       FPKMmean=0;
       for(round=0;round<NUM_ROUND;round++) FPKMmean += FPKMfinal[round][tid];
       FPKMmean/=NUM_ROUND;
       expected_readcount += FPKMmean*(adjEUMA[i]/1E3)*((double)TotalReadCount/1E6);
     }
     fprintf(outfile_aEUMA,"\t%d\t%f\n",ReadCount[i],expected_readcount);
  }
  fclose(outfile_aEUMA);
}


void generate_SC_ST(void){
  unsigned int i;
  int j,cid,sid;
  int tid_exist;

  SC = (inta*)malloc((max_sid+1)*sizeof(inta));
  ST = (inta*)malloc((max_sid+1)*sizeof(inta));
  if(SC==NULL) { fprintf(stderr,"Failed to allocate memory to SC.\n"); exit(1); }
  if(ST==NULL) { fprintf(stderr,"Failed to allocate memory to ST.\n"); exit(1); }

  for(sid=0;sid<=max_sid;sid++) { SC[sid].array=NULL; SC[sid].size=0; }
  for(sid=0;sid<=max_sid;sid++) { ST[sid].array=NULL; ST[sid].size=0; }

  for(cid=0;cid<=max_cid;cid++) { 
     sid=CS[cid];

     /* update SC */
     if(sid==-1) continue;

     if(SC[sid].array==NULL) {
        SC[sid].array = (int*)malloc(sizeof(int));
        SC[sid].array[0]=cid;
        SC[sid].size = 1;
     }
     else {
        SC[sid].array = realloc(SC[sid].array,(SC[sid].size+1)*sizeof(int));
        SC[sid].array[SC[sid].size]=cid;
        SC[sid].size++;
     }

     /* update ST */
     if(ST[sid].array==NULL) {
        ST[sid].array = (int*)malloc(CT[cid].size*sizeof(int));
        ST[sid].size=0;
        for(i=0;i<CT[cid].size;i++) {
          tid_exist=0;
          for(j=0;j<i;j++) if(ST[sid].array[j]==CT[cid].array[i] ) { tid_exist=1; break; }  // nonredundantly, because some cids contain multiple occurrences of the same tid's due to internal repeats.
          if(tid_exist==0) { ST[sid].array[j]=CT[cid].array[i]; ST[sid].size++; }
        }
        if(CT[cid].size!=ST[sid].size) ST[sid].array=realloc(ST[sid].array,ST[sid].size*sizeof(int));
     }
     else {
        for(j=0;j<CT[cid].size;j++) {
          tid_exist = 0;
          for(i=0;i<ST[sid].size;i++) if(ST[sid].array[i]==CT[cid].array[j]) { tid_exist = 1; break; }
          if(!tid_exist) {  /* no need to sort, but redundant tid's cannot be allowed. */
            ST[sid].array = realloc(ST[sid].array,(ST[sid].size+1)*sizeof(int));
            ST[sid].array[ST[sid].size]=CT[cid].array[j];
            ST[sid].size++;
          }
        }
     }
  }
}


void initialize_CS_TS(void)
{
      int i;

      CS=(int*)malloc((max_cid+1)*sizeof(int));  /* combination-set array (storing set ID) */
      TS=(int*)malloc((max_tid+1)*sizeof(int));  /* transcript-set array (storing set ID) */
      if(CS==NULL) { fprintf(stderr,"Failed to allocate memory to CS\n"); exit(1); }
      if(TS==NULL) { fprintf(stderr,"Failed to allocate memory to TS\n"); exit(1); }

      for(i=0;i<=max_cid;i++) CS[i]=-1;
      for(i=0;i<=max_tid;i++) TS[i]=-1;
}

void delete_CS_TS(void)
{
  //TS
  free(CS); CS=NULL;
  free(TS); TS=NULL; 
}

void delete_CT_SC_ST_EUMA(void)
{
  int cid,sid,tid;
  LIST *p,*q;

  //CT
  for(cid=0;cid<=max_cid;cid++){
   if(CT[cid].array!=NULL) { free(CT[cid].array); CT[cid].array=NULL; }
   if(CT[cid].array2!=NULL) { free(CT[cid].array2); CT[cid].array2=NULL; }
  }
  free(CT); CT=NULL;

  //SC
  for(sid=0;sid<=max_sid;sid++){
    if(SC[sid].array!=NULL) { free(SC[sid].array); SC[sid].array=NULL; }
  }
  free(SC); SC=NULL;

  //ST
  for(sid=0;sid<=max_sid;sid++){
    if(ST[sid].array!=NULL) { free(ST[sid].array); ST[sid].array=NULL; }
  }
  free(ST); ST=NULL;

  //adjEUMA
  free(adjEUMA); adjEUMA=NULL;

  //TC
  for(tid=0;tid<=max_tid;tid++){
    p=TC[tid].list;
    while(p!=NULL){
      q=p->next;
      free(p);
      p=q;
    }
  }
  free(TC); TC=NULL;
}

void print_SC(void) {
  int sid,i;
  printf("==SC==\n");
  for(sid=0;sid<=max_sid;sid++){
   if(SC[sid].array!=NULL) {
     for(i=0;i<SC[sid].size;i++) printf("sid=%d, SC[sid].array[%d]=%d\n",sid,i,SC[sid].array[i]);
   } else { printf("sid=%d, NULL\n",sid); }
  }
}



void print_CT(void) {
  int cid,i;
  printf("==CT==\n");
  for(cid=0;cid<=max_cid;cid++){
   if(CT[cid].array!=NULL) {
     for(i=0;i<CT[cid].size;i++) printf("cid=%d, CT[cid].array[%d]=%d\n",cid,i,CT[cid].array[i]);
   } else { printf("cid=%d, NULL\n",cid); }
  }
}



void print_ST(void) {
  int sid,i;
  printf("==ST==\n");
  for(sid=0;sid<=max_sid;sid++){
   if(ST[sid].array!=NULL) {
     for(i=0;i<ST[sid].size;i++) printf("sid=%d, ST[sid].array[%d]=%d\n",sid,i,ST[sid].array[i]);
   } else { printf("sid=%d, NULL\n",sid); }
  }
}







/*------- fragment length distribition -------*/
/*------- common -------*/

void parse_readlength_range(char* readlength_str){ // only for SE.
   int i, dash_ind=-1;
   for(i=0;i<strlen(readlength_str);i++) if(readlength_str[i]=='-') { dash_ind=i; break; }
   if(dash_ind==-1) { Readlengths.min=atoi(readlength_str); Readlengths.max=Readlengths.min; } // a single number, not a range.
   else { Readlengths.max=atoi(readlength_str+dash_ind+1); readlength_str[dash_ind]='\0'; Readlengths.min=atoi(readlength_str); }
   Fraglengths.max=Readlengths.max;
   Fraglengths.min=Readlengths.min;
   nFraglen = Fraglengths.max - Fraglengths.min +1;
}

void determine_fraglength_range(void){
   Fraglengths.min=Min_Fraglength>readlength?Min_Fraglength:readlength;
   Fraglengths.max=Max_Fraglength>=Fraglengths.min?Max_Fraglength:Fraglengths.min;
   nFraglen = Fraglengths.max - Fraglengths.min +1;
}

void print_FraglengthDist(char *fraglengthfile_name)
{
   int fl_ind;
   FILE * Fraglengthfile;

   Fraglengthfile = fopen(fraglengthfile_name,"w");
   if(Fraglengthfile==NULL){
      fprintf(stderr,"Can't write to fraglength file %s\n",fraglengthfile_name);
      exit(1);
   }

   // print out to file
   fprintf(Fraglengthfile, "Fragment.length\tObs.Counts\tnormalized.Fragment.length.sampling.prob\n");
   for(fl_ind=0;fl_ind<nFraglen;fl_ind++) fprintf(Fraglengthfile, "%d\t%d\t%lg\n",fl_ind+Fraglengths.min,FraglengthCounts[fl_ind+Fraglengths.min],Wf[fl_ind]);

   fclose(Fraglengthfile);
}






/*------- fragment length distribution (light version - default) -------*/

// use observed fraglength count as fraglen sampling probability
void transfer_fraglendist_to_Wf (void){
   int fl_ind;
   double sumWf;
   
   Wf=(double*)malloc(nFraglen*sizeof(double));
   if(Wf==NULL) { fprintf(stderr,"Failed to allocate memory to Wf\n"); exit(1); }

   sumWf=0;
   for(fl_ind=0;fl_ind<nFraglen;fl_ind++) { Wf[fl_ind]=FraglengthCounts[fl_ind+Fraglengths.min]; sumWf += Wf[fl_ind]; }
   for(fl_ind=0;fl_ind<nFraglen;fl_ind++) Wf[fl_ind] /= sumWf;
}



double compute_adjEUMA(int* EUMAarr)
{
   int i;
   double adjustedEUMA=0;
   for(i=0;i<nFraglen;i++) adjustedEUMA += Wf[i]*(double)EUMAarr[i];  /* fragment-length-adjustment */
   return(adjustedEUMA);
}





/*------- positional bias model -------*/

// pos is relative to the coordinate on tid not on seq_array.
char determine_poscategory(int tid, int pos, int fraglen){
   if(posmodel==0) return 0;  // no-bias
   else if(posmodel=='a'){  // arbitrary (primitive version mostly for testing & debugging)
       //fprintf(stderr,"determining poscategory : pos=%d, tid=%d, fraglen=%d\n", pos,tid,fraglen); //DEBUGGING
       if(pos<100) return 1; // 100 bp near 5'end
       if(pos>cuml[tid+1]-cuml[tid-1]-100) return 2; // 100 bp near 3'end
       return 0; // otherwise
   }
   else { fprintf(stderr,"ERROR: incorrect positional bias model.\n"); exit(1); }
}


void normalize_perpos_freq( void ){
   double sumPerpos_5=0;
   double sumPerpos_3=0;
   int i;
   for(i=0;i<perpos_freq_len;i++){
     sumPerpos_5+=perpos_freq_5[i];
     sumPerpos_3+=perpos_freq_3[i];
   }
   for(i=0;i<perpos_freq_len;i++){
     perpos_normfreq_5[i] = perpos_freq_5[i]/(sumPerpos_5-perpos_unavail_freq_5[i]);
     perpos_normfreq_3[i] = perpos_freq_3[i]/(sumPerpos_3-perpos_unavail_freq_3[i]);
   }
}


void determine_scaling_factor_for_perpos_prob (void){
   double impute,b_j; // imputed value
   int i,j,tlen;
   Y=(double*)malloc((max_tid+1)*sizeof(double));
   // impute value for middle range
   for(i=perpos_freq_len-perpos_freq_impute_len;i<perpos_freq_len;i++) impute+=perpos_normfreq_5[i]+perpos_normfreq_3[i];
   impute/=2*perpos_freq_len;

   for(i=0;i<=max_tid;i++){  // i : transcript id
      tlen = cuml[i+1]-cuml[i]-1;
      if(tlen==2*perpos_freq_len) Y[i] = 1;
      else if(tlen>2*perpos_freq_len) Y[i] = (tlen-2*perpos_freq_len)*impute;
      else {
        for(j=0;j<tlen;j++){ // j : pos
          if(j<tlen-perpos_freq_len) b_j = perpos_normfreq_5[j];
        }
      }
   }
}


void print_posbias(char* outfilename_posbias){
  unsigned int i;
  FILE *outfile_posbias;

  outfile_posbias = fopen(outfilename_posbias,"w");
  if(outfile_posbias==NULL){
      fprintf(stderr,"Can't write to output posbias file %s\n",outfilename_posbias);
      exit(1);
  }

  fprintf(outfile_posbias,"relative_position\t5-frag_count\t5-avail_count\t5-norm_frag_count\t3-frag_count\t3-avail_count\t3-norm_frag_count\n");
  for(i=0;i<perpos_freq_len;i++) {
     fprintf(outfile_posbias,"%d\t%f\t%f\t%f\t%f\t%f\t%f\n",i,perpos_freq_5[i],perpos_unavail_freq_5[i],perpos_normfreq_5[i],perpos_freq_3[i],perpos_unavail_freq_3[i],perpos_normfreq_3[i]);
  }

  fclose(outfile_posbias);

}



/*------- seq_array position handling & miscellaneous -------*/

/* return position of its rc (or its fc) */
int flip(int k){  /* k is the index on seq_array, not on sfa */
    return(seqlength - k - readlength);
}

// sf_i and sf_ib return transcript index for a given position on the concatenated sequence (k).
/* assumes that k < borderpos for PE or SS SE. */
int sf_ib(int k){
   int j;
   INTPAIR c;
   c=cumlI[(int)(k/bin)];
   if(c.start==c.end) return(c.start);
   else for(j=c.start;j<=c.end;j++) if(k<cuml[j+1]) return(j);
}

/* does not assumes that k < borderpos for PE or SS SE. */
int sf_i(int k){
   int j;
   INTPAIR c;

   if(k+readlength>borderpos) k = flip(k);  /* not necessary for PE or SS SE, but leave it here temporarily */
   c=cumlI[(int)(k/bin)];
   if(c.start==c.end) return(c.start);
   else for(j=c.start;j<=c.end;j++) if(k<cuml[j+1]) return(j);
}

// sf_pb and sf_p returns position relative to each transcript, given a position on the concatenated sequence (k) and transcript index (i)
/* assumes that k < borderpos for PE or SS SE. */
int sf_pb(int k,int i){
   return(k-cuml[i]);
}
/* does not assume that k < borderpos */
int sf_p(int k,int i){
   if(k+readlength>borderpos) k = flip(k);
   return(k-cuml[i]);
}


// marks all start positions whose substring of length readlength contains a non-ACGT character (eg. 'N' or '@' or '$')
void mark_noncanonical(void)
{
     int i,i2,i2start,i2end;
     int noncanonarr_size = bitsize(seqlength+1);
     noncanonarr = (char*)calloc(noncanonarr_size,sizeof(char));
     if(noncanonarr==NULL) { fprintf(stderr,"Failed to allocate memory to noncanonarr.\n"); exit(1); }

     i=0;
     while(i<=seqlength){
       if(is_noncanonical(seq_array+i,1)) {
          i2start = (i-readlength+1>=0)?(i-readlength+1):0;
          if(i2start<=i){
            for(i2=i2start;i2<=i;i2++) putb(noncanonarr,i2,1);
          }
       }
       i++;
     }

}


int strdiff_pe (char* s1, char* s2, int d) {
   int i;
   for(i=readlength+d;i>=d;i--){  // compare from the last base, because the comparison is on consecutive indices of sorted sfa, the first bases are likely identical.
      if(s1[i]!=s2[i]) return 1;
   }
   for(i=readlength-1;i>=0;i--){  // compare from the last base, because the comparison is on consecutive indices of sorted sfa, the first bases are likely identical.
      if(s1[i]!=s2[i]) return 1;
   }
   return 0; 
}

int strcmp_pe (char* s1, char* s2, int d) {
   int res;
   res = strncmp(s1,s2,readlength);
   if(res==0) return(strncmp(s1+d,s2+d,readlength));
   else return res;
}
int strdiff_se (char* s1, char* s2) {
   int i;
   for(i=readlength-1;i>=0;i--){
      if(s1[i]!=s2[i]) return 1;
   }
   return 0;
}



// return 1 if non-ACGT characters are included in the substring starting at position pos (and length k (eg. k=1 or k=readlength))
// non-ACGT characters include N , $, @ etc.
char is_noncanonical(char* pos, int k){
   int i;
   char base;
   for(i=0;i<k;i++){
     base = *(pos+i);
     if(base!='A' && base!='C' && base!='G' && base!='T') return 1;
   }
   return 0;
}


// both single-end and paired-end
// This function has been rewritten to avoid strtok, which could cause memory problems.
void parse_ensembl_header (char* header_str,char** parsed_header_str)
{
    int i, term=0;
    for(i=0;i<strlen(header_str);i++){
      if(header_str[i]=='\t' || header_str[i]==' ') { (*parsed_header_str)[i]='\0'; term=1; break; }
      else { (*parsed_header_str)[i]=header_str[i]; }
    } 
    if(term==0) (*parsed_header_str)[strlen(header_str)]='\0';
}
void parse_refseq_header(char* header_str, char** parsed_header_str)
{
    int i,numtab=0,j=0,term=0;
    for(i=0;i<strlen(header_str);i++){
      if(header_str[i]=='|'){ 
        numtab++; 
        if(numtab==4) { (*parsed_header_str)[j]='\0'; term=1; break; }
      }      
      else if(numtab==3) { (*parsed_header_str)[j]=header_str[i]; j++; }
    }
    if(term==0) (*parsed_header_str)[j]='\0';
}
char revcomp(char c)
{
   switch(c){
     case 'A': 
     case 'a':
       return 'T';
     case 'C':
     case 'c':
       return 'G';
     case 'G':
     case 'g':
       return 'C';
     case 'T':
     case 't':
       return 'A';
     case '@':
       return '@';
     case '$':
       return '$';
     default:
       return 'N'; 
   }
}

char uc(char c)
{
   switch(c){
     case 'A': 
     case 'a':
       return 'A';
     case 'C':
     case 'c':
       return 'C';
     case 'G':
     case 'g':
       return 'G';
     case 'T':
     case 't':
       return 'T';
     case '@':
       return '@';
     case '$':
       return '$';
     default:
       return 'N'; 
   }

}








/*------- Paired-end suffix array scanning & thread control -------*/

void run_process_mate1_cluster_by_mate_2_PE_threads_2(void) {
 int i;
 int increment,local_MAX_Thread;
 int a[MAX_Thread];

 if(verbose_flag) { fprintf(stdout,"processing mate1 suffix array... (using %d thread(s))\n",MAX_Thread); fflush(stdout); system("date +%m/%d,%T"); }
 //DEBUGGING

 if(MAX_Thread>1) for(i=0;i<MAX_Thread-1;i++) working_thread[i]=0;

 a[0]=0;
 if(MAX_Thread>1) for(i=1;i<MAX_Thread;i++) a[i]=i;

 nDone=0;  // number of finished tids (for progress bar)
 // create MAX_Thread-1 threads
 if(MAX_Thread>1) {
    for(i=0;i<MAX_Thread-1;i++){
       working_thread[i]=1;
       if(verbose_flag>1) { fprintf(stdout,"Creating thread %d...\n",i); fflush(stdout); }
       pthread_create(&pth[i],NULL,process_mate1_cluster_by_mate_3,(void*)(a+i));
    }
 }
 if(verbose_flag) { fprintf(stdout,"Running the unthreaded part...\n");; fflush(stdout); }
 process_mate1_cluster_by_mate_3((void*)(a+MAX_Thread-1)); // This part is unthreaded. (when MAX_Thread==1, only this part is run)

 // Wait until all the threads finish.
 if(MAX_Thread>1)
   for(i=0;i<MAX_Thread-1;i++) {
       pthread_join(pth[i],NULL);
       working_thread[i]=0;
   }

  max_cid = max_tid + rsh_size;
  
  if(verbose_flag>1) fprintf(stdout,"\r%3d%% done...", 100); // must be at least 1% increase to be printed.
}



void* process_mate1_cluster_by_mate_3(void* a){
   int sfa_mod = *((int*)a);
   int i,j,k,start,end,d,sfa_s_size,cmpres;
   int *sfa_s,*sfa_d;
   int nDone_incre=0;
   int sfa_size_1000 = sfa_size/1000;
   int min_tid,sfi;

   start=0;
   do {
    i=start; min_tid=sf_i(sfa[i]);
    do{ 
      if((sfi = sf_i(sfa[i]))<min_tid) min_tid=sfi; 
      i++; 
    } while(i<sfa_size && getb(sfa_m,i)==0);
    end=i-1;
    if(min_tid%MAX_Thread==sfa_mod){
  
      // create sfa_s and sfa_d arrays with identical mate1 and all possible mate2
      //if(MAX_Thread>1) pthread_mutex_lock(&mutex1); //DEBUGGING
      //fprintf(stderr,"creating mate2 arrays..(thread %d)\n", sfa_mod); //DEBUGGING
      sfa_s_size = nFraglen*(end-start+1);
      //fprintf(stderr,"sfa_s_size: %d, nFraglen= %d, end=%d, start=%d\n",sfa_s_size,nFraglen,end,start); //DEBUGGING
      if(MAX_Thread>1) pthread_mutex_lock(&mutex1); 
      sfa_s = (int*)malloc(sfa_s_size*sizeof(int));
      sfa_d = (int*)malloc(sfa_s_size*sizeof(int));
      if(sfa_s==NULL) { fprintf(stderr,"Failed to allocate memory to sfa_s.\n"); exit(1); }
      if(sfa_d==NULL) { fprintf(stderr,"Failed to allocate memory to sfa_d.\n"); exit(1); }
      if(MAX_Thread>1) pthread_mutex_unlock(&mutex1); 

      k=0;
      for(j=start;j<=end;j++)
       for(d=Fraglengths.min-readlength;d<=Fraglengths.max-readlength;d++){
         if(sfa[j]+d<seqlength && sf_i(sfa[j]+d)==sf_i(sfa[j]) && getb(noncanonarr,sfa[j]+d)==0 && !(sfa[j]<borderpos && sfa[j]+d>borderpos)){
           //if(j==0 && d==200) fprintf(stderr, "included: sfa[%d]=%d, d=%d, noncan=%d, seqlen=%d\n",j,sfa[j],d,getb(noncanonarr,sfa[j]+d),seqlength);
           if(library_strand_type!=0){ //strand-specific
             sfa_s[k]=sfa[j]+d; //mate2 start position 
             sfa_d[k]=d;
             k++;
           }else{ //unstranded
             cmpres = strcmp_pe(seq_array+sfa[j],seq_array+flip(sfa[j]+d),d);
             if((sfa[j]<borderpos && cmpres<=0) || (sfa[j]>borderpos && cmpres<0)){
               sfa_s[k]=sfa[j]+d; //mate2 start position 
               sfa_d[k]=d;
               k++;
             }
           }
         }
         //else fprintf(stderr,"excluded: sfa[%d]=%d d=%d\n",j,sfa[j],d); //DEBUGGING
       }
       if(k>0){
         if(k<sfa_s_size) { 
           sfa_s_size=k; 
           sfa_s=realloc(sfa_s,sfa_s_size*sizeof(int));
           if(sfa_s==NULL) { fprintf(stderr,"Failed to reallocate memory to sfa_s (size=%d).\n",sfa_s_size); exit(1); }
           sfa_d=realloc(sfa_d,sfa_s_size*sizeof(int));
           if(sfa_d==NULL) { fprintf(stderr,"Failed to reallocate memory to sfa_d (size=%d).\n",sfa_s_size); exit(1); }
         }
   
         // DEBUGGING
         /*
         for(k=0;k<sfa_s_size;k++){
           fprintf(stderr,"pos_on_seq_array=%d, tid=%d, pos=%d, d=%d, ",sfa_s[k],sf_i(sfa_s[k]),sf_p(sfa_s[k],sf_i(sfa_s[k])),sfa_d[k]);
           for(j=0;j<readlength;j++) fprintf(stderr,"%c",seq_array[sfa_s[k]+j]); //mate2
           fprintf(stderr,"\n");
         }
         */ 
       
   
         // sort the sfa_s and sfa_d arrays by mate2
         //fprintf(stderr,"sorting by mate2..(thread %d)\n", sfa_mod); //DEBUGGING
         quick_sort_suffixarray_by_mate_2(&sfa_s,&sfa_d,0,sfa_s_size-1);
   
         // DEBUGGING
         /*
         for(k=0;k<sfa_s_size;k++){
           fprintf(stderr,"k=%d, tid=%d, pos=%d, d=%d, ",k, sf_i(sfa_s[k]),sf_p(sfa_s[k],sf_i(sfa_s[k])),sfa_d[k]);
           //for(j=0;j<readlength;j++) fprintf(stderr,"%c",seq_array[sfa_s[k]+j]); //mate2
           //fprintf(stderr,"\n");
           //for(j=0;j<readlength;j++) fprintf(stderr,"%c",seq_array[sfa_s[k]+j-sfa_d[k]]); //mate1
           fprintf(stderr,"\n");
         } 
         */
       
   
         // construct rsh array 
         //fprintf(stderr,"constructing rshbucket..(thread %d)\n", sfa_mod); //DEBUGGING
         if(posmodel==0) construct_rshbucket_PE_3(&sfa_s,&sfa_d,sfa_s_size);
         else construct_rshbucket_PE_3_posbias(&sfa_s,&sfa_d,sfa_s_size);
         //fprintf(stderr,"finished constructing rshbucket..(thread %d)\n", sfa_mod); //DEBUGGING
       }   


       // free the new arrays
       free(sfa_s); sfa_s=NULL;
       free(sfa_d); sfa_d=NULL;
   
       //progress bar
       if(verbose_flag>1){
        nDone_incre+=end-start+1; 
        if(nDone_incre>30000) { 
           if(MAX_Thread>1) pthread_mutex_lock(&mutex1); 
           nDone+=nDone_incre; //progress bar
           fprintf(stdout,"\r%3.2lf%% done...", (double)nDone/(double)sfa_size_1000/10);
           if(MAX_Thread>1) pthread_mutex_unlock(&mutex1); 
           nDone_incre=0;
        }
       }
    }
    start=end+1;
   }while(start<sfa_size);
}









/* -------------- EMSAR sample run MLE for FPKM & iReadcount)-------------------- */

double Fp (inta C){
    double sum=0;
    double lamb,logP;
    int i,cid;
    for(i=0;i<C.size;i++){
        cid=C.array[i];
        if(EUMAps[cid]==0) continue;  /* skip if EUMA=0 */
        lamb = lambdap(cid);
        if(lamb == 0){
          if(ReadCount[cid] == 0)  logP = 0;   /* lamb==0 & R ==0 */
          else return(NEAR_LOWEST_NUMBER);                      /* lamb==0 & R!=0 */
        }
        else if(lamb < 0) return(NEAR_LOWEST_NUMBER);    /* lamb<0 (likelihood goes -Inf to avoid negative estimates)*/
        else logP = (double)ReadCount[cid] * log(lamb) - lamb ;
        sum+= logP;
    }
    if(sum<NEAR_LOWEST_NUMBER) sum=NEAR_LOWEST_NUMBER; /* very low but not as low as LOWEST_NUMBER */
    return sum;
}

double lambdap (int cid){
    double sum=0;
    int i,tid;
    for(i=0;i<CT[cid].size;i++){
       tid = CT[cid].array[i];
       if(FPKM[tid]<0) return -1;  /* contraint that all FPKM values must be >=0. This return will cause the F value to be NEAR_LOWEST_NUMBER. */
       sum+= FPKM[tid];
    }
    return EUMAps[cid]*sum;
}

void run_MLE_threads(void) {
 int i;
 int sid_increment,sid_division_starts[MAX_Thread],local_MAX_Thread;

 if(MAX_Thread>1) for(i=0;i<MAX_Thread-1;i++) working_thread[i]=0;

 /* Divide sids to roughly equal number of groups */
 sid_increment = max_sid/MAX_Thread+1; // add 1 to make sure the last group contains no more than sid_increment sid's.
 sid_division_starts[0]=0;
 local_MAX_Thread=MAX_Thread;
 if(MAX_Thread>1) for(i=1;i<MAX_Thread;i++) {
    sid_division_starts[i]=sid_division_starts[i-1]+sid_increment;
    if(sid_division_starts[i]>max_sid) { local_MAX_Thread = i; break; }
 }

 nDone=0;  // number of finished tids (for progress bar)
 // create MAX_Thread-1 threads
 if(local_MAX_Thread>1) {
    for(i=0;i<local_MAX_Thread-1;i++){
       working_thread[i]=1;
       pthread_create(&pth[i],NULL,MLE_range,sid_division_starts+i);
    }
 }
 MLE_range(sid_division_starts+local_MAX_Thread-1); // This part is unthreaded. (when MAX_Thread==1, only this part is run)

 // Wait until all the threads finish.
 if(local_MAX_Thread>1)
   for(i=0;i<local_MAX_Thread-1;i++) {
       pthread_join(pth[i],NULL);
       working_thread[i]=0;
   }

 if(verbose_flag>1) fprintf(stdout,"\r%3d%% done...",100);

}

void* MLE_range(void* arg) {
  int sid;
  int sid_start = *((int*)arg);
  int sid_end = sid_start + max_sid/MAX_Thread; //consistent with function run_MLE_Thread (end = start+increment-1)
  if(sid_start<=max_sid){
    if(sid_end>max_sid) sid_end=max_sid;
    for(sid=sid_start;sid<=sid_end;sid++){
       MLE(sid);
       if(MAX_Thread>1) pthread_mutex_lock(&mutex1);
       nDone+=ST[sid].size; MLE_progressbar();
       if(MAX_Thread>1) pthread_mutex_unlock(&mutex1);
    }
  }
}

void MLE_progressbar(void){
  if(verbose_flag>1) fprintf(stdout,"\r%3d%% done...",(nDone*100)/(max_tid+1));
}
/* option 'l' vs 'r' no longer exists */
/* returns 1 */
int MLE (int sid){
    inta C = SC[sid];
    double F,maxF,premaxF;
    int i,j,tid,cid;
    int sum;

    double stepSize[ST[sid].size];
    double init_max=100;
    double init_min=0;
    double init_stepSize=100;
    double acc=2.0;
    double max_stepSize;
    double candidate[5];
    int best;
    
    int nIter;
    int notconverged;
    int nLoop;


    /* in case of no read count at all in all of the cid's in S, all the FPKM values are 0. */
    sum=0;
    for(i=0;i<C.size;i++) sum+=ReadCount[C.array[i]];
    if(sum==0) {
       for(i=0;i<ST[sid].size;i++) { tid=ST[sid].array[i]; FPKM[tid]=0; }
       return(1);
    }

    /* For the obvious case where the single tid shares no sequences with any other transcripts */
    if(C.size==1 && ST[sid].size==1) {  // these may be cases in which C.size>1 but ST[sid].size=1 because one cid consists of repetion of a single tid (internal repeat).
       tid=ST[sid].array[0];
       cid=C.array[0];
       FPKM[tid]=(double)ReadCount[cid]/EUMAps[cid];
    }
    else {

      nLoop=0;
      do {
  
        /* number of iterations, if this is too large, then go back and choose a different initial values */
        nIter=0;
        notconverged=0; 
    
        /* random initial point betweeen init_min and init_max (FPKM = 'currentPoint') */
        for(i=0;i<ST[sid].size;i++) {
          tid = ST[sid].array[i];
          FPKM[tid] = rand()/(RAND_MAX+1.0)*(init_max-init_min)+init_min;
        }
    
        /* initial step size */
        for(i=0;i<ST[sid].size;i++) stepSize[i] = init_stepSize;
    
        do {
    
           nIter++;
           premaxF = Fp(C);  //'before'
    
           max_stepSize=0;
           for(i=0;i<ST[sid].size;i++){
             tid = ST[sid].array[i];
             if(stepSize[i]<CONVERGENCE_EPSILON_STEPSIZE) continue;
             candidate[0]=0; candidate[1]=-stepSize[i]*acc*2; candidate[2]=-stepSize[i]; candidate[3]=stepSize[i]; candidate[4]=stepSize[i]*acc*2;
             /* the first one must be 0, so that in case all five has a tie, the first one is selected. */
    
             best = -1;
             maxF = LOWEST_NUMBER; //'bestscore', initialized to the worst score so that anything can beat it.
             for(j=0;j<5;j++){
                FPKM[tid]=FPKM[tid]+candidate[j];
                F=Fp(C);  // 'temp'
                FPKM[tid]=FPKM[tid]-candidate[j];
                if(F>maxF){
                   maxF=F;
                   best=j;
                }
             }
             if(best==0){  stepSize[i]=stepSize[i]/acc;  } // candidate[best]==0
             else if(best==2||best==3) { FPKM[tid]=FPKM[tid]+candidate[best]; }
             else {
                FPKM[tid]=FPKM[tid]+candidate[best];
                stepSize[i]=stepSize[i]*acc;
             }
    
             if(stepSize[i] > max_stepSize) max_stepSize = stepSize[i];
           }
    
           if(nIter> MAX_NITER_MLE) { notconverged=1; fprintf(stdout,"FPKM not converging, reinitializing.. (sid=%d size=%d)\n",sid,ST[sid].size); break; }
        }while(Fp(C) - premaxF >= CONVERGENCE_EPSILON || max_stepSize > CONVERGENCE_EPSILON_STEPSIZE);
        nLoop++;
      }while(notconverged==1 && nLoop < MAX_NLOOP_MLE);

    }

    return(1);
}


void construct_FPKMfinal (int round){
   int tid;

   FPKMfinal[round]=(double*)malloc((max_tid+1)*sizeof(double));
   if(FPKMfinal[round]==NULL) { fprintf(stderr,"Failed to allocate memory to FPKMfinal[round=%d].\n",round); exit(1); }
   for(tid=0;tid<=max_tid;tid++) FPKMfinal[round][tid]=FPKM[tid];  //tid=cid for single-tid cid's.
}


void delete_FPKMfinal(void)
{
   int round;
   for(round=0;round<NUM_ROUND;round++){
     free(FPKMfinal[round]); FPKMfinal[round]=NULL;
   }
   free(FPKMfinal); FPKMfinal=NULL;
}


void construct_EUMAps (void){
   int cid;
   EUMAps=(double*)malloc((max_cid+1)*sizeof(double));
   if(EUMAps==NULL) { fprintf(stderr,"Failed to allocate memory to EUMAps.\n"); exit(1); }
   for(cid=0;cid<=max_cid;cid++) EUMAps[cid] = adjEUMA[cid] / 1E3 * ((double)TotalReadCount/1E6) * pow(10,DELTA);

}

void delete_index_table (void){
  int tid;
  for(tid=0;tid<=max_tid;tid++) { free(IndexTable[tid]); IndexTable[tid]=NULL; }
  free(IndexTable);
}


void print_FPKMfinal (char* filename){
   int round,tid;
   double sum,mean,sqsum,sd;
   FILE* file;
   double iReadcount,totalFPKM;
   int iReadcount_int, total_iReadcount=0;

    file=fopen(filename,"w");
    if(file == NULL){
      fprintf(stderr,"Can't write to FPKMfile %s",filename);
      exit(1);
    }

   totalFPKM=0;
   for(tid=0;tid<=max_tid;tid++) {
     for(round=0;round<NUM_ROUND;round++){
       totalFPKM+=FPKMfinal[round][tid]/NUM_ROUND;
     }
   }
  
   
   fprintf(file,"transcriptID\tFPKM\tsd.of.FPKM\teff.length\tiReadcount\tiReadcount.int\tTPM\n");

   for(tid=0;tid<=max_tid;tid++) {

     /* mean */
     sum=0;
     for(round=0;round<NUM_ROUND;round++) {
       sum+=FPKMfinal[round][tid];
     }
     mean = sum/NUM_ROUND;

     /* SD */
     sqsum=0;
     for(round=0;round<NUM_ROUND;round++) {
       sqsum+=pow((FPKMfinal[round][tid]-mean),2);
     }
     sd = sqrt( sqsum/(NUM_ROUND-1) ) / NUM_ROUND;  /* SD of Mean */

     /* inferred readcount */
     iReadcount = (iEUMA[tid]/1E3)*mean*((double)TotalReadCount/1E6);
     iReadcount_int = Round_off(iReadcount);
     total_iReadcount+=iReadcount_int;

     fprintf(file,"%s\t%lf\t%lf\t%lf\t%lf\t%d\t%lf\n",IndexTable[tid],mean,sd,iEUMA[tid],iReadcount,iReadcount_int,mean*1E6/totalFPKM);
   }

   fclose(file);
   if(verbose_flag>0) fprintf(stdout,"Total inferred readcount=%d\n",total_iReadcount);
}


int Round_off(double x){
   return( x-(int)x>=0.5?(int)x+1:(int)x );
}
void compute_iEUMA(void)
{
   int i,cid,tid;

   iEUMA = (double*)calloc((max_tid+1),sizeof(double)); //total effective length for each transcript
   if(iEUMA==NULL) { fprintf(stderr,"Failed to allocate memory to iEUMA.\n"); exit(1); }   
   for(cid=0;cid<=max_cid;cid++){
     if(CT[cid].array!=NULL) {
       for(i=0;i<CT[cid].size;i++) {
           tid=CT[cid].array[i];
           iEUMA[tid]+=adjEUMA[cid];
       }
     }
   }
}








/*------- Overall preprocessing (EUMA) control -------*/

void preprocess_SE(int taglen) {
    int i,j,numtags;
    char** tag;
    SFA_RANGE a;

    numtags = generate_seqtag(&tag,taglen);

    for(readlength=Readlengths.min; readlength<=Readlengths.max ; readlength++){

         /*-- generate suffixarray --*/

         if(verbose_flag>0) { fprintf(stdout, "Generating suffix array for read length %d...\n", readlength); fflush(stdout); fflush(stdout); system("date +%m/%d,%T"); }
         determine_sfa_size();

         for(i=0;i<numtags;i++){

           if(verbose_flag>0) { fprintf(stdout, "initializating suffix array... (w/ tag %s) :",tag[i]); fflush(stdout); system("date +%m/%d,%T"); }
           if(library_strand_type!=0) initialize_suffixarray_SS_4(tag[i]); //stranded
           else initialize_suffixarray_NS_5(tag[i]); //unstranded
           if(verbose_flag>0) fprintf(stdout,"sfa_size=%d (%d~%d)\n",sfa_size,local_sfa_start,local_sfa_end);

           //fprintf(stderr,"==sfa==\n");for(j=0;j<sfa_size;j++) fprintf(stderr,"sfa[%d]=%d\n",j,sfa[j]); // DEBUGGING

           if(verbose_flag>0) { fprintf(stdout, "sorting suffix array... :"); fflush(stdout); system("date +%m/%d,%T"); }
           if(MAX_Thread>1) for(j=0;j<MAX_Thread-1;j++) working_thread[j]=0;
           a.left=local_sfa_start; a.right=local_sfa_end; quick_sort_suffixarray_4((void*)&a);
           //a.left=local_sfa_start; a.right=local_sfa_end; quick_sort_suffixarray_4_nothread((void*)&a);  // no threading, temporary replacement
           if(verbose_flag>0) { fprintf(stdout, "sorting suffix array finished... :"); fflush(stdout); system("date +%m/%d,%T"); }

           if(print_sfa_flag==1) { print_sfa(); }


           /* generating rsh array */
           if(verbose_flag>0) { fprintf(stdout, "generating rsh array by scanning suffix array... :");fflush(stdout); system("date +%m/%d,%T"); }
           if(posmodel==0) construct_rshbucket_2(readlength); //readlength is the same as fragment length for single-end.
           else construct_rshbucket_2_posbias(readlength); // with positional bias model.

           if(verbose_flag>0) { fprintf(stdout, "deleting suffix array... :");fflush(stdout); system("date +%m/%d,%T"); }
           free(sfa); sfa=NULL;
         }

    }
    delete_seqtag(&tag, numtags);

    if(verbose_flag) { fprintf(stdout, "freeing sequence array... :"); fflush(stdout); system("date +%m/%d,%T"); }

    free(seq_array); seq_array=NULL;
}



void preprocess_PE(int taglen){
    int numtags,i,j,d,dmin,dmax;
    char** tag;
    SFA_RANGE a;

    numtags = generate_seqtag(&tag,taglen);

    if(verbose_flag>0) { fprintf(stdout, "Marking noncanonical positions...\n"); fflush(stdout); system("date +%m/%d,%T"); }
    mark_noncanonical();
    //for(i=0;i<seqlength;i++) fprintf(stderr,"%d",getb(noncanonarr,i)); fprintf(stderr,"\n"); //DEBUGGING

    /*-- generate suffixarray --*/
    if(verbose_flag>0) { fprintf(stdout, "Generating suffix array for read length %d...\n", readlength); fflush(stdout); fflush(stdout); system("date +%m/%d,%T"); }
    determine_sfa_size();

    for(i=0;i<numtags;i++){

      if(verbose_flag>0) { fprintf(stdout, "initializating suffix array... (w/ tag %s) :",tag[i]); fflush(stdout); system("date +%m/%d,%T"); }
      if(library_strand_type!=0) initialize_suffixarray_SS_PE_2(tag[i]); //stranded
      else initialize_suffixarray_NS_PE_2(tag[i]); //unstranded
      if(verbose_flag>0) fprintf(stdout,"sfa_size=%d (%d~%d)\n",sfa_size,local_sfa_start,local_sfa_end);

      //fprintf(stderr,"==sfa test==\n"); for(j=local_sfa_start;j<=local_sfa_end;j++) if(seq_array+sfa[j]=='*'||seq_array+sfa[j]+readlength-1=='*') fprintf(stderr,"lala\n"); // DEBUGGING

      if(verbose_flag>0) { fprintf(stdout, "sorting suffix array... :"); fflush(stdout); system("date +%m/%d,%T"); }
      if(MAX_Thread>1) for(j=0;j<MAX_Thread-1;j++) working_thread[j]=0;
      a.left=local_sfa_start; a.right=local_sfa_end; quick_sort_suffixarray_4((void*)&a);
      if(verbose_flag>0) { fprintf(stdout, "sorting suffix array finished... :"); fflush(stdout); system("date +%m/%d,%T"); }

    }
    delete_seqtag(&tag, numtags);

    // sometimes some parts of the sequence ('N') are skipped so the actual sfa size may be smaller than max_sfa_size. In this case, readjust the size.
    if(local_sfa_end+1!=sfa_size) { 
       sfa_size = local_sfa_end+1;
       sfa = realloc(sfa,sfa_size*sizeof(int));
       if(sfa==NULL) { fprintf(stderr,"Failed to reallocate memory to suffix array.\n"); exit(1); }
       if(verbose_flag>0) fprintf(stdout, "suffix array size readjusted to %d...\n", sfa_size); 
    }

    if(verbose_flag>0) { fprintf(stdout, "generating mark array for mate1 suffix array...\n"); fflush(stdout); system("date +%m/%d,%T"); } 
    mark_sfa_se();

    if(print_sfa_flag==1) { print_sfa(); }

    //fprintf(stderr,"In proprocess_PE: nFraglen=%d\n",nFraglen); //DEBUGGING
    run_process_mate1_cluster_by_mate_2_PE_threads_2();

    free(sfa); sfa=NULL;
    free(sfa_m); sfa_m=NULL;
    free(noncanonarr); noncanonarr=NULL;

    if(verbose_flag>0) fprintf(stdout, "freeing sequence array... :");
    free(seq_array); seq_array=NULL;
}








/*------- usage -------*/

void printusage(char* program)
{
   printf("Usage : %s <options> -x fastafile outdir outprefix alignmentfile|alignmentfilelist\n",program);
   printf("Usage2 : %s <options> -I rshfile outdir outprefix alignmentfile|alignmentfilelist\n",program);
   printf("Usage3 : bowtie command | %s <options> [-x fastafile][-I rshfile] outdir outprefix\n\n",program);
   printf("\tex : %s -p 4 -h R -x human.rna.fna RNAseq sample22 sample22.bowtieout\n",program);
   printf("\tex2 : %s -p 4 -B --PE -R -x human.rna.fna RNAseq sample22 sample22.BAM\n",program);
   printf("\tex3: %s -M -p 16 -B -I human.rna.rsh RNAseq samples samples.BAMlist   ## BAM list file with -M and -B options\n",program);
   printf("\tex4 : bowtie -v 2 -a -m 100 -p 4 human.rna sample22.fastq | %s human.rna.fna RNAseq sample22\n\n",program);

   printf("\tInput files\n\n");

   printf("\t  *Either a fasta file (-x) or an rsh file (-I) must be provided.\n");
   printf("\t   fasta file : A transcriptome sequence library file in fasta format. The header can either be in Refseq format('>xx|xx|xx|name|xx') in which name or xx does not contain '|' letter, or in Ensembl format '>name' or '>name xx' in which name does not contain whitespace(space or tab) and xx is separated by whitespace from the name. Name must be a unique identifier of a sequence.\n");
   printf("\t   rsh file : This file can be constructed by using -R option combined with -x fasta file. Once it is constructed, it can be used for other RNA-seq samples with the same read length for speedy calculation.\n\n");
   printf("\t  *Either an alignment file or a file listing alignment files must be provided.\n");
   printf("\t   alignmentfile : either a default bowtie output format (bowtieoutfile) or SAM/BAM format. SAM and BAM files must be used with -S and -B options, respectively. The SAM/BAM file should be either produced with the default sorting by bowtie or sorted by qname. We strongly recommend that bowtie is run with options -a -v 2 -m 100, without --best or --strata, without allowing any indels, for best results with EMSAR. The bowtie/SAM/BAM files can be streamed directly to EMSAR through a pipe. If mismatches were allowed, it is highly recommended that the SAM/BAM files contain the auxiliary MD flags that contains the match/mismatch information. Bowtie by default produces SAM/BAM files with the MD field.\n");
   printf("\t   alignmentfilelist : a text file containing one alignment file name per line. The alignment files can be one of the formats specified above. When a list file is used, -M (--multisample) option must be used along with the -B, -S options. All the files must be of the same read length and library type. Paired and single-end samples cannot be combined.\n\n");
  
   printf("\tOptions\n");
   printf("\t  -M, --multisample : multisample (default : single sample). When this option is used, a file containing a list of alignment files must be specified as input instead of an alignment file. The resulting expression values may be slightly different between multi-sample and single-sample runs, because a pooled fragment length distribution is used for a multi-sample run.\n");
   printf("\t  -P, --PE : paired-end data (default : single-end)\n");
   printf("\t  -s, --strand_type <strand_type> : set strand type ('ns','ssf','ssr' for single-end, 'ns','ssfr','ssrf' for paired-end). (default: ns(unstranded))\n");
   printf("\t  -S, --SAM : input file format is SAM (by default, default bowtie output)\n");
   printf("\t  -B, --BAM : input file format is BAM (by default, default bowtie output)\n");
   printf("\t  -R, --print_rsh : create an rsh file from the current fasta file (must be combined with -x not -I) and the bam file (that determines read length). This rsh file can be reused for later samples for fast calculation.\n");
   printf("\t  -p, --maxthread <num_thread>: number of threads to use simultaneously. Using multiple threads comes with a slight increase in memory usage. We recommend -p 4 for an optimal performance. (default : 1)\n");
   printf("\t  -F, --maxfraglen <Max_fraglen> : Maximum fragment length. Use a number that is safely large, unless you want to apply fragment length filtering. Default 400. Not applicable for SE.\n");
   printf("\t  -f, --minfraglen <Min_fraglen> : Minimum fragment length. Use a number that is safely small, unless you want to apply fragment length filtering. Default 1. Not applicable for SE\n");
   printf("\t  -h, --header <E|R> : fasta header option. E : Ensembl header(default), R : RefSeq header.\n");
   printf("\t  -k, --max_repeat <max_repeat> : the maximum number of alignments per read allowed. Reads exceeding this number are discarded.\n");
   printf("\t  -g, --print_segments : print out segment information.\n");
   printf("\t  -T, --print_sfa : print out (mate1) suffix array.\n");
   printf("\t  -v --verbose : make a (very) verbose output log.\n");
   printf("\t  -q --no_verbose : turn-off verbosity.\n\n");

   printf("\tAdvanced Options\n");
   printf("\t  -b, --binsize <binsize> : binsize for indexing transcriptome (default : 5000). This bin size can be made smaller to improve speed at the cost of more memory, or vise versa.\n");
   printf("\t  -t, --taglen <taglen> : length of short sequence tags to use for constructing suffix array on a subset of substrings only. This affects only speed and memory usage. Currently, three values are supported (1,2,3). (default 2)\n");
   printf("\t  -n, --nround <num_rounds> : number of MLE runs for computing mean and sd of FPKM. Default 4.\n");
   printf("\t  -e, --epsilon <epsilon> : epsilon value to check convergence of MLE (default : 1E-9)\n");
   printf("\t  -r, --precision <precision> : estimate precision (epsilon for step size) (default : 1E-15)\n");
   printf("\t  -i, --max_niter_mle <max_niter_mle> : maximum number of iteration before determining not converging and reinitializing for MLE (default : 10000)\n");
   printf("\t  -d, --delta <delta> : delta offset for MLE (default : 0)\n");
}


void printusage_build(char* program)
{
   printf("Usage : %s <options> fastafile readlength(range) outdir outprefix\n",program);
   printf("\tex (SE, unstranded, readlength 76 bp) : %s Hsa.GRCh37.65.gtf.fa 76 rsh human.rna.fna.SE.l76\n",program);
   printf("\tex (SE, forward-stranded, readlength 50-60 bp) : %s -s ssf Hsa.GRCh37.65.gtf.fa 50-60 rsh human.rna.fna.SE.l50-60.ssf\n",program);
   printf("\tex (PE, unstranded, readlength 101 bp, fraglength range 1-400 bp, use 4 threads) : %s --PE -p 8 Hsa.GRCh37.65.gtf.fa 101 rsh human.rna.fna.PE.l101.F1-400\n",program);
   printf("\tex (PE, unstranded, readlength 101 bp, fraglength range 1-500 bp, Refseq format, use 8 threads) : %s --PE -h R -p 8 -F 500 human.rna.fna 101 rsh human.rna.fna.PE.l101.F1-500\n",program);


   printf("\tInput files\n\n");

   printf("\t   fasta file : A transcriptome sequence library file in fasta format. The header can either be in Refseq format('>xx|xx|xx|name|xx') in which name or xx does not contain '|' letter, or in Ensembl format '>name' or '>name xx' in which name does not contain whitespace(space or tab) and xx is separated by whitespace from the name. Name must be a unique identifier of a sequence.\n\n");
  
   printf("\tOptions\n");
   printf("\t  -P, --PE : paired-end data (default : single-end)\n");
   printf("\t  -s, --strand_type <strand_type> : set strand type ('ns','ssf','ssr' for single-end, 'ns','ssfr','ssrf' for paired-end). (default: ns(unstranded))\n");
   printf("\t  -p, --maxthread <num_thread>: number of threads to use simultaneously. Using multiple threads comes with a slight increase in memory usage. We recommend -p 4 for an optimal performance. (default : 1)\n");
   printf("\t  -F, --maxfraglen <Max_fraglen> : Maximum fragment length. Use a number that is safely large, unless you want to apply fragment length filtering. Default 400. Not applicable for SE.\n");
   printf("\t  -f, --minfraglen <Min_fraglen> : Minimum fragment length. Use a number that is safely small, unless you want to apply fragment length filtering. Default 1. Not applicable for SE\n");
   printf("\t  -h, --header <E|R> : fasta header option. E : Ensembl header(default), R : RefSeq header.\n");
   printf("\t  -k, --max_repeat <max_repeat> : the maximum number of alignments per read allowed. Reads exceeding this number are discarded.\n");
   printf("\t  -T, --print_sfa : print out (mate1) suffix array.\n");
   printf("\t  -v --verbose : make a (very) verbose output log.\n");
   printf("\t  -q --no_verbose : turn-off verbosity.\n\n");

   printf("\tAdvanced Options\n");
   printf("\t  -b, --binsize <binsize> : binsize for indexing transcriptome (default : 5000). This bin size can be made smaller to improve speed at the cost of more memory, or vise versa.\n");
   printf("\t  -t, --taglen <taglen> : length of short sequence tags to use for constructing suffix array on a subset of substrings only. This affects only speed and memory usage. Currently, three values are supported (1,2,3). (default 2)\n");
}
