#include "readgenerator.h"

/* read raw fasta file (can be multi-line) and store it as a concatenated seq_array and creates an index_table file. */
void generate_reads(char* fastafilename, int readlength, int fraglength, int numreads,char strand_specific_option, char* header_prefix, char* outfilepath, char* outfilepath2, char pe)
{
  FILE* fastafile=0;
  FILE *outfile=0, *outfile2=0;
  int fasta_size=INIT_FASTA_SIZE;
  int i,j,jr,jf;
  int c,prev_c;
  char mode;
  char *seq_array;
  int seqlength,pos;

  fastafile = fopen((const char*)fastafilename,"r");
  /*Check for validity of the file.*/
  if(fastafile == 0)
  {
     printf("can't open fasta file.\n");
     exit(1);
  }

  j=0; // position on seq_array.
  mode='h'; // start with a header mode.
  prev_c=0;

  seq_array=(char*)malloc(fasta_size*sizeof(char));
  fprintf(stdout,"initialized seq_array..\n");  //DEBUGGING


  while((c = fgetc(fastafile)) != EOF) {
     if(prev_c==0 && c!='>') { fprintf(stderr,"ERROR: wrong fasta file format.\n"); exit(1); } //firstline not starting with '>'
     if(c=='>') {
        mode='h'; // header
        if(j!=0) { 
           seq_array[j]='@'; // delimiter (getting ready for new sequence)
           j++; 
           if(j>=fasta_size) {
             fasta_size+=FASTA_SIZE_ADD; 
             seq_array = realloc(seq_array,fasta_size*sizeof(char)); 
           }
        }
     }
     else if(prev_c=='\n' && mode=='h') mode='s'; //sequence

     //fprintf(stdout,"current mode = %c, i=%d,j=%d, c=%c, max_nm_index=%d, fasta_size=%d, header_str=%s\n",mode,i,j,c,max_nm_index,fasta_size,header_str);  //DEBUGGING
     if(mode=='s'){
       if(c!='\n' && c!=' ' && c!='\t') { 
          seq_array[j]=c; 
          j++;
          if(j>=fasta_size) {
            fasta_size+=FASTA_SIZE_ADD; 
            seq_array = realloc(seq_array,fasta_size*sizeof(char));
          }
       }
     }
     prev_c=c;
  }
  fclose(fastafile);

  // rc
  if(strand_specific_option==0){
    seq_array[j]='$'; // delimiter between fw and rc.
    fasta_size= (j+1)*2;
    seq_array = realloc(seq_array,(fasta_size)*sizeof(char));
  
    jr=j+1; 
    for(jf=j-1;jf>=0;jf--){
      seq_array[jr]=revcomp(seq_array[jf]); jr++;
    }
    seq_array[jr]='$';  // the concatenated sequence is in the format: f0@f1@f2@f3$r3@r2@r1@r0$
    seq_array[jr+1]='\0';
    seqlength=jr; // last '$' position.
  }
  else {
    seq_array[j]='$';
    seq_array[j+1]='\0';
    fasta_size=j+1;
    seq_array = realloc(seq_array,(fasta_size)*sizeof(char));
    seqlength=j;
  }
  //seq_array is done
  //fprintf(stdout,"seqarray=\n%s\n",seq_array);  //DEBUGGING

  // generate and print reads
  outfile = fopen((const char*)outfilepath,"w");
  if(pe) outfile2 = fopen((const char*)outfilepath2,"w");

  for(i=0;i<numreads;i++){
    if(numreads>100 && i%(numreads/100)==0) fprintf(stderr,"\r%3d%% done...",i/(numreads/100));  /* progress log */
    else if(i==numreads-1) fprintf(stderr,"\r%3d%% done...",100); /* progress log*/

    do { pos=rand()%(seqlength-fraglength+1); }while(checkvalid(seq_array+pos,fraglength)==0);
    if(pe) {
       fprintf(outfile, ">%s%d/1\n",header_prefix,i);
       fprintf(outfile2, ">%s%d/2\n",header_prefix,i);
       for(j=0;j<readlength;j++) fprintf(outfile,"%c",seq_array[pos+j]);
       for(j=readlength-1;j>=0;j--) fprintf(outfile2,"%c",revcomp(seq_array[pos+fraglength-readlength+j]));
       fprintf(outfile,"\n");
       fprintf(outfile2,"\n");
    }
    else {
       fprintf(outfile, ">%s%d\n",header_prefix,i);
       for(j=0;j<readlength;j++) fprintf(outfile,"%c",seq_array[pos+j]);
       fprintf(outfile,"\n");
    }
  }

  fprintf(stdout,"Done.\n"); 

  fclose(outfile);
  if(pe) fclose(outfile2);
  free(seq_array);
}


// returns 1 if valid, 0 if not.
char checkvalid(char* s,int n){
   int i;
   if(strlen(s)<n) return 0;
   for(i=0;i<n;i++){
      if(s[i]=='@' || s[i]=='$') return 0;
   }
   return 1;
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

