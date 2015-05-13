#include "emsar.h"
#include <getopt.h>

int main(int argc, char* argv[])
{
           char *segmentfile_name,*rshfile_name,*outdir,*outprefix; //sfafile_name is defined globally.
           char fastafilename[FILENAMEMAX],alnfilename[FILENAMEMAX],alnfilenamelistfile[FILENAMEMAX],input_rshfile_name[FILENAMEMAX]; 
           FILE *alnfilenamelistfileH;
           int nAlnfiles;
           char **alnfilenames;
           char *fraglengthfile_name;
           char *fpkmfile_name;
           char *posbiasfile_name;
           char mkdir_command[FILENAMEMAX+10];
           unsigned int i,j,j2,k,T;
           LIST *lst;
           int print_segments;
           int c,option_index;
           char fasta_option, heavy, bamflag, print_rsh_flag;
           int taglen;
           int round;
           char strand_type_str[5];
           char reiterate;
           char multisample;
           char mode; // for multi-sample, 0: first bam, 1: another bam (continued), 2: another bam (new)
           char line[FILENAMEMAX];

           if(argc<3) { printusage(argv[0]); return(0); }

           /* options */
           static struct option long_options[]= {
              {"rsh",no_argument, 0, 'I'}, /* user-specified rsh file as input */
              {"fasta",no_argument, 0, 'x'}, /* user-specified transcriptome fasta file as input */
              {"print_segments",no_argument, 0, 'g' }, /* generate an extra file that contains all the segment information */
              {"print_sfa",no_argument, 0, 'T'}, /* print suffix array */
              {"print_rsh",no_argument, 0, 'R'}, /* print rsh bucket structure */
              {"BAM",no_argument, 0, 'B' }, /* BAM input */
              {"SAM",no_argument, 0, 'S' }, /* SAM input */
              {"PE", no_argument, 0, 'P'}, /* paired-end */
              {"strand_type", required_argument, 0, 's'}, // SE: "ns","ssf","ssr", PE: "nsfr","nsff","ssfr","ssff","ssrr","ssrf"
              {"multisample", required_argument, 0, 'M'}, /* multisample */
              {"bias_model", required_argument, 0, 'm'}, /* model category for positional bias */
              {"posbias_training_len", required_argument, 0, 'W'}, /* length from 5' and 3' to train bias */
              {"posbias_impute_len", required_argument, 0, 'w'}, /* length toward 5' and 3' to use for imputing positional probability for untrained region */
              {"binsize", required_argument, 0, 'b' },
              {"maxthread", required_argument, 0, 'p' },
              {"header", required_argument,0,'h'}, /* fasta header option */
              {"taglen", required_argument,0,'t'}, /* tag length */
              {"maxfraglen",required_argument,0,'F' }, /* maximum fragment length */
              {"minfraglen",required_argument,0,'f' }, /* minimum fragment length */
              {"max_repeat",required_argument, 0, 'k'}, /* maximum number of locations a read is mapped to */
              {"nround",required_argument, 0, 'n'},
              {"epsilon",required_argument, 0, 'e'},
              {"precision",required_argument, 0, 'r'},  /* This is the same as CONVERGENCE_EPSILON_STEPSIZE */
              {"delta",required_argument, 0, 'd'},
              {"max_niter_mle",required_argument, 0, 'i'}, /* maximum number of iteration before determining not converging and reinitializing for MLE */
              {"verbose", no_argument, 0, 'v' },
              {"no_verbose", no_argument, 0, 'q' },
              {0,0,0,0}
           };

           /* default values */
           bin=5000;
           MAX_Thread=1;
           fasta_option='E';
           taglen=2;
           strcpy(strand_type_str,"ns"); /* because SE is defult */
           multisample=0;
           Max_Fraglength = 400; /* some safe number */
           Min_Fraglength = 1; /* some safe number */
           MAX_REPEAT=100;
           verbose_flag=1;
           print_segments=0;
           pe=0;
           posmodel=0;
           bamflag=0;
           perpos_freq_len=1000;
           perpos_freq_impute_len=200;
           print_sfa_flag=0;
           print_rsh_flag=0;
           strcpy(input_rshfile_name,""); /* no file specified by default */
           strcpy(fastafilename,""); /* no file specified by default */
 
           //parameters for MLE
           NUM_ROUND=4;
           CONVERGENCE_EPSILON = 1E-9;
           CONVERGENCE_EPSILON_STEPSIZE = 1E-15;  /* precision */
           DELTA=0; /* scaling factor for lambda for convenient calculation. */
           MAX_NITER_MLE = 200000;

           /* some initializations */
           nThread=0;
           EUMAcut=0;

           /* parsing options */
           while(1){
             option_index=0;
             c= getopt_long (argc, argv, "vqPs:b:p:h:t:F:f:n:e:r:p:d:gm:MHBSW:w:k:i:TRI:x:", long_options, &option_index);

             if(c==-1) break;
             switch(c){
                case 0:  /* ? */
                   if (long_options[option_index].flag!=0) break;
                   fprintf(stdout,"option %s", long_options[option_index].name);
                   if(optarg) fprintf(stdout," with arg %s",optarg);
                   fprintf(stdout,"\n");
                   break;
                case 'I':
                   strcpy(input_rshfile_name,optarg);
                   break;
                case 'x':
                   strcpy(fastafilename,optarg);
                   break;
                case 'P':
                   pe=1;
                   break;
                case 's':
                   strcpy(strand_type_str,optarg);
                   break;
                case 'b':
                   bin=atoi(optarg);
                   break;
                case 'p':
                   MAX_Thread= atoi(optarg);
                   if(MAX_Thread<1) { fprintf(stderr, "error: number of threads must be at least 1 (option -p).\n"); exit(1); }
                   if(MAX_Thread>1){
                     pth = (pthread_t*)malloc(sizeof(pthread_t)*(MAX_Thread-1));
                     working_thread = (int*)malloc(sizeof(int)*(MAX_Thread-1));
                   }
                   break;
                case 'h':
                   fasta_option=optarg[0];
                   break;
                case 't':
                   taglen=atoi(optarg);
                   if(taglen!=1 && taglen!=2 && taglen!=3) { fprintf(stderr,"error: currently taglength (-t) up to 3 is supported.\n"); exit(1); }
                   break;
                case 'F':
                   Max_Fraglength=atoi(optarg);
                   break;
                case 'f':
                   Min_Fraglength=atoi(optarg);
                   break;
                case 'k':
                   MAX_REPEAT=atoi(optarg);
                   break;
                case 'n':
                   NUM_ROUND = atoi(optarg);
                   if(NUM_ROUND<=0) { fprintf(stderr, "option -n must be a natural number.\n"); return(0); }
                   break;
                case 'e':
                   CONVERGENCE_EPSILON = atof(optarg);
                   if(CONVERGENCE_EPSILON<=0) { fprintf(stderr, "option -e must be positive.\n"); return(0); }
                   break;
                case 'r':
                   CONVERGENCE_EPSILON_STEPSIZE = atof(optarg);
                   if(CONVERGENCE_EPSILON_STEPSIZE <= 0 ) { fprintf(stderr, "option -p must be positive.\n"); return(0); }
                   break;
                case 'i':
                   MAX_NITER_MLE = atoi(optarg);
                   if(MAX_NITER_MLE <= 0 ) { fprintf(stderr, "option -i must be positive.\n"); return(0); }
                   break;
                case 'd':
                   DELTA = atoi(optarg);
                   break;
                case 'g':
                   print_segments=1;
                   break;
                case 'm':
                   posmodel= optarg[0];
                   if(posmodel=='0') posmodel=0; // no-bias
                   else if(posmodel=='1') posmodel=1; // bias
                   else { fprintf(stderr, "Invalid positional bias model (-m). Put either 0 or 1.\n"); return(0); }
                   break;
                case 'M':
                   multisample=1;
                   break;
                case 'H':
                   heavy=1;
                   break;
                case 'B':
                   if(bamflag=='s') { fprintf(stderr,"error: Options -B(--BAM) and -S(--SAM) cannot be used simultaneously.\n"); return(0); }
                   bamflag='b';
                   break;
                case 'S':
                   if(bamflag=='b') { fprintf(stderr,"error: Options -B(--BAM) and -S(--SAM) cannot be used simultaneously.\n"); return(0); }
                   bamflag='s';
                   break;
                case 'W':
                   perpos_freq_len = atoi(optarg);
                   if(perpos_freq_len<=0 || perpos_freq_len>=10000) { fprintf(stderr,"error: Option -W(--posbias_training_len) must be between 1 and 10000.\n"); return(0); }
                   break;
                case 'w':
                   perpos_freq_impute_len = atoi(optarg);
                   if(perpos_freq_impute_len<=0 || perpos_freq_impute_len>perpos_freq_len) { fprintf(stderr,"error: Option -w(--posbias_impute_len) must be between 1 and posbias_training_len.\n"); return(0); }
                   break;
                case 'T':
                   print_sfa_flag=1;
                   break;
                case 'R':
                   print_rsh_flag=1;
                   break;
                case 'v':
                   verbose_flag=2;
                   break;
                case 'q':
                   verbose_flag=0;
                   break;
                case '?':
                   fprintf(stderr,"error: unknown option?\n"); /*?*/
                default:
                   return(0);
             }
           }
           if(strlen(input_rshfile_name)==0 && strlen(fastafilename)==0) { fprintf(stderr,"error: either fasta file or an rsh file must be used as an input.\n"); exit(1); }
           if(Min_Fraglength>Max_Fraglength||Min_Fraglength<1||Max_Fraglength<1) { fprintf(stderr,"error: invalid fragment length range.\n"); exit(1); }
           if(set_library_strand_type(strand_type_str,pe)<0) { fprintf(stderr,"error: invalid strand type.\n"); exit(1); }

           if(verbose_flag>0) fprintf(stdout,"input fastafile name= %s\n",fastafilename);
           if(verbose_flag>0) fprintf(stdout,"input rshfile name= %s\n",input_rshfile_name);
           if(verbose_flag>0) fprintf(stdout,"Input type= %s\n",bamflag==0?"default bowtie output":(bamflag=='s'?"SAM":"BAM"));
           if(verbose_flag>0) fprintf(stdout,"Paired-end= %c\n",pe?'y':'n');
           if(verbose_flag>0) fprintf(stdout,"strand type= %s\n",strand_type_str);
           if(verbose_flag>0) fprintf(stdout,"Multisample= %c\n",multisample?'y':'n');
           if(verbose_flag>0) fprintf(stdout,"Max_Fraglen= %d\n",Max_Fraglength);
           if(verbose_flag>0) fprintf(stdout,"Min_Fraglen= %d\n",Min_Fraglength);
           if(verbose_flag>0) fprintf(stdout,"MAX_REPEAT= %d\n",MAX_REPEAT);
           if(verbose_flag>0) fprintf(stdout,"bias model= %d %s\n",posmodel,posmodel==0?"(no bias model)":"");
           if(verbose_flag>0) fprintf(stdout,"positional bias training length= %d\n",perpos_freq_len);
           if(verbose_flag>0) fprintf(stdout,"positional bias impute training length= %d\n",perpos_freq_impute_len);
           if(verbose_flag>0) fprintf(stdout,"fasta header option= %c\n",fasta_option);
           if(verbose_flag>0) fprintf(stdout,"MAX_Thread= %d\n",MAX_Thread);
           if(verbose_flag>0) fprintf(stdout,"NUM_ROUND= %d\n",NUM_ROUND);
           if(verbose_flag>0) fprintf(stdout,"CONVERGENCE_EPSILON= %g\n",CONVERGENCE_EPSILON);
           if(verbose_flag>0) fprintf(stdout,"CONVERGENCE_EPSILON_STEPSIZE= %g\n",CONVERGENCE_EPSILON_STEPSIZE);
           if(verbose_flag>0) fprintf(stdout,"MAX_NITER_MLE= %d\n",MAX_NITER_MLE);
           if(verbose_flag>0) fprintf(stdout,"binsize = %d\n",bin);
           if(verbose_flag>0) fprintf(stdout,"taglen = %d\n",taglen);
           if(verbose_flag>0) fprintf(stdout,"print segments = %c\n",print_segments?'y':'n');
           if(verbose_flag>0) fprintf(stdout,"print suffix aray = %c\n",print_sfa_flag==1?'y':'n');
           if(verbose_flag>0) fprintf(stdout,"print rsh structure = %c\n",print_rsh_flag==1?'y':'n');

          

           /* non-option arguments */
           if(optind+1<argc) {
              outdir = (char*)argv[optind];
              outprefix = (char*)argv[optind+1];
              if(optind+2<argc) strcpy(alnfilenamelistfile,argv[optind+2]); else strcpy(alnfilenamelistfile,"");

              /* alignment file names */
              nAlnfiles=0;
              alnfilenames = (char**)malloc(MAX_nALNFILES*sizeof(char*));
              if( multisample==0 ) { 
                alnfilenames[0] = (char*)malloc(FILENAMEMAX*sizeof(char));
                strcpy(alnfilenames[0],alnfilenamelistfile); nAlnfiles=1; 
              }
              else {
                if( (alnfilenamelistfileH = fopen(alnfilenamelistfile,"r"))!=0 ){
                  while(fgets(line,sizeof(line),alnfilenamelistfileH)){
                    if(line[strlen(line)-1]=='\n') line[strlen(line)-1]='\0'; // chomp
                    alnfilenames[nAlnfiles] = (char*)malloc(FILENAMEMAX*sizeof(char));
                    strcpy(alnfilenames[nAlnfiles],line); nAlnfiles++;
                  }
                  fclose(alnfilenamelistfileH);
                }
                else { fprintf(stderr,"Can't open alignment list file.\n"); return(1); }
              }
              if(nAlnfiles==0) { fprintf(stderr,"No alignment files in the alignment list\n"); return(1); }
              if(nAlnfiles<MAX_nALNFILES) alnfilenames = realloc(alnfilenames,nAlnfiles*sizeof(char*));
              if(alnfilenames==NULL) { fprintf(stderr,"failed in reallocation of alnfilenames.\n"); return(1); } 
           }
           else{ printusage(argv[0]); return(0); }

           if(verbose_flag>0) fprintf(stdout,"finished reading options and arguments..\n");

           sprintf(mkdir_command,"mkdir -p %s",outdir);
           fflush(stdout); system(mkdir_command);

           /* rshbucket pointer setup (in case an alternative function is used later with added features */
           update_rshbucket_PTR = &update_rshbucket;
           update_rshbucket_single_PTR = &update_rshbucket_single; 
           clear_readcounts_in_rshbucket_PTR = &clear_readcounts_in_rshbucket;

           //output sfa filename
           sfafile_name = (char*)malloc((strlen(outdir)+strlen(outprefix)+6)*sizeof(char));
           sprintf(sfafile_name,"%s/%s.sfa",outdir,outprefix);

           //output rsh filename
           rshfile_name = (char*)malloc((strlen(outdir)+strlen(outprefix)+6)*sizeof(char));
           sprintf(rshfile_name,"%s/%s.rsh",outdir,outprefix);

           if(strlen(input_rshfile_name)==0){ // no rshfile specified. learn it from fasta file.
  
             /* read fasta file and store it to fasta array */
             if(verbose_flag>0) { fprintf(stdout, "reading fasta file... :"); fflush(stdout); system("date +%m/%d,%T"); }
             read_raw_fasta(fastafilename,fasta_option);
  

             /* determine read length */
             // for pe, read just the first line from bam file. We assume the read length is fixed.
             // for se, read the whole bam file to get the read length range as Fraglength.min - Fraglengths.max.
             if(pe){
               if(bamflag!=0) read_BAM_get_readlength (alnfilenames[0], bamflag);
               else read_bowtie_get_readlength (alnfilenames[0]);
             } else {
               if(bamflag!=0) read_BAM_get_readlengths_se (alnfilenames[0], bamflag);
               else read_bowtie_get_readlengths_se (alnfilenames[0]);
             }  


             /* determine fragment length */
             if(pe==1) { 
                if(verbose_flag>0) { fprintf(stdout, "reading fragment length range...:"); fflush(stdout); system("date +%m/%d,%T"); }
                determine_fraglength_range(); 
                if(verbose_flag>0) { fprintf(stdout, "fragment length range : %d - %d\n",Fraglengths.min, Fraglengths.max); }
                if(verbose_flag>0) { fprintf(stdout, "read length : %d\n",readlength); }
             } else {
                if(verbose_flag>0) { fprintf(stdout, "read length range : %d - %d\n",Fraglengths.min, Fraglengths.max); }
             }
  
             /* initialize rshbucket */ 
             if(verbose_flag>0) { fprintf(stdout, "initializing rsh array... :");fflush(stdout); system("date +%m/%d,%T"); }
             initialize_rshbucket(-1);
  
             //dinucleotide tag filter
             local_sfa_end=-1; local_sfa_start=0; //initialize
      
             //initializing, sorting suffix array and constructing rshbucket
             if(pe==0) preprocess_SE(taglen);
             else preprocess_PE(taglen);
             
             /* freeing some memory */
             free(cuml); cuml=NULL;
             free(cumlI); cumlI=NULL;
             if(pe){
               free(sfa_m);
             }

           } else {
             if(verbose_flag>0) { fprintf(stdout, "reading rsh array... :");fflush(stdout); system("date +%m/%d,%T"); }
             construct_rsh_from_rshfile(input_rshfile_name); 
           }

           

           /* Initialize posbias arrays */
           perpos_freq_5 = (double*)calloc(perpos_freq_len,sizeof(double));
           perpos_freq_3 = (double*)calloc(perpos_freq_len,sizeof(double));
           perpos_unavail_freq_5 = (double*)calloc(perpos_freq_len,sizeof(double));
           perpos_unavail_freq_3 = (double*)calloc(perpos_freq_len,sizeof(double));
           perpos_normfreq_5 = (double*)calloc(perpos_freq_len,sizeof(double));
           perpos_normfreq_3 = (double*)calloc(perpos_freq_len,sizeof(double));
           

           /* normalize positional probability */
           if(posmodel==1){
            normalize_perpos_freq();
            posbiasfile_name = (char*)malloc((strlen(outdir)+strlen(outprefix)+10)*sizeof(char));
            sprintf(posbiasfile_name,"%s/%s.posbias",outdir,outprefix);
            print_posbias(posbiasfile_name);
            //return(0);  // DEBUGGING

            //determine_scaling_factor_for_perpos_prob();
           }




           /* read alignment file  */
           if(verbose_flag>0) { fprintf(stdout, "reading alignment file(s)... :"); fflush(stdout); system("date +%m/%d,%T"); }
           mode=2;
           for(i=0;i<nAlnfiles;i++){

             /* set read counts to be zero */
             FraglengthCounts=(int*)calloc(Max_Fraglength+1,sizeof(int));
             (*clear_readcounts_in_rshbucket_PTR)();
             fprintf(stdout,"alnfile[%d]=%s\n",i,alnfilenames[i]); /* DEBUGGING */
             if(bamflag==0) { // default bowtieout
                if(pe) read_bowtie_PE(alnfilenames[i], mode);  /* add '-' option later to SE main. The function can handle it. */ /* differnet option is needed for SS PE (later) */
                else read_bowtie_SE(alnfilenames[i], mode);  /* add '-' option later to SE main. The function can handle it. */ /* differnet option is needed for SS PE (later) */
             }
             else {
                if(pe) read_BAM_PE(alnfilenames[i], bamflag, mode);
                else read_BAM_SE(alnfilenames[i], bamflag, mode);
             }


             transfer_fraglendist_to_Wf();  // compute_Wf


             /* print rsh here if requested */
             if(print_rsh_flag==1) print_rsh(rshfile_name);

             /* rsh is deleted during this process */
             if(verbose_flag>0) { fprintf(stdout, "\nscanning rsh array and constructing EUMA, ReadCount and CT array... :");fflush(stdout); system("date +%m/%d,%T");}
             scan_rshbucket();
  
  
             /* dividing transcripts into disjoint sets */
             if(verbose_flag>0) { fprintf(stdout, "constructing transcript-combination array... :");fflush(stdout); system("date +%m/%d,%T"); }
             build_TC_from_CT_2();
             if(verbose_flag>0) { fprintf(stdout, "constructing disjoint sets... :");fflush(stdout); system("date +%m/%d,%T"); }
             do {
               initialize_CS_TS();
               T=0; reiterate=0; 
               for(j=0;j<=max_cid;j++) {
                 ntid_per_sid=0;
                 if(propagate_2(j,T)) T++;
                 if(ntid_per_sid>MAX_NTID_PER_SID) {
                   EUMAcut+=EUMACUT_INCREMENT;
                   delete_CS_TS();
                   reiterate=1;
                   if(verbose_flag>0) { fprintf(stdout, "module size too big (%d). EUMAcut is readjusted to %.0f\n",ntid_per_sid,EUMAcut);fflush(stdout); system("date +%m/%d,%T"); }
                   break; // break from for loop
                 }
               }
             }while(reiterate==1);
  
             if(verbose_flag>0) { fprintf(stdout, "finished constructing disjoint sets... :");fflush(stdout); system("date +%m/%d,%T"); }
  
             max_sid=T-1;
             generate_SC_ST();
             if(verbose_flag>0) { fprintf(stdout, "finished constructing SC and ST arrays... :");fflush(stdout); system("date +%m/%d,%T"); }
  
  
             /* preparation */
             if(verbose_flag>0) fprintf(stdout, "constructing EUMAps array...\n");
             construct_EUMAps();
             //for(i=0;i<=max_cid;i++) fprintf(stderr,"EUMAps[%d]=%lf\n",i,EUMAps[i]); //DEBUGGING
             FPKMfinal = (double**)malloc(NUM_ROUND*sizeof(double*));
             FPKM=(double*)malloc((max_tid+1)*sizeof(double));
             if(FPKM==NULL) { fprintf(stderr,"Failed to allocate memory to FPKM\n"); exit(1); }
             srand(time(NULL));
           
             /* running MLE */
             for(round=0;round<NUM_ROUND;round++){
                if(verbose_flag>0) { fprintf(stdout,"round %d/%d...\n",round+1,NUM_ROUND); fflush(stdout); system("date +%m/%d,%T"); } /* progress log */
                run_MLE_threads();
                if(verbose_flag>0) fprintf(stdout,"MLE finished...\n"); /* progress log */
                construct_FPKMfinal(round);
                if(verbose_flag>0) fprintf(stdout,"FPKMfinal finished...\n"); /* progress log */
             }
  
             // computing iEUMA
             if(verbose_flag>0) {fprintf(stdout, "computing effective length for read count inference...\n");fflush(stdout); system("date +%m/%d,%T");} 
             compute_iEUMA();
           
             if(verbose_flag>0) fprintf(stdout, "printing FPKM file...\n");
             //output fpkm file name
             fpkmfile_name = (char*)malloc((strlen(outdir)+strlen(outprefix)+7+5)*sizeof(char));
             sprintf(fpkmfile_name,"%s/%s.%d.fpkm",outdir,outprefix,i);
             print_FPKMfinal(fpkmfile_name);

             //output fraglengtheffect file name
             fraglengthfile_name = (char*)malloc((strlen(outdir)+strlen(outprefix)+20+5)*sizeof(char));
             sprintf(fraglengthfile_name,"%s/%s.%d.fraglength_effect",outdir,outprefix,i);
             print_FraglengthDist(fraglengthfile_name);

             //print segment file  
             segmentfile_name=(char*)malloc((strlen(outdir)+strlen(outprefix)+11+5)*sizeof(char));
	     sprintf(segmentfile_name,"%s/%s.%d.segments",outdir,outprefix,i);
             if(print_segments) print_aEUMA_3(segmentfile_name,T-1);

             //finishing log.
             fprintf(stdout,"Complete: Output file :\n  %s\n  %s\n",fpkmfile_name,fraglengthfile_name);
             if(print_segments) fprintf(stdout,"  %s\n",segmentfile_name);
             fflush(stdout); system("date +%m/%d,%T");
            
             //freeing some space
             delete_CS_TS();
             free(EUMAps); EUMAps=NULL;
             free(FPKM); FPKM=NULL;
             delete_FPKMfinal();
             free(iEUMA); iEUMA=NULL;
             free(ReadCount); ReadCount=NULL;
             delete_CT_SC_ST_EUMA();
             free(FraglengthCounts);
             free(Wf);

           }

         
           // freeing memory 

           //allocated memory for threading.
           if(MAX_Thread>1){
             free(pth); pth=NULL;
             free(working_thread); working_thread=NULL;
           }

           if(verbose_flag>0) { fprintf(stdout, "freeing rsh array ... :");fflush(stdout); system("date +%m/%d,%T"); }
           delete_rshbucket();

           delete_index_table();
           delete_treenode(tname_tree); tname_tree=NULL;

           free(segmentfile_name); segmentfile_name=NULL;
           free(fraglengthfile_name); fraglengthfile_name=NULL;
           free(fpkmfile_name); fpkmfile_name=NULL;


           return 0;
}



