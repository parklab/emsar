#include "emsar.h"
#include <getopt.h>

int main(int argc, char* argv[])
{
           char *fastafilename,*readlength_str,*outdir,*outprefix,*rshfile_name; //sfafile_name is defined globally.
           char mkdir_command[FILENAMEMAX+10];
           unsigned int i,j,k;
           int c,option_index;
           char fasta_option, bamflag, print_rsh_flag;
           int taglen;
           char strand_type_str[5];
           char line[FILENAMEMAX];

           if(argc<4) { printusage_build(argv[0]); return(0); }

           /* options */
           static struct option long_options[]= {
              {"print_sfa",no_argument, 0, 'T'}, /* print suffix array */
              {"PE", no_argument, 0, 'P'}, /* paired-end */
              {"strand_type", required_argument, 0, 's'}, // SE: "ns","ssf","ssr", PE: "nsfr","nsff","ssfr","ssff","ssrr","ssrf"
              {"bias_model", required_argument, 0, 'm'}, /* model category for positional bias */
              {"posbias_training_len", required_argument, 0, 'W'}, /* length from 5' and 3' to train bias */
              {"posbias_impute_len", required_argument, 0, 'w'}, /* length toward 5' and 3' to use for imputing positional probability for untrained region */
              {"binsize", required_argument, 0, 'b' },
              {"maxthread", required_argument, 0, 'p' },
              {"max_repeat",required_argument, 0, 'k'}, /* maximum number of locations a read is mapped to */
              {"header", required_argument,0,'h'}, /* fasta header option */
              {"taglen", required_argument,0,'t'}, /* tag length */
              {"maxfraglen",required_argument,0,'F' }, /* maximum fragment length */
              {"minfraglen",required_argument,0,'f' }, /* minimum fragment length */
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
           Max_Fraglength = 400; /* some safe number */
           Min_Fraglength = 1; /* some safe number */
           MAX_REPEAT=100;
           verbose_flag=1;
           pe=0;
           posmodel=0;
           perpos_freq_len=1000;
           perpos_freq_impute_len=200;
           print_sfa_flag=0;
           nThread=0;

           /* parsing options */
           while(1){
             option_index=0;
             c= getopt_long (argc, argv, "vqPs:b:p:h:t:F:f:p:m:W:w:Tk:", long_options, &option_index);

             if(c==-1) break;
             switch(c){
                case 0:  /* ? */
                   if (long_options[option_index].flag!=0) break;
                   fprintf(stdout,"option %s", long_options[option_index].name);
                   if(optarg) fprintf(stdout," with arg %s",optarg);
                   fprintf(stdout,"\n");
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
                case 'k':
                   MAX_REPEAT=atoi(optarg);
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
                case 'm':
                   posmodel= optarg[0];
                   if(posmodel=='0') posmodel=0; // no-bias
                   else if(posmodel=='1') posmodel=1; // bias
                   else { fprintf(stderr, "Invalid positional bias model (-m). Put either 0 or 1.\n"); return(0); }
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
           if(Min_Fraglength>Max_Fraglength||Min_Fraglength<1||Max_Fraglength<1) { fprintf(stderr,"error: invalid fragment length range.\n"); exit(1); }
           if(set_library_strand_type(strand_type_str,pe)<0) { fprintf(stderr,"error: invalid strand type.\n"); exit(1); }

           if(verbose_flag>0) fprintf(stdout,"Paired-end= %c\n",pe?'y':'n');
           if(verbose_flag>0) fprintf(stdout,"strand type= %s\n",strand_type_str);
           if(verbose_flag>0) fprintf(stdout,"Max_Fraglen= %d\n",Max_Fraglength);
           if(verbose_flag>0) fprintf(stdout,"Min_Fraglen= %d\n",Min_Fraglength);
           if(verbose_flag>0) fprintf(stdout,"MAX_REPEAT= %d\n",MAX_REPEAT);
           if(verbose_flag>0) fprintf(stdout,"bias model= %d %s\n",posmodel,posmodel==0?"(no bias model)":"");
           if(verbose_flag>0) fprintf(stdout,"positional bias training length= %d\n",perpos_freq_len);
           if(verbose_flag>0) fprintf(stdout,"positional bias impute training length= %d\n",perpos_freq_impute_len);
           if(verbose_flag>0) fprintf(stdout,"fasta header option= %c\n",fasta_option);
           if(verbose_flag>0) fprintf(stdout,"MAX_Thread= %d\n",MAX_Thread);
           if(verbose_flag>0) fprintf(stdout,"binsize = %d\n",bin);
           if(verbose_flag>0) fprintf(stdout,"taglen = %d\n",taglen);
           if(verbose_flag>0) fprintf(stdout,"print suffix aray = %c\n",print_sfa_flag==1?'y':'n');

          

           /* non-option arguments */
           if(optind+3<argc) {
              fastafilename = (char*)argv[optind];
              readlength_str = (char*)argv[optind+1];
              if(pe==1) readlength = atoi(readlength_str);
              else parse_readlength_range(readlength_str);  /* parse it as a range */
              outdir = (char*)argv[optind+2];
              outprefix = (char*)argv[optind+3];
           }
           else{ printusage_build(argv[0]); return(0); }

           if(verbose_flag>0) fprintf(stdout,"finished reading options and arguments..\n");
           if(verbose_flag>0) fprintf(stdout,"input fastafile name= %s\n",fastafilename);
           if(verbose_flag>0) if(pe==1) fprintf(stdout,"readlength = %d\n",readlength);
           if(verbose_flag>0) if(pe==0) fprintf(stdout,"readlength range = %d - %d\n",Readlengths.min,Readlengths.max);

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

  
           //calculate nFraglen for PE
           if(pe==1) determine_fraglength_range();


           /* read fasta file and store it to fasta array */
           if(verbose_flag>0) { fprintf(stdout, "reading fasta file... :"); fflush(stdout); system("date +%m/%d,%T"); }
           read_raw_fasta(fastafilename,fasta_option);


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


           /* print rsh */
           if(verbose_flag>0) { fprintf(stdout, "writing to an rsh file ... :");fflush(stdout); system("date +%m/%d,%T"); }
           print_rsh(rshfile_name);


           /* freeing / deleting stuff */
           if(verbose_flag>0) { fprintf(stdout, "freeing memory ... :");fflush(stdout); system("date +%m/%d,%T"); }
           if(MAX_Thread>1){
             free(pth); pth=NULL;
             free(working_thread); working_thread=NULL;
           }

           delete_rshbucket();
           delete_index_table();
           delete_treenode(tname_tree); tname_tree=NULL;
           free(sfafile_name);
           free(rshfile_name);

           return 0;
}



