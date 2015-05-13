#include "readgenerator.h"
#include <getopt.h>
//This is a minimal RNA-seq simulator just for a simple test-run of EMSAR. This script is not meant to be used for evaluation of RNA-seq analysis methods.


int main(int argc, char* argv[])
{
           char *fastafilename,*outdir,*outfile,*outfilepath,*outfilepath2;
           char mkdir_command[FILENAMEMAX+10];
           int strand_specific_option;
           int readlength, numreads, fraglen;
           int c,option_index;
           char header_prefix[MAX_HEADER_PREFIX_LEN]="";
           char pe;

           if(argc<5) { printf("Usage: %s <options> fastafile readlength numreads outdir outfilename\n",argv[0]); return(0); }

           /* options */
           static struct option long_options[]= {
              {"PE", no_argument, 0, 'P' },
              {"ss", no_argument, 0, 's' },
              {"f", required_argument, 0, 'f' },
              {"header_prefix",required_argument, 0, 'h' },
              {0,0,0,0}
           };

           /* default values */
           pe=0; // default, SE
           strand_specific_option=0;
           fraglen=0;

           /* parsing options */
           while(1){
             option_index=0;
             c= getopt_long (argc, argv, "Psf:h:", long_options, &option_index);

             if(c==-1) break;
             switch(c){
                case 'P':
                   pe=1;
                   break;
                case 's':
                   strand_specific_option=1;
                   break;
                case 'f':
                   fraglen=atoi(optarg);
                   break;
                case 'h':
                   strcpy(header_prefix,optarg);
                   break;
                case '?':
                   fprintf(stderr,"error: unknown option?\n"); /*?*/
                default:
                   return(0);
             }
           }

           if(pe && fraglen==0) { fprintf(stderr,"fraglen must be provided for PE.(-f)\n"); exit(1); }

           /* non-option arguments */
           if(optind+4<argc) {
              fastafilename = (char*)argv[optind];
              readlength = atoi(argv[optind+1]);
              numreads = atoi(argv[optind+2]);
              outdir = (char*)argv[optind+3];
              outfile = (char*)argv[optind+4];
           }
           else{ printf("Usage: %s <options> fastafile readlength numreads outdir outfilename\n",argv[0]); return(0); }


           sprintf(mkdir_command,"mkdir -p %s",outdir);
           fflush(stdout); system(mkdir_command);

           //output filename
           if(pe){
             outfilepath=(char*)malloc((strlen(outdir)+strlen(outfile)+5)*sizeof(char));
	     sprintf(outfilepath,"%s/%s.R1",outdir,outfile);
             outfilepath2=(char*)malloc((strlen(outdir)+strlen(outfile)+5)*sizeof(char));
	     sprintf(outfilepath2,"%s/%s.R2",outdir,outfile);
           }
           else {
             outfilepath=(char*)malloc((strlen(outdir)+strlen(outfile)+2)*sizeof(char));
	     sprintf(outfilepath,"%s/%s",outdir,outfile);
             outfilepath2=NULL;
           }

           generate_reads(fastafilename, readlength, fraglen, numreads, strand_specific_option, header_prefix, outfilepath, outfilepath2, pe);
           free(outfilepath); outfilepath=NULL;
           if(pe) { free(outfilepath2); outfilepath2=NULL; }

           return 0;
}



