#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <limits.h>
#include <math.h>
#include "stringhash.h"
#include "alignment.h"
#include "sam.h"

#define SFA_LINE_MAX 200 /* maximum length of a suffix array line */
#define SINGLE_LINE_MAX 1000   /* maximum length of a bowtie output line */
#define RSH_SINGLE_LINE_MAX 100000 /* maximum length of a rsh main line */
#define ISOFORM_MAX 500000   /* maximum number of isoforms */
#define FILENAMEMAX 1000  /* maximum length of output file path */
#define INIT_RSHBUCKET_MAX_T_SIZE 10 /* initial maximum rshbucket max t_size. (actual size of the rshbucket is this number - 1 */
#define MAX_NTID_PER_SID 5000 /* maximum size of sequence-sharing module (set) */
#define EUMACUT_INCREMENT 2

#define LOWEST_NUMBER -1E308 /* for MLE */
#define NEAR_LOWEST_NUMBER -9.9E307 /* for MLE */

#define MAX_NLOOP_MLE 10000 /* maximum number of loops of iterations before giving up for MLE */
#define MAX_nALNFILES 1000 /* maximum number of bam files to process in a single multisample run */


#define flip_bool(x) x?0:1; /* 0->1, 1->0 */


/* The following have been modified:
variables: node1, EUMA_array
function declaration: update_rshbucket, newNode1, scan_rshbucket */


typedef struct DIST_RANGE
{
     int min;
     int max;
} D_RANGE;

typedef struct _SFA_RANGE
{
     int left;
     int right;
} SFA_RANGE;

typedef struct INTVARPAIR
{
     int start;
     int end;
} INTPAIR;

typedef struct intarray {
   int size;
   int* array;
} inta;

typedef struct intarray2 {
   int size;
   int* array;
   char* array2;
} inta2;

struct iLIST
{
    unsigned int val;
    struct iLIST* next;
};
typedef struct iLIST LIST;

typedef struct TCelement
{
     LIST* list;
     LIST* last;
} TCE;

typedef struct rshbucket_node
{
    struct rshbucket_node *next;
    int *EUMA;  // array of int values instead of a single double value
    int ReadCount; // total read count for each segment
    //char *poscat; //positional category
    int tarr[1];  /* variable size */
} node1;






/*-- global_variables --*/
// global parameters
char pe;
char library_strand_type;
double EUMAcut;
int ntid_per_sid; 
int MAX_REPEAT; /* maximum number of occurrences of the same sequence of length readlength that is counted into EUMA */
int MAX_NITER_MLE; /* maximum number of iteration before determining not converging and reinitializing for MLE */
char print_sfa_flag; /* whether to print out suffix array */
char* sfafile_name;

//sequence and sequence index
char *seq_array;
int *cuml;
INTPAIR *cumlI; /* index for cuml */
int bin; /* bin for cumlI */
int borderpos; /* the position of the border '$'. = the length of concatenated fw sequence. */
int seqlength; /* length of concatenated sequence. = the position of the last '$'. */
int fasta_size; /* size of the seq_array */
char* noncanonarr; /* marking noncanonical substrings. bitstring. */
char remove_polyA;  /* 1: remove polyA from fasta, 0: don't remove polyA */

//transcript name indexing
unsigned int max_tid;
treenode* tname_tree;
char **IndexTable;

// suffix array and suffix array index
int *sfa;  
char *sfa_m; /* suffix array for frp class */
char *sfa_m_pe, *sfa_f;
int sfa_size;
int max_sfa_size;  /* when sfa is constructed by tags, this is the max sfa size. */
int local_sfa_start; /* for initializing and sorting sfa (esp for SE) */
int local_sfa_end; /* for initializing and sorting sfa (esp for SE) */


//tid-cid-sid mapping
inta2 *CT;  /* CT (combination-transcript) array */
TCE *TC;  /* TC (transcript-combination) array */
int *CS,*TS; /* combination-set array & transcript-set array */
inta *SC;
inta *ST;
int CT_size;
int max_sid;
int max_cid;

//rshbucket, EUMA & read count
node1 ***rshbucket;
int rshbucket_max_t_size;
node1 **rshbucket_single;
int *ReadCount;
double *adjEUMA;
int rsh_size;
int rsh_size_single;


//readlengths & Fragment lengths
int Max_Fraglength;
int Min_Fraglength;
int *FraglengthCounts;
D_RANGE Fraglengths, Readlengths;
int readlength;
double *Wf;  // fragment length sampling probability
int nFraglen;

//positional bias
char posmodel; // 0 : no bias model. 'a' : arbitrary (primitive version for testing & debugging)

//MLE & FPKM
double *FPKM;
double **FPKMfinal;
double *EUMAps;
int TotalReadCount;
double CONVERGENCE_EPSILON;  /* precision of ML function */
double CONVERGENCE_EPSILON_STEPSIZE;  /* precision of the estimate */
int NUM_ROUND;
double LOGCONST;
double DELTA;
double *iEUMA; /* total effective transcript length for each tid */


//positional bias
double *perpos_freq_5;
double *perpos_freq_3;
double *perpos_unavail_freq_5;
double *perpos_unavail_freq_3;
double *perpos_normfreq_5;
double *perpos_normfreq_3;
double *Y; // scale factor for perpos probability for individual transcript
int perpos_freq_len;
int perpos_freq_impute_len;


//threading
int MAX_Thread,nThread;
pthread_t *pth;
int *working_thread;
int nDone; /* number of finished tid's (for progress bar for MLE incorporating threading) or suffix array scanning (for progress bar for PE suffix array scanning) */
double pDone; /* percentage of finished tid's (for progress bar for MLE incorporating threading) or suffix array scanning (for progress bar for PE suffix array scanning) */
int pDone_int; /* integer part of pDone */

// others
int verbose_flag; /* 0 : not verbose, 1: default, 2: very verbose */
static char base[4]={'A','C','G','T'};



/* function declaration */
char set_library_strand_type(char*,char);

//rshbucket & fraglength (light version)
node1* newNode1(int*,int,node1*,char,int,char*);
node1* newNode2 (int *, int, node1*, int*, char*);
short int cmptarr(int*,int*,int);
char update_rshbucket(int,int*,char,int,char*);
char update_rshbucket_single(int, char, int, char*);
void initialize_rshbucket(int);
void delete_rshbucket(void);
void scan_rshbucket(void);
void transfer_fraglendist_to_Wf (void);
void clear_readcounts_in_rshbucket(void);
void parse_rsh_indexline (char*);
void parse_rsh_headerline (char*);
node1* parse_rsh_mainline (char*,node1*);
void construct_rsh_from_rshfile(char*);


char (*update_rshbucket_PTR)(int,int*,char,int,char*);  // may refer to either update_rshbucket or update_rshbucket_heavy
char (*update_rshbucket_single_PTR)(int, char, int, char*); // may refer to either update_rshbucket_single or update_rshbucket_single_heavy
void (*clear_readcounts_in_rshbucket_PTR)(void);

double compute_adjEUMA(int*);
void renormalize_Wf(void);

/* common */
void initialize_tname_indextable(void);
void read_raw_fasta(char*,char);
void determine_fraglength_range(void);
void determine_sfa_size(void);
void update_ReadCounts(alignment_list*);
void print_FraglengthDist(char*);
void build_TC_from_CT_2(void);
void delete_CS_TS(void);
void initialize_CS_TS(void);
int propagate_2(unsigned int,unsigned int);
void print_aEUMA_3(char*,int);
int sf_i(int);
int sf_ib(int);
int sf_p(int,int);
int sf_pb(int,int);
int flip(int);
void *quick_sort_suffixarray_4(void *arg);
int partition_suffixarray_3(int, int, int);
char revcomp(char);  // A->T
char uc(char); // a->A
void parse_ensembl_header(char*,char**);
void parse_refseq_header(char*,char**);
int generate_seqtag(char***, int);
void delete_seqtag(char***, int);
void printusage(char*);
void generate_SC_ST(void);
void read_rsh(char*);


//positional bias
char determine_poscategory(int, int, int);
void normalize_perpos_freq(void);
void determine_scaling_factor_for_perpos_prob(void);
void print_posbias(char*);


/* se-specific */
void preprocess_SE(int);
void read_bowtie_SE (char*,char);
alignment* parse_bowtieline(char*, char**);
void construct_rshbucket_2(int fraglength);
void construct_rshbucket_2_posbias(int fraglength);
void initialize_suffixarray_NS_5(char*);
void initialize_suffixarray_SS_4(char*);


/* pe-specific */
void preprocess_PE(int);
void read_bowtie_PE (char*, char);
alignment* parse_bowtieline_PE (char*, char*, char**);
char check_mate_strand(char,char);
void initialize_suffixarray_NS_PE_2(char*);
void initialize_suffixarray_SS_PE_2(char*);
void run_process_mate1_cluster_by_mate_2_PE_threads_2(void);
void* process_mate1_cluster_by_mate_3(void*);
void quick_sort_suffixarray_by_mate_2(int**, int**, int, int);
int partition_suffix_array_by_mate_2(int**, int**, int, int, int);
void construct_rshbucket_PE_3(int**, int**, int);
void construct_rshbucket_PE_3_posbias(int**, int**, int);
void mark_sfa(void);
void mark_sfa_pe(int);
void mark_noncanonical(void);
char is_noncanonical(char*,int);
int strcmp_pe(char*,char*,int);
int strcmp_se(char*,char*);


/* MLE & FPKM */
int MLE (int);
double Fp (inta);
double lambdap (int);
void construct_FPKMfinal (int);
void construct_EUMAps (void);
void print_FPKMfinal (char*);
void print_SC(void);
void print_CT(void);
void print_ST(void);
void run_MLE_threads(void);
void* MLE_range(void*);
void delete_CT_SC_ST_EUMA(void);
void delete_FPKMfinal(void);
void compute_iEUMA(void);
void delete_index_table (void);
int Round_off(double);
void MLE_progressbar(void);
void print_sfa (void);
void print_rsh (char*);


static char NUMBERS_STR[10]="0123456789";

void read_BAM_get_readlength (char*, char);
void read_bowtie_get_readlength (char*);
int parse_bowtieline_get_readlength (char*);
void read_BAM_SE (char*, char, char);
alignment* convert_bam_alignment_2_alignment (bam1_t*, char*);
void read_BAM_PE (char*, char, char);
alignment* convert_bam_alignment_2_alignment_PE (bam1_t*,bam1_t*, char*);
int parse_SAM_mmstr(char*);


