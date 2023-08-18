#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "fasta-genome-io.h"
#include "sam-parse.h"
#include "kmer.h"

#define TRUE (1)
#define FALSE (0)
#define DEBUG (0)

static int KLEN = 8;
static int MIN_MQ = 0;
static unsigned long MIN_READ_LEN = 0;
static unsigned long MAX_READ_LEN = 250000000;
static int OPT_M = FALSE;
static int OPT_T = FALSE;


/* Get the reverse complement of a sequence
Args: const char* seq - input sequence that we want to get the
                           the reverse complement of
      char* rvcmp - pointer to char array where the reverse
                         complemented sequence is stored
      size_t seq_len - length of the input sequence */ 
void do_rvcmp(const char* seq, char* rvcmp, size_t seq_len) {
    for (int i = seq_len-1; i >= 0; i--) {
        char base = seq[i];
        if ( (base == 'A') || (base == 'a') ) {
            rvcmp[seq_len-(i+1)] = 'T';
        }
        else if ( (base == 'C') || (base == 'c') ) {
            rvcmp[seq_len-(i+1)] = 'G';
        }
        else if ( (base == 'G') || (base == 'g') ) {
            rvcmp[seq_len-(i+1)] = 'C';
        }
        else if ( (base == 'T') || (base == 't') ) {
            rvcmp[seq_len-(i+1)] = 'A';
        }
        else {
            rvcmp[seq_len-(i+1)] = base;
        }
    }
}


/* Check if read length is between values specified by -l and -L
   Args: const int seq_len - length of the read
   Returns: 1 if passes; 0 otherwise */
int read_len_ok(const int seq_len) {
    if ( (seq_len >= MIN_READ_LEN) && 
         (seq_len <= MAX_READ_LEN) ) {
        return 1;
    }
    return 0;
}


/* Check if CIGAR string only contains a number (that equals 
   the aligned read length) and an M 
   Args: const int seq_len - length of read whose CIGAR we want to check 
         const char* cigar - input CIGAR string
   Returns: 1 if cigar only consists of a number equal to seq_len
            followed by the character 'M' (no other characters allowed);
            0 otherwise */
int cigar_ok(const int seq_len, const char* cigar) {
    int n_digits = snprintf(NULL, 0, "%d", seq_len); // get number of digits in seq_len
    char buf[n_digits + 2]; // extra space for 'M' and NULL
    snprintf(buf, n_digits + 2, "%d", seq_len); // store seq_len as string in buf
    buf[n_digits] = 'M';
    buf[n_digits+1] = '\0';
    if (strcmp(buf, cigar) == 0) {
       return 1;
    }
    return 0;
}


/* Call samtools view to convert bam to sam 
   Args: const char* bam_fn - pointer to BAM filename
   Returns: pointer to sam stdout */
FILE* bam_to_sam(const char* bam_fn) {
    char cmd_buf[512];
    sprintf(cmd_buf, "samtools view %s", bam_fn);
    FILE* sam_out = popen(cmd_buf, "r"); // get file pointer to sam stdout
    if (sam_out == NULL) {
        fprintf( stderr, "Error: Unable to open %s with samtools view.\n", bam_fn );
        exit(1);
    }
    return sam_out;
}


/* Get only part of a string
   Args: const char* str - pointer to original string
         char* substr - pointer to char array where substring will be stored
         size_t start - index of where the substring starts in str
         int len - length of substring */
void get_substring(const char* str, char* substr, size_t start, int len) {
    strncpy(substr, str+start, len);
}


/* Takes a pointer to a single alignment in SAM format and performs checks to 
   determine if it passes filters for various fields; if it does, updates the
   counts for the kmers observed in the genome, centered on the 5' & 3'
   fragmentation points of the alignment, to the 5' & 3' Kmers count structures
   (in the case of merged reads).
   If the alignment's reverse flag is on, extract the genome sequence where the
   read maps to (including part of the kmers extending from both ends), reverse
   complement this sequence, then add to the counts.

   Args: KSP fpks - pointer to Kmers structure where 5' kmer counts will be added
         KSP tpks - pointer to Kmers structure where 3' kmer counts will be added
         Genome* genome - pointer to reference genome sequence
         Saml* sp - pointer to alignment
   Returns: 0 if counts added successfully; 
            1 if the reference sequence cannot be found from alignment info;
            2 if the alignment does not pass filters */
int process_aln(KSP fpks, KSP tpks, Genome* genome, Saml* sp) {
    
    Seq* ref = find_seq(genome, sp->rname);
    if (ref == NULL) {
        return 1;
    }

    unsigned long aln_start = (sp->pos)-1; // pos in sam is 1-based
    unsigned long aln_end = aln_start + sp->seq_len - 1;

    // number of context bases outside alignment is always k/2;
    // if k is odd, count one more base inside alignment
    unsigned int ok = KLEN/2;
    unsigned int ik = KLEN-ok;

    if ( aln_start-(KLEN/2) >= 0 && // check for presence of k/2 external bases
         aln_end+(KLEN/2) <= ref->len-1 && 
         sp->mapq >= MIN_MQ &&
         read_len_ok(sp->seq_len) &&
         cigar_ok(sp->seq_len, sp->cigar) &&
         sp->unmap == FALSE &&
         sp->secondary == FALSE &&
         sp->qc_failed == FALSE &&
         sp->duplicate == FALSE &&
         sp->supplementary == FALSE ) {
        
        // process merged/non-paired reads
        if ( sp->paired == FALSE ) {

            int add_5p, add_3p;
            if ( sp->reverse == TRUE ) {
                // allocate space to store part of genome where the alignment
                // maps + k external bases + NULL
                char sub_ref[sp->seq_len+(KLEN+1)];
                get_substring(ref->seq, sub_ref, aln_start-ok, sp->seq_len+KLEN);

                // do reverse complement of that substring
                char rvcmp_sub_ref[sp->seq_len+(KLEN+1)];
                do_rvcmp(sub_ref, rvcmp_sub_ref, sp->seq_len+KLEN);
                
                // the start of the genome substring is already the beginning
                // of the 5' kmer context, so 5' index is 0
                add_5p = add_to_ksp(&(rvcmp_sub_ref[0]), fpks);
                // for 3', starting index = number of context bases outside alignment +
                // length of alignment - number of context bases inside alignment
                add_3p = add_to_ksp(&(rvcmp_sub_ref[ok+sp->seq_len-ik]), tpks);
                if ( (add_5p == 0) && (add_3p == 0) ) {
                    return 0;
                }
                return -1;
            } 
            
            else {
                // for forward reads, use original genome sequence
                add_5p = add_to_ksp(&(ref->seq[aln_start-ok]), fpks);
                add_3p = add_to_ksp(&(ref->seq[aln_start+sp->seq_len-ik]), tpks);
                if ( (add_5p == 0) && (add_3p == 0) ) {
                    return 0;
                }
                return -1;
            }
        }
        
        // process paired reads
        // use 1st read in pair for 5' kmer context, and 2nd read for 3'
        else if ( sp->paired == TRUE &&
                  OPT_M == TRUE &&
                  sp->proper_pair == TRUE &&
                  sp->munmap == FALSE ) {
        
            if ( sp->reverse == TRUE ) {
                char sub_ref[sp->seq_len+(KLEN+1)];
                char rvcmp_sub_ref[sp->seq_len+(KLEN+1)];
                get_substring(ref->seq, sub_ref, aln_start-ok, sp->seq_len+KLEN);
                do_rvcmp(sub_ref, rvcmp_sub_ref, sp->seq_len+KLEN);
    
                if ( sp->read1 == TRUE ) {
                    return add_to_ksp(&(rvcmp_sub_ref[0]), fpks);
                }
                else if ( sp->read2 == TRUE ) {
                    return add_to_ksp(&(rvcmp_sub_ref[ok+sp->seq_len-ik]), tpks);
                }
            }
            else {
                if ( sp->read1 == TRUE ) {
                    return add_to_ksp(&(ref->seq[aln_start-ok]), fpks);
                }
                else if ( sp->read2 == TRUE ) {
                    return add_to_ksp(&(ref->seq[aln_start+sp->seq_len-ik]), tpks);
                }
            }
        }
    }
    return 2;
}


/* Generate all possible kmers of length k and print their counts
   Args: const char bases[] - set of possible bases (A, C, G, T)
         char kmer_tmp[] - string representing the kmer, which starts
                           as empty and has bases recursively added until
                           reaching the length of k, at which point the
                           function returns and prints out the count from
                           the input Kmers structure
         int small_k - integer representing k; with each recursive call, 
                       this number decreases until reaching 0 and the function
                       returns (base case)
         KSP fpks - pointer to 5' Kmers structure
         KSP tpks - pointer to 3' Kmers structure */
void print_kmer_counts(const char bases[], char kmer_tmp[], int small_k, KSP fpks, KSP tpks) {
    if (small_k == 0) {
        kmer_tmp[KLEN] = '\0';
        unsigned int fp_count = kmer2count(&(kmer_tmp[0]), fpks);
        unsigned int tp_count = kmer2count(&(kmer_tmp[0]), tpks);
        printf("%s\t%u\t%u\n", kmer_tmp, fp_count, tp_count);
        return;
    }
    int i, j;
    for (i = 0; i < 4; i++) {
        char new_tmp[100];
        for (j = 0; j < KLEN-small_k; j++) {
            new_tmp[j] = kmer_tmp[j];
        }
        new_tmp[j] = bases[i];
        // recursive call to continue adding bases until length = k
        print_kmer_counts(bases, new_tmp, small_k-1, fpks, tpks);
    }
}


/*** MAIN PROGRAM ***/
int main(int argc, char* argv[]) {
    
    //clock_t start, end;
    //double time_elapsed;
    //start = clock();

    int option;
    char* opts = ":F:B:k:l:L:q:t:m";
    char* fasta_fn, *bam_fn, *ptr1, *ptr2, *utag;
    char user_cmd[MAX_FIELD_WIDTH + 1];
    while ( (option = getopt(argc, argv, opts)) != -1 ) {
        switch (option) {
            case 'F':
                fasta_fn = strdup(optarg);
                break;
            case 'B':
                bam_fn = strdup(optarg);
                break;
            case 'k':
                KLEN = atoi(optarg);
                break;
            case 'l':
                MIN_READ_LEN = strtoul(optarg, &ptr1, 10);
                break;
            case 'L':
                MAX_READ_LEN = strtoul(optarg, &ptr2, 10);
                break;
            case 'q':
                MIN_MQ = atoi(optarg);
                break;
            case 't':
                OPT_T = TRUE;
                utag = strdup(optarg);
                break;
            case 'm':
                OPT_M = TRUE;
                break;
            case ':':
                fprintf( stderr, "Please enter required argument for option -%c.\n", optopt );
                exit(0);
            case '?':
                if ( isprint(optopt) ) {
                    fprintf( stderr, "Unknown option -%c.\n", optopt );
                }
                else {
                    fprintf( stderr, "Unknown option character \\x%x.\n", optopt );
                }
                break;
            default:
                fprintf( stderr, "Error parsing command-line options.\n" );
                exit(0);
        }
    }
    for (int i = optind; i < argc; i++) {
            fprintf( stderr, "Non-option argument %s\n", argv[i] );
    }

    if (!fasta_fn || !bam_fn) {
        fprintf(stderr, "fragkon: Program for describing kmer-based genomic sequence\n");
        fprintf(stderr, "contexts around the fragmentation points of aligned reads.\n" );
        fprintf(stderr, "-F <reference FASTA (required)>\n");
        fprintf(stderr, "-B <input BAM (required)>\n");
        fprintf(stderr, "-k <kmer length (default: 8)>\n");
        fprintf(stderr, "-l <minimum length of read to report (default: 0)>\n");
        fprintf(stderr, "-L <maximum length of read to report (default: 250000000)>\n");
        fprintf(stderr, "-q <map quality filter of read to report (default: 0)>\n" );
        fprintf(stderr, "-t <only consider reads with this optional field (in TAG:TYPE:VALUE format)>\n");
        fprintf(stderr, "-m <only consider merged reads>\n");
        exit(1);
    }

    strcpy(user_cmd, argv[0]);
    strcat(user_cmd, " ");
    for (int i = 1; i < argc; i++) {
        strcat(user_cmd, argv[i]);
        strcat(user_cmd, " ");
    }
    fprintf( stderr,  "# Entered command: %s\n", user_cmd );

    fprintf( stderr, "Input kmer length = %d.\n", KLEN);
    if ((KLEN & 1) != 0) {
        fprintf( stderr, "    *** k is odd - counting %d bases outside %d bases inside of alignment.\n", KLEN/2, (KLEN/2)+1);
    }
    fprintf( stderr, "Reading genome sequence from: %s\n", fasta_fn);
    Genome* genome = init_genome(fasta_fn);
    fprintf( stderr, "Finished loading genome.\nCounting kmer contexts for: %s\n", bam_fn );
    
    KSP fpks = init_KSP(KLEN);
    KSP tpks = init_KSP(KLEN);

    FILE* sam_out = bam_to_sam(bam_fn);
    char saml_buf[MAX_LINE_LEN + 1];
    Saml* sp = malloc(sizeof(Saml));
     
    while ( fgets(saml_buf, MAX_LINE_LEN + 1, sam_out) ) {
        int parse_status = line2saml(saml_buf, sp);
        if (parse_status) {
            if (DEBUG) {
                fprintf( stderr, "Problem parsing alignment, continuing to next entry...\n" );
            }
            continue;
        }

        if ( OPT_T == TRUE ) {
            if ( !has_tag(sp, utag) ) {
                if (DEBUG) {
                    fprintf( stderr, "%s does not have specified tag %s\n", sp->qname, utag);
                }
                continue;
            }
        }

        int process_status = process_aln(fpks, tpks, genome, sp);
        if (process_status == 0) {
            continue;
        }
        else if ( (process_status == 1) && (DEBUG == TRUE) ) {
            fprintf( stderr, "%s: Unable to find sequence %s in genome.\n", sp->qname, sp->rname );
        }
        else if ( (process_status == -1) && (DEBUG == TRUE) ) {
            fprintf( stderr, "%s: Failed to add context counts for this alignment.\n", sp->qname );
        }
        else if ( (process_status == 2) && (DEBUG == TRUE) ) {
            fprintf( stderr, "%s: Alignment did not pass filters.\n", sp->qname );
        }
    }

    char bases[5] = "ACGT";
    char kmer_tmp[100];
    printf("### fragkon.c v0.3\n### %s\n### %s\n", fasta_fn, bam_fn);
    printf("# KMER\t5' CONTEXT COUNTS\t3' CONTEXT COUNTS\n");
    print_kmer_counts(bases, kmer_tmp, KLEN, fpks, tpks);
    free(fasta_fn);
    free(bam_fn);
    free(utag);
    free(sp);
    fclose(sam_out);
    
    destroy_genome(genome);
    destroy_KSP(fpks);
    destroy_KSP(tpks);
    fprintf(stderr, "Done.\n");
    
    //end = clock();
    //time_elapsed = ( (double)(end-start)  / CLOCKS_PER_SEC ) / 60;
    //fprintf( stderr, "Execution time: %f minutes.\n", time_elapsed );
    
    return 0;
}