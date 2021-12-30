#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "fasta-genome-io.h"
#include "sam-parse.h"

#define TRUE (1)
#define FALSE (0)
#define DEBUG (0)

static int REGION_LEN = 15;
static unsigned long MIN_READ_LEN = 0;
static unsigned long MAX_READ_LEN = 250000000;
static int MIN_MQ = 0;
static char* UPSTR_BASE_CNTXT = "ACGT";
static char* DWNSTR_BASE_CNTXT = "ACGT";
static unsigned int MERGED_ONLY = FALSE;


/* Initialize the forward and reverse count matrices with 0
   Returns: A 2D unsigned long matrix with size (REGION_LEN+2) x 16;
            all cells initialized to 0 */
unsigned long** init_count_mtrx() {

    // extra 2 bases upstream/downstream of aln start */
    unsigned long** count_mtrx = malloc( (REGION_LEN+2) * sizeof(unsigned long*) );
    for (int i = 0; i < REGION_LEN+2; i++) {
        unsigned long* count_arr = malloc( 16 * sizeof(unsigned long) ); // 16 total base substitutions
        for (int j = 0; j < 16; j++) {
            count_arr[j] = (unsigned long)0;
        }
        count_mtrx[i] = count_arr;
    }
    return count_mtrx;
}


/* Get the reverse complement of a sequence
Args: const char* seq - input sequence that we want to get the
                           the reverse complement of
      char* revcomp - pointer to char array where the reverse
                         complemented sequence is stored
      size_t seq_len - length of the input sequence */ 
void do_reverse_complement(const char* seq, char* revcomp, size_t seq_len) {
    for (int i = seq_len-1; i >= 0; i--) {
        char base = seq[i];
        if ( (base == 'A') || (base == 'a') ) {
            revcomp[seq_len-(i+1)] = 'T';
        }
        else if ( (base == 'C') || (base == 'c') ) {
            revcomp[seq_len-(i+1)] = 'G';
        }
        else if ( (base == 'G') || (base == 'g') ) {
            revcomp[seq_len-(i+1)] = 'C';
        }
        else if ( (base == 'T') || (base == 't') ) {
            revcomp[seq_len-(i+1)] = 'A';
        }
        else {
            revcomp[seq_len-(i+1)] = base;
        }
    }
}


/* Check if read length is between values specified by -l and -L 
   and whether is larger or equal to REGION_LEN
   Args: const int seq_len - length of the read
   Returns: 1 if passes the above conditions; 0 otherwise */
int read_len_ok(const int seq_len) {
    if ( (seq_len >= MIN_READ_LEN) && 
         (seq_len <= MAX_READ_LEN) &&
         (seq_len >= REGION_LEN) ) {
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


/* Check if 1st base upstream/downstream of alignment is specified by -U/-D 
   Args: const char* genome_seq - pointer to genome sequence of this alignment with
                                  2 unaligned bases at the beginning and 2 unaligned
                                  bases at the end
         int read_len - length of the alignment 
   Returns: 1 if the first base before alignment is in UPSTR_BASE_CONTXT and
            first base after alignment is in DWNSTR_BASE_CNTXT;
            0 otherwise */
int context_bases_ok(const char* genome_seq, int read_len) {
    char fb_upstr = genome_seq[1];
    char fb_dwnstr = genome_seq[read_len+2];
    if ( strchr(UPSTR_BASE_CNTXT, fb_upstr) &&
         strchr(DWNSTR_BASE_CNTXT, fb_dwnstr) ) {
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


/* Add counts from the 2 context bases upstream/downstream of alignment 
   Args: unsigned long** count_mtrx - pointer to matrix where the counts will be added
         const char fcb - first context base either before or after alignment
         const char scb - second context base either before or after alignment 
   Returns: 0 if counts added successfully */
void add_context_counts(unsigned long** count_mtrx, const char fcb, const char scb) {
    char context_bases[3] = {scb, fcb}; // second context base is 1st row
    for (int i = 0; i < 2; i++) {
        switch (context_bases[i]) {
            case 'A':
                count_mtrx[i][0] += 1;
                break;
            case 'C':
                count_mtrx[i][5] += 1;
                break;
            case 'G':
                count_mtrx[i][10] += 1;
                break;
            case 'T':
                count_mtrx[i][15] += 1;
                break;
            default:
                continue;
        }
    }
}


/* Add base substitution counts from interior of alignment in the forward direction 
   Args: unsigned long** fwd_counts - pointer to matrix where forward counts will be added 
         const char* genome_seq - pointer to genome sequence with 2 unaligned bases at 
                                  the beginning and 2 unaligned bases at the end
         const char* read_seq - pointer to read sequence
   Returns: 0 if counts added successfully */
void add_fwd_counts(unsigned long** fwd_counts, const char* genome_seq, \
                   const char* read_seq) {
    for (int i = 0; i < REGION_LEN; i++) {

        char ref_base = genome_seq[2+i];
        char read_base = toupper(read_seq[i]);
        char pair[3] = {read_base, ref_base};

        if (strcmp(pair, "AA") == 0) {
            fwd_counts[i+2][0] += 1;
        }
        else if (strcmp(pair, "AC") == 0) {
            fwd_counts[i+2][1] += 1;
        }
        else if (strcmp(pair, "AG") == 0) {
            fwd_counts[i+2][2] += 1;
        }
        else if (strcmp(pair, "AT") == 0) {
            fwd_counts[i+2][3] += 1;
        }
        else if (strcmp(pair, "CA") == 0) {
            fwd_counts[i+2][4] += 1;
        }
        else if (strcmp(pair, "CC") == 0) {
            fwd_counts[i+2][5] += 1;
        }
        else if (strcmp(pair, "CG") == 0) {
            fwd_counts[i+2][6] += 1;
        }
        else if (strcmp(pair, "CT") == 0) {
            fwd_counts[i+2][7] += 1;
        }
        else if (strcmp(pair, "GA") == 0) {
            fwd_counts[i+2][8] += 1;
        }
        else if (strcmp(pair, "GC") == 0) {
            fwd_counts[i+2][9] += 1;
        }
        else if (strcmp(pair, "GG") == 0) {
            fwd_counts[i+2][10] += 1;
        }
        else if (strcmp(pair, "GT") == 0) {
            fwd_counts[i+2][11] += 1;
        }
        else if (strcmp(pair, "TA") == 0) {
            fwd_counts[i+2][12] += 1;
        }
        else if (strcmp(pair, "TC") == 0) {
            fwd_counts[i+2][13] += 1;
        }
        else if (strcmp(pair, "TG") == 0) {
            fwd_counts[i+2][14] += 1;
        }
        else if (strcmp(pair, "TT") == 0) {
            fwd_counts[i+2][15] += 1;
        }
        else {
            continue;
        }
    }
}


/* Add base substitution counts from interior of alignment in the reverse direction 
   Args: unsigned long** rev_counts - pointer to matrix where reverse counts will be added 
         const char* genome_seq - pointer to genome sequence with 2 unaligned bases at 
                                  the beginning and 2 unaligned bases at the end
         const char* read_seq - pointer to the aligned read sequence
         int read_len - length of read
   Returns: 0 if counts added successfully */
void add_rev_counts(unsigned long** rev_counts, const char* genome_seq, \
                    const char* read_seq, int read_len) {
    for (int i = 0; i < REGION_LEN; i++) {

        char ref_base = genome_seq[read_len+1-i];
        char read_base = toupper(read_seq[read_len-1-i]);
        char pair[3] = {read_base, ref_base};

        if (strcmp(pair, "AA") == 0) {
            rev_counts[i+2][0] += 1;
        }
        else if (strcmp(pair, "AC") == 0) {
            rev_counts[i+2][1] += 1;
        }
        else if (strcmp(pair, "AG") == 0) {
            rev_counts[i+2][2] += 1;
        }
        else if (strcmp(pair, "AT") == 0) {
            rev_counts[i+2][3] += 1;
        }
        else if (strcmp(pair, "CA") == 0) {
            rev_counts[i+2][4] += 1;
        }
        else if (strcmp(pair, "CC") == 0) {
            rev_counts[i+2][5] += 1;
        }
        else if (strcmp(pair, "CG") == 0) {
            rev_counts[i+2][6] += 1;
        }
        else if (strcmp(pair, "CT") == 0) {
            rev_counts[i+2][7] += 1;
        }
        else if (strcmp(pair, "GA") == 0) {
            rev_counts[i+2][8] += 1;
        }
        else if (strcmp(pair, "GC") == 0) {
            rev_counts[i+2][9] += 1;
        }
        else if (strcmp(pair, "GG") == 0) {
            rev_counts[i+2][10] += 1;
        }
        else if (strcmp(pair, "GT") == 0) {
            rev_counts[i+2][11] += 1;
        }
        else if (strcmp(pair, "TA") == 0) {
            rev_counts[i+2][12] += 1;
        }
        else if (strcmp(pair, "TC") == 0) {
            rev_counts[i+2][13] += 1;
        }
        else if (strcmp(pair, "TG") == 0) {
            rev_counts[i+2][14] += 1;
        }
        else if (strcmp(pair, "TT") == 0) {
            rev_counts[i+2][15] += 1;
        }
        else {
            continue;
        }
    }
}


/* Takes a genomic sequence with 2 bases of context upstream and downstream,
   the read sequence, and the count matrices. Goes through the beginning and
   end of the alignment and updates the counts.
   The genome and read sequence must be aligned without gaps. That is, the
   CIGAR string must have been checked to make sure that there are no
   D, I, H, or S characters, only M.
   If the alignment was the reverse complement of the genome sequence, then
   the read and genome sequence are assumed to already be reverse complemented.

   Args: unsigned long** fwd_counts - pointer to matrix where Genome/Read base
                                      counts from the beginning of this alignment
                                      will be added
         unsigned long** rev_counts - pointer to matrix where Genome/Read base
                                      counts from the end of this alignment
                                      will be added
         const char* genome_seq - pointer to genome sequence of this alignment with
                                  2 unaligned bases at the beginning and 2 unaligned
                                  bases at the end
         const char* read_seq - pointer to the aligned read sequence
         int read_len - length of the read */
void add_counts_from_aligned_seq(unsigned long** fwd_counts, unsigned long** rev_counts, \
                                const char* genome_seq, const char* read_seq, int read_len) {
    char fb_upstr = genome_seq[1];
    char sb_upstr = genome_seq[0];
    char fb_dwnstr = genome_seq[read_len+2];
    char sb_dwnstr = genome_seq[read_len+3];

    add_context_counts(fwd_counts, fb_upstr, sb_upstr);
    add_context_counts(rev_counts, fb_dwnstr, sb_dwnstr);
    add_fwd_counts(fwd_counts, genome_seq, read_seq);
    add_rev_counts(rev_counts, genome_seq, read_seq, read_len);
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
   determine if it passes filters for various fields; if it does, updates counts
   to both the forward and reverse matrices (in the case of merged reads) using 
   information from the genome and aligned read sequences. If the alignment's 
   reverse flag is on, reverse complements both the genome and read sequences 
   before adding counts.
   
   Args: unsigned long** fwd_counts - pointer to matrix where counts from the
                                      beginning of alignment will be added
         unsigned long** rev_counts - pointer to matrix where counts from the
                                      end of alignment will be added
         Genome* genome - pointer to reference genome sequence
         Saml* sp - pointer to alignment
   Returns: 0 if counts added successfully; 
            1 if the reference sequence cannot be found from alignment info;
            2 if the alignment does not pass filters */
int process_aln(unsigned long** fwd_counts, unsigned long** rev_counts, \
                Genome* genome, Saml* sp) {

    // allocate char array that will contain only aligned region of genome
    // and 2 context bases upstream/downstream
    char genome_seq[sp->seq_len+5]; // 4 context bases total + NULL
    Seq* ref = find_seq(genome, sp->rname);
    if (ref == NULL) {
        return 1;
    }

    char* ref_seq = ref->seq;

    unsigned long aln_start = (sp->pos)-1; // pos in sam is 1-based
    unsigned long aln_end = aln_start + sp->seq_len - 1;

    if ( aln_start-2 >= 0 &&   // check for presence of context bases
         aln_end+2 <= ref->len-1 && 
         sp->mapq >= MIN_MQ &&
         read_len_ok(sp->seq_len) &&
         cigar_ok(sp->seq_len, sp->cigar) &&
         sp->unmap == FALSE &&
         sp->secondary == FALSE &&
         sp->qc_failed == FALSE &&
         sp->duplicate == FALSE &&
         sp->supplementary == FALSE ) {
        
        // copy aligned segment of ref_seq + 4 context bases into genome_seq
        get_substring(ref_seq, genome_seq, aln_start-2, sp->seq_len+4);
        
        // process merged/non-paired reads
        if ( sp->paired == FALSE ) {

            if ( sp->reverse == TRUE ) {
                // allocate char array to hold reverse complement of genome_seq
                char revcomp_genome_seq[sp->seq_len+5];
                do_reverse_complement(genome_seq, revcomp_genome_seq, sp->seq_len+4);

                if ( context_bases_ok(revcomp_genome_seq, sp->seq_len) ) {  
                    char revcomp_read_seq[sp->seq_len+1];
                    do_reverse_complement(sp->seq, revcomp_read_seq, sp->seq_len);
                    add_counts_from_aligned_seq(fwd_counts, rev_counts, revcomp_genome_seq,
                                                revcomp_read_seq, sp->seq_len);
                    return 0;
                }     
            }
            else if ( sp->reverse == FALSE &&
                      context_bases_ok(genome_seq, sp->seq_len) ) {
                add_counts_from_aligned_seq(fwd_counts, rev_counts,
                                            genome_seq, sp->seq, sp->seq_len);
                return 0;
            }
        }
        
        // process paired reads
        else if ( sp->paired == TRUE &&
                  MERGED_ONLY == FALSE &&
                  sp->proper_pair == TRUE &&
                  sp->munmap == FALSE ) {
            
            if ( sp->reverse == TRUE ) {
                char revcomp_genome_seq[sp->seq_len+5];
                do_reverse_complement(genome_seq, revcomp_genome_seq, sp->seq_len+4);
                
                // if 1st read in pair, only add forward & upstream
                // context counts
                if ( sp->read1 == TRUE &&
                     strchr(UPSTR_BASE_CNTXT, revcomp_genome_seq[1]) ) {
                    char revcomp_read_seq[sp->seq_len+1];
                    do_reverse_complement(sp->seq, revcomp_read_seq, sp->seq_len);
                    add_context_counts(fwd_counts, revcomp_genome_seq[1],
                                       revcomp_genome_seq[0]);
                    add_fwd_counts(fwd_counts, revcomp_genome_seq, revcomp_read_seq);
                    return 0;
                }

                // if 2nd read in pair, only add reverse & downstream
                // context counts
                else if ( sp->read2 == TRUE &&
                          strchr(DWNSTR_BASE_CNTXT, revcomp_genome_seq[sp->seq_len+2]) ) {
                    char revcomp_read_seq[sp->seq_len+1];
                    do_reverse_complement(sp->seq, revcomp_read_seq, sp->seq_len);
                    add_context_counts(rev_counts, revcomp_genome_seq[sp->seq_len+2],
                                       revcomp_genome_seq[sp->seq_len+3]);
                    add_rev_counts(rev_counts, revcomp_genome_seq, revcomp_read_seq, sp->seq_len);
                    return 0;
                }
            }
            
            // paired & non-reverse-complemented reads
            else if ( sp->reverse == FALSE &&
                      sp->read1 == TRUE &&
                      strchr(UPSTR_BASE_CNTXT, genome_seq[1]) ) {
                    add_context_counts(fwd_counts, genome_seq[1], genome_seq[0]);
                    add_fwd_counts(fwd_counts, genome_seq, sp->seq);
                    return 0;
            }

            else if ( sp->reverse == FALSE &&
                      sp->read2 == TRUE &&
                      strchr(DWNSTR_BASE_CNTXT, genome_seq[sp->seq_len+2]) ) {
                    add_context_counts(rev_counts, genome_seq[sp->seq_len+2], 
                                       genome_seq[sp->seq_len+3]);
                    add_rev_counts(rev_counts, genome_seq, sp->seq, sp->seq_len);
                    return 0;
            }
        }
    }
    return -1;
}


/* Print count matrices to stdout 
   Args: const char* fasta_fn - pointer to FASTA filename
         const char* bam_fn - pointer to BAM filename 
         unsigned long** fwd_counts - pointer to forward count matrix 
         unsigned long** rev_counts - pointer to reverse count matrix */
void print_output(const char* fasta_fn, const char* bam_fn, \
                    unsigned long** fwd_counts, unsigned long** rev_counts) {
    printf("### pss-bam.c v1.1\n### %s\n### %s\n", fasta_fn, bam_fn);
    puts("### Format of table:");
    puts("### Counts of how often a read base and genome base were seen at");
    puts("### each position in the aligned reads.");
    puts("### First base is what was seen in the read.");
    puts("### Second base is what was in the genome at that position.");
    puts("### POS AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT");
    puts("### Forward read substitution counts and base context");

    for (int i = -2; i < REGION_LEN; i++) {
        printf("%d\t", i);
        for (int j = 0; j < 16; j++) {
            printf("%lu\t", fwd_counts[i+2][j]);
        }
        printf("\n");
    }
    puts("\n\n### Reverse read substitution counts and base context");

    // print all rows except context rows in rev_counts in reverse order
    for (int i = REGION_LEN-1; i >= 0; i--) {
        printf("%d\t", i);
        for (int j = 0; j < 16; j++) {
            printf( "%lu\t", rev_counts[i+2][j] );
        }
        printf( "\n" );
    }
    
    // print context rows in rev_counts
    for (int i = 1; i < 3; i++) {
        printf("%d\t", i);
        for (int j = 0; j < 16; j++) {
            printf("%lu\t", rev_counts[2-i][j]);
        }
        printf( "\n" );
    }
}


/* Free memory allocated by init_count_mtrx() */
int destroy_count_mtrx(unsigned long** count_mtrx) {
    if (!count_mtrx) {
        return 0;
    }
    for (int i = 0; i < REGION_LEN+2; i++) {
        free(count_mtrx[i]);
    }
    free(count_mtrx);
    return 0;
}


/*** MAIN PROGRAM ***/
int main(int argc, char* argv[]) {
    
    //clock_t start, end;
    //double time_elapsed;
    //start = clock();

    int option;
    char* opts = ":F:B:r:l:L:q:U:D:m";
    char* fasta_fn, *bam_fn, *ptr1, *ptr2;
    while ( (option = getopt(argc, argv, opts)) != -1 ) {
        switch (option) {
            case 'F':
                fasta_fn = strdup(optarg);
                break;
            case 'B':
                bam_fn = strdup(optarg);
                break;
            case 'r':
                REGION_LEN = atoi(optarg);
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
            case 'U':
                UPSTR_BASE_CNTXT = optarg;
                break;
            case 'D':
                DWNSTR_BASE_CNTXT = optarg;
                break;
            case 'm':
                MERGED_ONLY = TRUE;
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
        fprintf(stderr, "pss-bam: Program for describing base context and counting\n");
        fprintf(stderr, "the number of matches/mismatches in aligned reads to a genome.\n" );
        fprintf(stderr, "-F <reference FASTA (required)>\n");
        fprintf(stderr, "-B <input BAM (required)>\n");
        fprintf(stderr, "-r <length in basepairs into the interior of alignments to report on (default: 15)>\n");
        fprintf(stderr, "-l <minimum length of read to report (default: 0)>\n");
        fprintf(stderr, "-L <maximum length of read to report (default: 250000000)>\n");
        fprintf(stderr, "-q <map quality filter of read to report (default: 0)>\n" );
        fprintf(stderr, "-U <upstream context base filter; first base before alignment must be one of these (default: ACGT)>\n");
        fprintf(stderr, "-D <downstream context base filter; first base before alignment must be one of these (default: ACGT)>\n");
        fprintf(stderr, "-m <only consider merged reads>\n");
        exit(1);
    }

    if (MERGED_ONLY == TRUE) {
        fprintf( stderr,
                 "Full command: %s -F %s -B %s -r %d -l %lu -L %lu -q %d -U %s -D %s -m\n", 
                 argv[0], fasta_fn, bam_fn, REGION_LEN, MIN_READ_LEN, MAX_READ_LEN, MIN_MQ, UPSTR_BASE_CNTXT, DWNSTR_BASE_CNTXT);
    }

    else {
        fprintf( stderr,
                 "Full command: %s -F %s -B %s -r %d -l %lu -L %lu -q %d -U %s -D %s\n", 
                 argv[0], fasta_fn, bam_fn, REGION_LEN, MIN_READ_LEN, MAX_READ_LEN, MIN_MQ, UPSTR_BASE_CNTXT, DWNSTR_BASE_CNTXT);
    }

    fprintf( stderr, "Reading genome sequence from:\n%s\n", fasta_fn);
    Genome* genome = init_genome(fasta_fn);
    fprintf( stderr, "Finished loading genome.\nCounting matches/mismatches from:\n%s\n", bam_fn );
    
    unsigned long** fwd_counts = init_count_mtrx();
    unsigned long** rev_counts = init_count_mtrx();

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

        int process_status = process_aln(fwd_counts, rev_counts, genome, sp);
        if (process_status == 0) {
            continue;
        }
        else if ( (process_status == 1) && (DEBUG == TRUE) ) {
            fprintf( stderr, "%s: Unable to find sequence %s in genome.\n", sp->qname, sp->rname );
        }
        else if ( (process_status == -1) && (DEBUG == TRUE) ) {
            fprintf( stderr, "%s: Alignment did not pass filters.\n", sp->qname );
        }
    }

    print_output(fasta_fn, bam_fn, fwd_counts, rev_counts);
    
    free(fasta_fn);
    free(bam_fn);
    free(sp);
    fclose(sam_out);
    
    destroy_genome(genome);
    destroy_count_mtrx(fwd_counts);
    destroy_count_mtrx(rev_counts);

    fprintf(stderr, "Done.\n");
    
    //end = clock();
    //time_elapsed = ( (double)(end-start)  / CLOCKS_PER_SEC ) / 60;
    //fprintf( stderr, "Execution time: %f minutes.\n", time_elapsed );
    
    return 0;
}