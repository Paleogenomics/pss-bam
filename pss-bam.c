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


/* Initialize the forward and reverse count matrices with 0 */
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


/* Return the reverse complement of the input sequence */
void do_reverse_complement(const char* seq, char* rev_buf, size_t seq_len) {
    for (size_t i = seq_len-1; i >= 0; i--) {
        char base = seq[i];
        if ( (base == 'A') || (base == 'a') ) {
            rev_buf[seq_len-(i+1)] = 'T';
        }
        else if ( (base == 'C') || (base == 'c') ) {
            rev_buf[seq_len-(i+1)] = 'G';
        }
        else if ( (base == 'G') || (base == 'g') ) {
            rev_buf[seq_len-(i+1)] = 'C';
        }
        else if ( (base == 'T') || (base == 't') ) {
            rev_buf[seq_len-(i+1)] = 'A';
        }
        else {
            rev_buf[seq_len-(i+1)] = base;
        }
    }
}


/* Check if read length is between values specified by -l and -L 
   and that it's >= REGION_LEN */
int read_len_ok(const int seq_len) {
    if ( (seq_len >= MIN_READ_LEN) && 
         (seq_len <= MAX_READ_LEN) &&
         (seq_len >= REGION_LEN) ) {
        return 1;
    }
    return 0;
}


/* Check if CIGAR string only contains a number (that equals 
   the aligned read length) and an M */
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


/* Check if 1st base upstream/downstream of aln is specified in -U/-D */
int context_bases_ok(const char* genome_seq, int read_len) {
    char fb_upstr = genome_seq[1];
    char fb_dwnstr = genome_seq[read_len+2];
    if ( strchr(UPSTR_BASE_CNTXT, fb_upstr) &&
         strchr(DWNSTR_BASE_CNTXT, fb_dwnstr) ) {
        return 1;
    }
    return 0;
}


/* Call samtools view to convert bam to sam */
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


/* Count the 2 context bases upstream/downstream of aln */
int add_context_counts(unsigned long** count_mtrx, const char fcb, const char scb) {
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
    return 0;
}


/* Count base substitutions in the interior of aln in forward direction */
int add_fwd_counts(unsigned long** fwd_counts, const char* genome_seq, \
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
    return 0;
}


/* Count base substitutions in the interior of aln in reverse direction */
int add_rev_counts(unsigned long** rev_counts, const char* genome_seq, \
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
    return 0;
}


/* add_counts_from_aligned_seq
   Takes a genomic sequence with 2 bases of context upstream and downstream,
   the read sequence, and the counts matrices. Goes through the beginning and
   end of the alignment and updates the counts.
   The genome and read sequence must be aligned without gaps. That is, the
   CIGAR string must have been checked to make sure that there are no
   D, I, H, or S characters, only M.
   If the alignment was to the reverse complement of the genome sequence, then
   the read and genome sequence are assumed to already by reverse complemented
   here.
   Args: unsigned long** fwd_counts - pointer to array of Genome/Read base
                                      counts to be added to by this alignment
                                      from the beginning of the alignment
         unsigned long** rev_counts - pointer to array of Genome/Read base
                                      counts to be added to by this alignment
                                      from the end of the alignment
         char* genome_seq - pointer to genome sequence of this alignment with
                            two unaligned bases at the beginning and two unaligned
                            bases at the end
         char* read_seq -   pointer to the aligned read sequence
         int read_len - length of the read sequence

Returns 0 if the alignment was analyzed successfully, adding info to
        fwd_counts and rev_counts.
Returns 1 if there was a problem. */
int add_counts_from_aligned_seq(unsigned long** fwd_counts, unsigned long** rev_counts, \
                                const char* genome_seq, const char* read_seq, int read_len) {
    char fb_upstr = genome_seq[1];
    char sb_upstr = genome_seq[0];
    char fb_dwnstr = genome_seq[read_len+2];
    char sb_dwnstr = genome_seq[read_len+3];

    int cfc = add_context_counts(fwd_counts, fb_upstr, sb_upstr);
    int crc = add_context_counts(rev_counts, fb_dwnstr, sb_dwnstr);
    int fc = add_fwd_counts(fwd_counts, genome_seq, read_seq);
    int rc = add_rev_counts(rev_counts, genome_seq, read_seq, read_len);

    return 0;
}


void get_substring(const char* str, char* substr, size_t start, int len) {
    strncpy(substr, str+start, len);
}


/* Add counts to output matrices from a single read */
int process_aln(unsigned long** fwd_counts, unsigned long** rev_counts, \
                Genome* genome, Saml* sp) {
    
    int read_len = sp->seq_len;

    // allocate char array that will contain only aligned region of genome
    // and 2 context bases upstream/downstream
    char genome_seq[read_len+5]; // 4 context bases total + NULL
    Seq* ref = find_seq(genome, sp->rname);
    if (ref == NULL) {
        return 1;
    }

    char* ref_seq = ref->seq;

    unsigned long aln_start = (sp->pos)-1; // pos in sam is 1-based
    unsigned long aln_end = aln_start+(sp->seq_len)-1;

    if ( (aln_start-2 >= 0) &&   // check for presence of context bases
         (aln_end+2 <= ref->len-1) && 
         (sp->mapq >= MIN_MQ) &&
         (read_len_ok(read_len)) &&
         (cigar_ok(sp->seq_len, sp->cigar)) &&
         (sp->read_paired == FALSE) && 
         (sp->read_unmapped == FALSE) &&
         (sp->not_primary == FALSE) &&
         (sp->not_passing_filters == FALSE) &&
         (sp->is_duplicate == FALSE) &&
         (sp->supplementary_alignment == FALSE) ) {
        
        // copy aligned segment of ref_seq + 4 context bases into genome_seq
        get_substring(ref_seq, genome_seq, aln_start-2, read_len+4);
        
        if (sp->read_reverse == TRUE) {
            // allocate char array to hold reverse complement of genome_seq
            char revcomp_genome_seq[read_len+5];
            do_reverse_complement(genome_seq, revcomp_genome_seq, read_len+4);

            if ( context_bases_ok(revcomp_genome_seq, read_len) ) {  
                char revcomp_read_seq[read_len+1];
                do_reverse_complement(sp->seq, revcomp_read_seq, read_len);
                int add_status = add_counts_from_aligned_seq(fwd_counts, rev_counts, revcomp_genome_seq,
                                                             revcomp_read_seq, read_len);
                return (add_status);
            }     
        }
            
        else if ( context_bases_ok(genome_seq, read_len) ) {
            int add_status = add_counts_from_aligned_seq(fwd_counts, rev_counts,
                                                         genome_seq, sp->seq, read_len);
            return (add_status);
        }
    }

    else {
        return 2;
    }
    return 0;
}


/* Output count matrices to stdout */
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
    char* opts = ":F:B:r:l:L:q:U:D";
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
        exit(1);
    }

    fprintf( stderr, "Reading genome sequence from:\n%s\n", fasta_fn );
    Genome* genome = init_genome(fasta_fn);
    fprintf( stderr, "Finished loading genome.\nCounting matches/mismatches from aligned reads in:\n%s\n", bam_fn );
    
    unsigned long** fwd_counts = init_count_mtrx();
    unsigned long** rev_counts = init_count_mtrx();

    FILE* sam_out = bam_to_sam(bam_fn);
    char saml_buf[MAX_LINE_LEN+1];
    Saml* sp = malloc(sizeof(Saml));
     
    while ( fgets(saml_buf, MAX_LINE_LEN+1, sam_out) ) {
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
            fprintf( stderr, "%s: Unable to find sequence with name %s in genome.\n", sp->qname, sp->rname );
        }
        else if ( (process_status == 2) && (DEBUG == TRUE) ) {
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
