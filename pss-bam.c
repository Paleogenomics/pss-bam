#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <regex.h>

#include "fasta-genome-io.h"
#include "sam-parse.h"

#define TRUE 1
#define FALSE 0

static int REGION_LEN = 15;
static unsigned long MIN_READ_LEN = 0;
static unsigned long MAX_READ_LEN = 250000000;
static int MIN_MQ = 0;
static char* UPSTR_BASE_CNTXT = "ACGT";
static char* DWNSTR_BASE_CNTXT = "ACGT";
static int MERGED_ONLY = 0;


/* Initialize the forward and reverse count matrices with 0 in all cells */
unsigned long** init_count_mtrx() {
    //extra 2 bases upstream/downstream of aln start */
    unsigned long** count_mtrx = malloc((REGION_LEN+2) * sizeof(unsigned long*));
    for (int i = 0; i < REGION_LEN+2; i++) {
        unsigned long* count_arr = malloc(16 * sizeof(unsigned long)); //16 total base substitutions
        for (int j = 0; j < 16; j++) {
            count_arr[j] = (unsigned long)0;
        }
        count_mtrx[i] = count_arr;
    }
    return count_mtrx;
}


/* Return the reverse complement of the input sequence
   EG NOTE: This works only if sequence is all upper-case
   You may want to add cases at the end for lower-case
   sequence characters.
   EG NOTE2: A better design may be to add a field to 
   sam-parse.h called char revcom_seq[MAX_FIELD_WIDTH +1]
   Then, just put the reverse complement sequence in that.
   In this way, you don't have to do a gazillion memory
   allocations (which are slow).
 */
char* do_reverse_complement(const char* seq) {
    unsigned long seq_len = strlen(seq);
    char* rev_comp = malloc(seq_len * sizeof(char)); //remember to free
    for (int i = 0; i < seq_len; i++) {
        char base = seq[i];
        switch (base) {
            case 'A':
                rev_comp[i] = 'T';
                break;
            case 'C':
                rev_comp[i] = 'G';
                break;
            case 'G':
                rev_comp[i] = 'C';
                break;
            case 'T':
                rev_comp[i] = 'A';
                break;
            default:
                rev_comp[i] = base;
        }
    }
    return rev_comp;
}


/* Check if read length is between values specified by -l and -L
   EG NOTE: input argument int seq_len should be const int seq_len
 */
int read_len_ok(int seq_len) {
    if ( (seq_len >= MIN_READ_LEN) && (seq_len <= MAX_READ_LEN) ) {
        return 1;
    }
    return 0;
}


/* Check if CIGAR string only contains a number and an M
   EG NOTE: Using the regular expression engine here is probably overkill and
   slow. I think a simpler approach might be better.
   For example, if you pass this function the length of the sequence as text
   ("73" for example), then the CIGAR string *must* be exactly this string:
   "73M". In that way, you could do a simple strcmp() to determine if it has
   the appropriate CIGAR string
 */
int cigar_ok(const char* cigar) {
    regex_t regex;
    const char* pattern = "^[0-9]+[M]$"; //only a number and an M
    int rc = regcomp(&regex, pattern, REG_EXTENDED);
    if (rc) {
        fprintf(stderr, "Error: Unable to compile regular expression for checking CIGAR string.\n");
        exit(1);
    }
    int status = regexec(&regex, cigar, 0, NULL, 0);
    if (status == 0) {
        regfree(&regex);
        return 1;
    }
    regfree(&regex);
    return 0;
}


/* Check if 1st base upstream/downstream of aln is specified in -U/-D;
    also check if 2 extra bases are present upstream/downstream of aln
    EG NOTE: I think there is a bug here. If I understand what you're doing, 
    you're checked to see if there is a NULL character in the first or second
    position upstream or downstream. If the alignment is at the very beginning
    or end of a chromosome, the characters just outside of that might *NOT*
    be NULL characters! There's no telling what's in those memory positions.
    A better approach would be to check the genome coordinates of those 
    positions. They cannot by less than 0 (for 0-indexed coordinates) 
    or greater than seq->len. In fact, you should check this *before* you
    fetch these characters from the genome. Otherwise, you're pulling in
    data from an unallocated memory region. This is *probably* where your
    segfault is coming from!
 */
int context_base_ok(const char fb_upstr, const char sb_upstr, \
                    const char fb_dwnstr, const char sb_dwnstr) {
    if ( (fb_upstr == '\0')
        | (sb_upstr == '\0') 
        | (fb_dwnstr == '\0')
        | (sb_dwnstr == '\0') ) {
            return 0;
        }
    if ( strchr(UPSTR_BASE_CNTXT, fb_upstr)
        && strchr(DWNSTR_BASE_CNTXT, fb_dwnstr) ) {
        return 1;
    }
    return 0;
}


/* Call samtools view to convert bam to sam */
FILE* bam_to_sam(const char* bam_fn) {
    char cmd_buf[512];
    sprintf(cmd_buf, "samtools view %s", bam_fn);
    FILE* sam_out = popen(cmd_buf, "r"); //get file pointer to sam stdout
    if (sam_out == NULL) {
        fprintf(stderr, "Error: Unable to open %s with samtools view.\n", bam_fn);
        exit(1);
    }
    return sam_out;
}


/* Count the 2 context bases upstream/downstream of aln */
void context_count(unsigned long** count_mtrx, const char fcb, const char scb) {
    char context_bases[3] = {scb, fcb}; //second context base is 1st row
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


/* Count base substitutions in the interior of aln in forward direction */
void fwd_sub_count(unsigned long** fwd_counts, const char* ref_seq, \
                    const char* read_seq, unsigned long aln_start) {
    for (int i = 0; i < REGION_LEN; i++) {
        char pair_bases[3] = {read_seq[i], ref_seq[aln_start+i]};

        if (strcmp(pair_bases, "AA") == 0) {
            fwd_counts[i+2][0] += 1;
        }
        else if (strcmp(pair_bases, "AC") == 0) {
            fwd_counts[i+2][1] += 1;
        }
        else if (strcmp(pair_bases, "AG") == 0) {
            fwd_counts[i+2][2] += 1;
        }
        else if (strcmp(pair_bases, "AT") == 0) {
            fwd_counts[i+2][3] += 1;
        }
        else if (strcmp(pair_bases, "CA") == 0) {
            fwd_counts[i+2][4] += 1;
        }
        else if (strcmp(pair_bases, "CC") == 0) {
            fwd_counts[i+2][5] += 1;
        }
        else if (strcmp(pair_bases, "CG") == 0) {
            fwd_counts[i+2][6] += 1;
        }
        else if (strcmp(pair_bases, "CT") == 0) {
            fwd_counts[i+2][7] += 1;
        }
        else if (strcmp(pair_bases, "GA") == 0) {
            fwd_counts[i+2][8] += 1;
        }
        else if (strcmp(pair_bases, "GC") == 0) {
            fwd_counts[i+2][9] += 1;
        }
        else if (strcmp(pair_bases, "GG") == 0) {
            fwd_counts[i+2][10] += 1;
        }
        else if (strcmp(pair_bases, "GT") == 0) {
            fwd_counts[i+2][11] += 1;
        }
        else if (strcmp(pair_bases, "TA") == 0) {
            fwd_counts[i+2][12] += 1;
        }
        else if (strcmp(pair_bases, "TC") == 0) {
            fwd_counts[i+2][13] += 1;
        }
        else if (strcmp(pair_bases, "TG") == 0) {
            fwd_counts[i+2][14] += 1;
        }
        else if (strcmp(pair_bases, "TT") == 0) {
            fwd_counts[i+2][15] += 1;
        }
        else {
            continue;
        }
    }
}


/* Count base substitutions in the interior of aln in reverse direction */
void rev_sub_count(unsigned long** rev_counts, const char* ref_seq, \
                    const char* read_seq, int seq_len, unsigned long aln_end) {
    for (int i = 0; i < REGION_LEN; i++) {
        char pair_bases[3] = {read_seq[seq_len-(i+1)], ref_seq[aln_end-i]};

        if (strcmp(pair_bases, "AA") == 0) {
            rev_counts[i+2][0] += 1;
        }
        else if (strcmp(pair_bases, "AC") == 0) {
            rev_counts[i+2][1] += 1;
        }
        else if (strcmp(pair_bases, "AG") == 0) {
            rev_counts[i+2][2] += 1;
        }
        else if (strcmp(pair_bases, "AT") == 0) {
            rev_counts[i+2][3] += 1;
        }
        else if (strcmp(pair_bases, "CA") == 0) {
            rev_counts[i+2][4] += 1;
        }
        else if (strcmp(pair_bases, "CC") == 0) {
            rev_counts[i+2][5] += 1;
        }
        else if (strcmp(pair_bases, "CG") == 0) {
            rev_counts[i+2][6] += 1;
        }
        else if (strcmp(pair_bases, "CT") == 0) {
            rev_counts[i+2][7] += 1;
        }
        else if (strcmp(pair_bases, "GA") == 0) {
            rev_counts[i+2][8] += 1;
        }
        else if (strcmp(pair_bases, "GC") == 0) {
            rev_counts[i+2][9] += 1;
        }
        else if (strcmp(pair_bases, "GG") == 0) {
            rev_counts[i+2][10] += 1;
        }
        else if (strcmp(pair_bases, "GT") == 0) {
            rev_counts[i+2][11] += 1;
        }
        else if (strcmp(pair_bases, "TA") == 0) {
            rev_counts[i+2][12] += 1;
        }
        else if (strcmp(pair_bases, "TC") == 0) {
            rev_counts[i+2][13] += 1;
        }
        else if (strcmp(pair_bases, "TG") == 0) {
            rev_counts[i+2][14] += 1;
        }
        else if (strcmp(pair_bases, "TT") == 0) {
            rev_counts[i+2][15] += 1;
        }
        else {
            continue;
        }
    }
}


/* Add counts to output matrices from a single read */
void process_aln(unsigned long** fwd_counts, unsigned long** rev_counts, \
                    Genome* genome, Saml* sp) {
    
    Seq* rseq = find_seq(genome, sp->rname);
    char* ref_seq = rseq->seq;

    unsigned long aln_start = (sp->pos)-1; //pos in sam file is 1-based
    unsigned long aln_end = aln_start+(sp->seq_len)-1;

    if ( (sp->mapq >= MIN_MQ) 
        && (read_len_ok(sp->seq_len))
        && (cigar_ok(sp->cigar)) 
        && (sp->read_paired == FALSE)  
        && (sp->read_unmapped == FALSE) 
        && (sp->not_primary == FALSE) 
        && (sp->not_passing_filters == FALSE)
        && (sp->is_duplicate == FALSE) 
        && (sp->supplementary_alignment == FALSE) ) {
        
        if (sp->read_reverse == TRUE) {
            char* rev_ref = do_reverse_complement(ref_seq);
            char rfb_upstr = rev_ref[aln_start-1]; //reverse 1st base upstream
            char rsb_upstr = rev_ref[aln_start-2]; //reverse 2nd base upstream
            char rfb_dwnstr = rev_ref[aln_end+1];  //downstream
            char rsb_dwnstr = rev_ref[aln_end+2];

            if (context_base_ok(rfb_upstr, rsb_upstr, rfb_dwnstr, rsb_dwnstr)) {
                context_count(fwd_counts, rfb_upstr, rsb_upstr); //count upstream context bases
                context_count(rev_counts, rfb_dwnstr, rsb_dwnstr); //count downstream context bases
                char* rev_read = do_reverse_complement(sp->seq);
                fwd_sub_count(fwd_counts, rev_ref, rev_read, aln_start);
                rev_sub_count(rev_counts, rev_ref, rev_read, sp->seq_len, aln_end);
                free(rev_read);
            }
            free(rev_ref);
        }
            
        else {
            char fb_upstr = ref_seq[aln_start-1];
            char sb_upstr = ref_seq[aln_start-2];
            char fb_dwnstr = ref_seq[aln_end+1];
            char sb_dwnstr = ref_seq[aln_end+2];

            if (context_base_ok(fb_upstr, sb_upstr, fb_dwnstr, sb_dwnstr)) {
                context_count(fwd_counts, fb_upstr, sb_upstr);
                context_count(rev_counts, fb_dwnstr, sb_dwnstr);
                fwd_sub_count(fwd_counts, ref_seq, sp->seq, aln_start);
                rev_sub_count(rev_counts, ref_seq, sp->seq, sp->seq_len, aln_end);
            }
        }
    }

    else {
        printf("Read %s did not pass filters\n", sp->qname);
    }
}


/* Output count tables to to stdout */
void print_output(const char* fasta_fn, const char* bam_fn, \
                    unsigned long** fwd_counts, unsigned long** rev_counts) {
    printf("### pss-bam.c v1.0\n### %s\n### %s\n", fasta_fn, bam_fn);
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
    
    //for rev_counts, print rows in backward order
    for (int i = REGION_LEN-1; i >= 0; i--) {
        printf("%d\t", i);
        for (int j = 0; j < 16; j++) {
            printf("%lu\t", rev_counts[(REGION_LEN+1)-i][j]);
        }
        printf("\n");
    }
    
    //print dowmstream context rows
    for (int i = 0; i < 2; i++) {
        printf("%d\t", i+1);
        for (int j = 0; j < 16; j++) {
            printf("%lu\t", rev_counts[REGION_LEN+i][j]);
        }
        printf("\n");
    }
}


/* Free memory allocated by init_count_table */
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


/* Go through every alignment in bam file to populate count matrices
   EG NOTE: I would have this in main(). This is the main control loop
   of the program.
*/
void do_count(const char* fasta_fn, const char* bam_fn) {
    Genome* genome = init_genome(fasta_fn);
    unsigned long** fwd_counts = init_count_mtrx();
    unsigned long** rev_counts = init_count_mtrx();

    FILE* sam_out = bam_to_sam(bam_fn);
    char saml_buf[MAX_LINE_LEN+1];
    Saml* sp = malloc(sizeof(Saml));
     
    while (fgets(saml_buf, MAX_LINE_LEN+1, sam_out)) {
        int status = line2saml(saml_buf, sp);
        if (status) {
            fprintf(stderr, "Problem parsing alignment entry, continuing to next entry...\n");
            continue;
        }
        process_aln(fwd_counts, rev_counts, genome, sp);
    }
    print_output(fasta_fn, bam_fn, fwd_counts, rev_counts);

    destroy_genome(genome);
    fclose(sam_out);
    free(sp);
    destroy_count_mtrx(fwd_counts);
    destroy_count_mtrx(rev_counts);
}



/*** MAIN PROGRAM ***/
int main(int argc, char* argv[]) {
    int option;
    char* opts = ":F:B:r:l:L:q:U:D";
    char* fasta_fn, *bam_fn, *ptr1, *ptr2;
    while ((option = getopt(argc, argv, opts)) != -1) {
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
                fprintf(stderr, "Please enter required argument for option -%c.\n", optopt);
                exit(0);
            case '?':
                if (isprint(optopt)) {
                    fprintf(stderr, "Unknown option -%c.\n", optopt);
                }
                else {
                    fprintf(stderr, "Unknown option character \\x%x.\n", optopt);
                }
                break;
            default:
                fprintf(stderr, "Error parsing command-line options.\n");
                exit(0);
        }
    }
    for (int i = optind; i < argc; i++) {
            printf("Non-option argument %s\n", argv[i]);
    }

    if (!fasta_fn || !bam_fn) {
        fprintf(stderr, "pss-bam: Program for describing base context and counting\n");
        fprintf(stderr, "the number of matches/mismatches in aligned reads to a genome.\n" );
        fprintf(stderr, "-F <reference FASTA>\n");
        fprintf(stderr, "-B <input BAM>\n");
        fprintf(stderr, "-r <length in basepairs into the interior of alignments to report on (default: 15)>\n");
        fprintf(stderr, "-l <minimum length of read to report (default: 0)>\n");
        fprintf(stderr, "-L <maximum length of read to report (default: 250000000)>\n");
        fprintf(stderr, "-q <map quality filter of read to report (default: 0)>\n" );
        fprintf(stderr, "-U <upstream context base filter; first base before alignment must be one of these (default: ACGT)>\n");
        fprintf(stderr, "-D <downstream context base filter; first base before alignment must be one of these (default: ACGT)>\n");
        exit(1);
    }

    do_count(fasta_fn, bam_fn);
    free(fasta_fn);
    free(bam_fn);
    return 0;
}
