#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "fasta-genome-io.h"
#include "sam-parse.h"

#define DEBUG (0)

static int REGION_LEN = 15;
static unsigned long MIN_READ_LEN = 0;
static unsigned long MAX_READ_LEN = 250000000;
static int MIN_MQ = 0;
static char* UP_CTX = "ACGT";
static char* DOWN_CTX = "ACGT";
static unsigned int MERGED_ONLY = 0;


/* Initialize the forward and reverse count matrices with 0
   Returns: A 2D unsigned long matrix with size (REGION_LEN+2) x 16;
            all cells initialized to 0 */
unsigned long** init_count_mtrx() {
    // extra 2 bases upstream/downstream of aln start */
    unsigned long** count_mtrx = malloc( (REGION_LEN+2) * sizeof(unsigned long*) );
    for (int i = 0; i < REGION_LEN+2; i++) {
        unsigned long* count_arr = malloc( 16 * sizeof(unsigned long) ); // 16 total nuc. combinations
        for (int j = 0; j < 16; j++) {
            count_arr[j] = (unsigned long)0;
        }
        count_mtrx[i] = count_arr;
    }
    return count_mtrx;
}


/* Initialize the forward and reverse rate matrices with 0
   Returns: A 2D double matrix with size REGION_LEN x 12;
            all cells initialized to 0 */
double** init_rate_mtrx() {
    double** rate_mtrx = malloc( REGION_LEN * sizeof(double*) );
    for (int i = 0; i < REGION_LEN; i++) {
        double* rate_arr = malloc( 12 * sizeof(double) ); // 12 different-base substitutions
        for (int j = 0; j < 12; j++) {
            rate_arr[j] = (double)0;
        }
        rate_mtrx[i] = rate_arr;
    }
    return rate_mtrx;
}


/* Get the reverse complement of a sequence
Args: const char* seq - input sequence that we want to get the
                           the reverse complement of
      char* revcomp - pointer to char array where the reverse
                         complemented sequence is stored
      size_t seq_len - length of the input sequence */ 
void do_revcomp(const char* seq, char* revcomp, size_t seq_len) {
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


/* Convert string to uppercase (in place)
   Args: char* str - pointer to string to be converted */
void to_uppercase(char* str) {
    while (*str) {
        *str = toupper((unsigned char) *str);
        str++;
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
            first base after alignment is in DOWN_CTX;
            0 otherwise */
int context_bases_ok(const char* genome_seq, int read_len) {
    char fb_upstr = genome_seq[1];
    char fb_dwnstr = genome_seq[read_len+2];
    if ( strchr(UP_CTX, fb_upstr) &&
         strchr(DOWN_CTX, fb_dwnstr) ) {
        return 1;
    }
    return 0;
}


/* Call samtools view to convert bam to sam 
   Args: const char* bam_fn - pointer to BAM filename
   Returns: pointer to sam stdout */
FILE* bam_to_sam(const char* bam_fn, char* read_group) {
    char cmd_buf[512];
    if ( read_group == NULL ) {
      sprintf(cmd_buf, "samtools view %s", bam_fn);
    }
    else {
      sprintf(cmd_buf, "samtools view -r %s %s", read_group, bam_fn);
    }
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
         const char scb - second context base either before or after alignment */
void add_ctx_counts(unsigned long** count_mtrx, const char fcb, const char scb) {
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
         const char* read_seq - pointer to read sequence */
void add_fwd_counts(unsigned long** fwd_counts, const char* genome_seq, \
                   const char* read_seq) {
    for (int i = 0; i < REGION_LEN; i++) {

        char ref_base = genome_seq[2+i];
        char read_base = read_seq[i];
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
         int read_len - length of read */
void add_rev_counts(unsigned long** rev_counts, const char* genome_seq, \
                    const char* read_seq, int read_len) {
    for (int i = 0; i < REGION_LEN; i++) {

        char ref_base = genome_seq[read_len+1-i];
        char read_base = read_seq[read_len-1-i];
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
void add_counts_from_seq(unsigned long** fwd_counts, unsigned long** rev_counts, \
                         const char* genome_seq, const char* read_seq, int read_len) {
    char fb_upstr = genome_seq[1];
    char sb_upstr = genome_seq[0];
    char fb_dwnstr = genome_seq[read_len+2];
    char sb_dwnstr = genome_seq[read_len+3];

    add_ctx_counts(fwd_counts, fb_upstr, sb_upstr);
    add_ctx_counts(rev_counts, fb_dwnstr, sb_dwnstr);
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

    Seq* ref = find_seq(genome, sp->rname);
    if (ref == NULL) {
        return 1;
    }
    char* ref_seq = ref->seq;
                    
    // allocate char array that will contain only aligned region of genome
    // and 2 context bases upstream/downstream
    int seq_len = abs(sp->isize);
    char genome_seq[seq_len+5]; // 4 context bases total + NULL
    long aln_start = (sp->pos)-1; // pos in sam is 1-based
    long aln_end = aln_start + seq_len - 1;

    // apply filters
    if ( (aln_start-2 < 0) ||
         (aln_end+2 > ref->len-1) ||
         (sp->mapq < MIN_MQ) ||
         !read_len_ok(seq_len) ||
         !cigar_ok(seq_len, sp->cigar) ||
         sp->unmap ||
         sp->secondary ||
         sp->qc_failed ||
         sp->duplicate ||
         sp->supplementary ||
         (MERGED_ONLY && sp->paired) ) {
        
        return -1;
    }
        
    // copy aligned segment of ref_seq + 4 context bases into genome_seq
    get_substring(ref_seq, genome_seq, aln_start-2, seq_len+4);
    to_uppercase(genome_seq);
    to_uppercase(sp->seq);
        
    // process merged/non-paired reads
    if ( !sp->paired ) {

        if ( sp->reverse ) {
            // allocate char array to hold reverse complement of genome_seq
            char rc_genome_seq[seq_len+5];
            do_revcomp(genome_seq, rc_genome_seq, seq_len+4);
            if ( context_bases_ok(rc_genome_seq, seq_len) ) {  
                char rc_read_seq[seq_len+1];
                do_revcomp(sp->seq, rc_read_seq, seq_len);
                add_counts_from_seq(fwd_counts, rev_counts, rc_genome_seq,
                                    rc_read_seq, seq_len);
                return 0;
            }
        }

        else if ( context_bases_ok(genome_seq, sp->seq_len) ) {
            add_counts_from_seq(fwd_counts, rev_counts, genome_seq, sp->seq, seq_len);
            return 0;
        }
    }
        
    // process paired reads
    else if ( sp->paired &&
              sp->proper_pair &&
              !sp->munmap ) {
            
        if ( sp->reverse ) {
            char rc_genome_seq[seq_len+5];
            do_revcomp(genome_seq, rc_genome_seq, seq_len+4);
                
            // if 1st read in pair, only add forward & upstream
            // context counts
            if ( sp->read1 &&
                 strchr(UP_CTX, rc_genome_seq[1]) ) {
                char rc_read_seq[seq_len+1];
                do_revcomp(sp->seq, rc_read_seq, seq_len);
                add_ctx_counts(fwd_counts, rc_genome_seq[1], rc_genome_seq[0]);
                add_fwd_counts(fwd_counts, rc_genome_seq, rc_read_seq);
                return 0;
            }

            // if 2nd read in pair, only add reverse & downstream
            // context counts
            else if ( sp->read2 &&
                      strchr(DOWN_CTX, rc_genome_seq[seq_len+2]) ) {
                char rc_read_seq[seq_len+1];
                do_revcomp(sp->seq, rc_read_seq, seq_len);
                add_ctx_counts(rev_counts, rc_genome_seq[seq_len+2], rc_genome_seq[seq_len+3]);
                add_rev_counts(rev_counts, rc_genome_seq, rc_read_seq, seq_len);
                return 0;
            }
        }

        // paired & non-reverse-complemented reads
        else if ( sp->read1 &&
                  strchr(UP_CTX, genome_seq[1]) ) {
                add_ctx_counts(fwd_counts, genome_seq[1], genome_seq[0]);
                add_fwd_counts(fwd_counts, genome_seq, sp->seq);
                return 0;
        }
        else if ( sp->read2 &&
                  strchr(DOWN_CTX, genome_seq[seq_len+2]) ) {
                add_ctx_counts(rev_counts, genome_seq[seq_len+2], genome_seq[seq_len+3]);
                add_rev_counts(rev_counts, genome_seq, sp->seq, seq_len);
                return 0;
        }
    }
    return -1;
}


/* Calculate substitution rates for each different-base substitution (so no AA/CC/GG/TT),
   using counts from the input count matrix
   Args: unsigned long** count_mtrx - pointer to count matrix
         double** rate_mtrx         - pointer to rate matrix where substitution
                                      rates will be added */
void find_sub_rates(unsigned long** count_mtrx, double** rate_mtrx) {
    double n_A = 0, n_C = 0, n_G = 0, n_T = 0;
    for (int i = 0; i < REGION_LEN; i++) {
        int ci = i+2; // skip the 2 context rows in count_mtrx
        n_A = count_mtrx[ci][0] + count_mtrx[ci][4] + count_mtrx[ci][8] + count_mtrx[ci][12];
        n_C = count_mtrx[ci][1] + count_mtrx[ci][5] + count_mtrx[ci][9] + count_mtrx[ci][13];
        n_G = count_mtrx[ci][2] + count_mtrx[ci][6] + count_mtrx[ci][10] + count_mtrx[ci][14];
        n_T = count_mtrx[ci][3] + count_mtrx[ci][7] + count_mtrx[ci][11] + count_mtrx[ci][15];
        if ( (n_A == 0) || (n_C == 0) || (n_G == 0) || (n_T == 0) ) {
            continue;
        }
        rate_mtrx[i][0] = count_mtrx[ci][1] / n_C;
        rate_mtrx[i][1] = count_mtrx[ci][2] / n_G;
        rate_mtrx[i][2] = count_mtrx[ci][3] / n_T;
        rate_mtrx[i][3] = count_mtrx[ci][4] / n_A;
        rate_mtrx[i][4] = count_mtrx[ci][6] / n_G;
        rate_mtrx[i][5] = count_mtrx[ci][7] / n_T;
        rate_mtrx[i][6] = count_mtrx[ci][8] / n_A;
        rate_mtrx[i][7] = count_mtrx[ci][9] / n_C;
        rate_mtrx[i][8] = count_mtrx[ci][11] / n_T;
        rate_mtrx[i][9] = count_mtrx[ci][12] / n_A;
        rate_mtrx[i][10] = count_mtrx[ci][13] / n_C;
        rate_mtrx[i][11] = count_mtrx[ci][14] / n_G;
    }
    return;
}


/* Print count matrices to .pss.counts.txt
   Args: const char* fasta_fn - pointer to FASTA filename
         const char* bam_fn - pointer to BAM filename
         const char* out_prefix - pointer to output filename prefix
         unsigned long** fwd_counts - pointer to forward count matrix 
         unsigned long** rev_counts - pointer to reverse count matrix */
int print_counts(const char* fasta_fn, const char* bam_fn, const char* out_prefix, \
                  unsigned long** fwd_counts, unsigned long** rev_counts) {

    char out_fn[MAX_FN_LEN];
    snprintf(out_fn, sizeof(out_fn), "%s.pss.counts.txt", out_prefix);
    FILE *fp = fopen(out_fn, "w");
    if (!fp) {
        fprintf( stderr, "ERROR: Cannot write to file %s\n.", out_fn);
        return 1;
    }

    fprintf( fp, "### pss-bam.c v1.2.1:\n### FASTA: %s\n### BAM: %s\n### OUT: %s\n", fasta_fn, bam_fn, out_fn );
    fprintf( fp, "### Format of table:\n" );
    fprintf( fp, "### Counts of how often a read base and genome base were seen at\n" );
    fprintf( fp, "### each position in the aligned reads.\n" );
    fprintf( fp, "### First base is what was seen in the read.\n" );
    fprintf( fp, "### Second base is what was in the genome at that position.\n" );
    fprintf( fp, "### POS AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT\n" );
    fprintf( fp, "### Forward read substitution counts and base context\n" );

    for (int i = -2; i < REGION_LEN; i++) {
        fprintf( fp, "%d\t", i );
        for (int j = 0; j < 16; j++) {
            fprintf( fp, "%lu\t", fwd_counts[i+2][j] );
        }
        fprintf( fp, "\n" );
    }

    fprintf( fp, "\n\n### Reverse read substitution counts and base context\n" );
    // print all rows except context rows in rev_counts in reverse order
    for (int i = REGION_LEN-1; i >= 0; i--) {
        fprintf( fp, "%d\t", i);
        for (int j = 0; j < 16; j++) {
            fprintf( fp, "%lu\t", rev_counts[i+2][j] );
        }
        fprintf( fp, "\n" );
    }
    
    // print context rows in rev_counts
    for (int i = 1; i < 3; i++) {
        fprintf( fp, "%d\t", i );
        for (int j = 0; j < 16; j++) {
            fprintf( fp, "%lu\t", rev_counts[2-i][j] );
        }
        fprintf( fp, "\n" );
    }
    fclose(fp);
    return 0;
}


/* Print rate matrices to .pss.rates.txt
   Args: const char* fasta_fn - pointer to FASTA filename
         const char* bam_fn - pointer to BAM filename
         const char* out_prefix - pointer to output filename prefix
         unsigned long** fwd_rates - pointer to forward rate matrix 
         unsigned long** rev_rates - pointer to reverse rate matrix */
int print_rates(const char* fasta_fn, const char* bam_fn, const char* out_prefix, \
                double** fwd_rates, double** rev_rates) {

    char out_fn[MAX_FN_LEN];
    snprintf(out_fn, sizeof(out_fn), "%s.pss.rates.txt", out_prefix);
    FILE *fp = fopen(out_fn, "w");
    if (!fp) {
        fprintf( stderr, "ERROR: Cannot write to file %s\n.", out_fn);
        return 1;
    }

    fprintf( fp, "### pss-bam.c v1.2\n### FASTA: %s\n### BAM: %s\n### OUT: %s\n", fasta_fn, bam_fn, out_fn );
    fprintf( fp, "### Format of table:\n" );
    fprintf( fp, "### Substitution rates for all possible nucleotide substitutions at\n" );
    fprintf( fp, "### each position in the aligned reads.\n" );
    fprintf( fp, "### First base is what was seen in the read.\n" );
    fprintf( fp, "### Second base is what was in the genome at that position.\n" );
    fprintf( fp, "### POS AC AG AT CA CG CT GA GC GT TA TC TG\n" );
    fprintf( fp, "### Forward read substitution rates\n" );

    for (int i = 0; i < REGION_LEN; i++) {
        fprintf( fp, "%d\t", i );
        for (int j = 0; j < 12; j++) {
            fprintf( fp, "%.5e\t", fwd_rates[i][j] );
        }
        fprintf( fp, "\n" );
    }
    
    fprintf( fp, "\n\n### Reverse read substitution rates\n" );
    for (int i = REGION_LEN-1; i >= 0; i--) {
        fprintf( fp, "%d\t", i);
        for (int j = 0; j < 12; j++) {
            fprintf( fp, "%.5e\t", rev_rates[i][j] );
        }
        fprintf( fp, "\n" );
    }
    fclose(fp);
    return 0;
}


/* Free memory allocated by init_count_mtrx() & init_rate_mtrx() */
int destroy_mtrx(void** mtrx, size_t n_rows) {
    if (!mtrx) {
        return 0;
    }
    for (int i = 0; i < n_rows; i++) {
        free(mtrx[i]);
    }
    free(mtrx);
    return 0;
}


/*** MAIN PROGRAM ***/
int main(int argc, char* argv[]) {

    int option;
    char* opts = ":F:B:o:R:r:l:L:q:U:D:m";
    char* fasta_fn, *bam_fn, *out_prefix, *ptr1, *ptr2;
    char* read_group = NULL;
    while ( (option = getopt(argc, argv, opts)) != -1 ) {
        switch (option) {
            case 'F':
                fasta_fn = strdup(optarg);
                break;
            case 'B':
                bam_fn = strdup(optarg);
                break;
            case 'o':
                out_prefix = strdup(optarg);
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
                UP_CTX = optarg;
                break;
            case 'D':
                DOWN_CTX = optarg;
                break;
            case 'm':
                MERGED_ONLY = 1;
                break;
	        case 'R' :
	            read_group = strdup(optarg);
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

    if ( !fasta_fn || !bam_fn || !out_prefix ) {
        fprintf(stderr, "pss-bam v1.2.1: Program for describing base context and counting\n");
        fprintf(stderr, "the number of matches/mismatches in aligned reads to a genome.\n" );
        fprintf(stderr, "-F <reference FASTA (required)>\n");
        fprintf(stderr, "-B <input BAM (required)>\n");
        fprintf(stderr, "-o <output filename prefix (required)>\n");
        fprintf(stderr, "-r <length in basepairs into the interior of alignments to report on (default: 15)>\n");
        fprintf(stderr, "-l <minimum length of read to report (default: 0)>\n");
        fprintf(stderr, "-L <maximum length of read to report (default: 250000000)>\n");
        fprintf(stderr, "-q <map quality filter of read to report (default: 0)>\n" );
        fprintf(stderr, "-R <read group name to restrict analysis to (default: all reads)>\n" );
        fprintf(stderr, "-U <upstream context base filter; first base before alignment must be one of these (default: ACGT)>\n");
        fprintf(stderr, "-D <downstream context base filter; first base before alignment must be one of these (default: ACGT)>\n");
        fprintf(stderr, "-m <only consider merged reads>\n");
        exit(1);
    }

    if (MERGED_ONLY && read_group) {
        fprintf( stderr,
                 "Full command: %s -F %s -B %s -o %s -r %d -l %lu -L %lu -q %d -R %s -U %s -D %s -m\n", 
                  argv[0], fasta_fn, bam_fn, out_prefix, REGION_LEN, MIN_READ_LEN, MAX_READ_LEN, 
                  MIN_MQ, read_group, UP_CTX, DOWN_CTX);
    }
    else if (MERGED_ONLY) {
        fprintf( stderr,
                 "Full command: %s -F %s -B %s -o %s -r %d -l %lu -L %lu -q %d -U %s -D %s -m\n", 
                 argv[0], fasta_fn, bam_fn, out_prefix, REGION_LEN, MIN_READ_LEN, MAX_READ_LEN, MIN_MQ, UP_CTX, DOWN_CTX);
    }
    else if (read_group) {
        fprintf( stderr,
                 "Full command: %s -F %s -B %s -o %s -r %d -l %lu -L %lu -q %d -R %s -U %s -D %s\n", 
                 argv[0], fasta_fn, bam_fn, out_prefix, REGION_LEN, MIN_READ_LEN, MAX_READ_LEN, 
                 MIN_MQ, read_group, UP_CTX, DOWN_CTX);
    }
    else {
        fprintf( stderr,
                 "Full command: %s -F %s -B %s -o %s -r %d -l %lu -L %lu -q %d -U %s -D %s\n", 
                 argv[0], fasta_fn, bam_fn, out_prefix, REGION_LEN, MIN_READ_LEN, MAX_READ_LEN, MIN_MQ, UP_CTX, DOWN_CTX);
    }

    fprintf( stderr, "Reading genome sequence from:\n%s\n", fasta_fn);
    Genome* genome = init_genome(fasta_fn);
    fprintf( stderr, "Finished loading genome.\nCounting matches/mismatches from:\n%s\n", bam_fn );
    
    unsigned long** fwd_counts = init_count_mtrx();
    unsigned long** rev_counts = init_count_mtrx();
    double** fwd_rates = init_rate_mtrx();
    double** rev_rates = init_rate_mtrx();

    FILE* sam_out = bam_to_sam(bam_fn, read_group);
    char saml_buf[MAX_LINE_LEN + 1];
    Saml* sp = malloc(sizeof(Saml));
     
    while ( fgets(saml_buf, MAX_LINE_LEN + 1, sam_out) ) {
        int parse_status = line2saml(saml_buf, sp);
        if (parse_status && DEBUG) {
            fprintf( stderr, "Problem parsing alignment, continuing to next entry...\n" );
            continue;
        }

        int process_status = process_aln(fwd_counts, rev_counts, genome, sp);
        if (process_status == 0) {
            continue;
        }
        else if ( (process_status == 1) && DEBUG ) {
            fprintf( stderr, "%s: Unable to find sequence %s in genome.\n", sp->qname, sp->rname );
        }
        else if ( (process_status == -1) && DEBUG ) {
            fprintf( stderr, "%s: Alignment did not pass filters.\n", sp->qname );
        }
    }
    find_sub_rates(fwd_counts, fwd_rates);
    find_sub_rates(rev_counts, rev_rates);

    print_counts(fasta_fn, bam_fn, out_prefix, fwd_counts, rev_counts);
    print_rates(fasta_fn, bam_fn, out_prefix, fwd_rates, rev_rates);
    
    free(fasta_fn);
    free(bam_fn);
    free(out_prefix);
    free(read_group);
    free(sp);
    fclose(sam_out);
    
    destroy_genome(genome);
    destroy_mtrx( (void**)fwd_counts, REGION_LEN+2 );
    destroy_mtrx( (void**)rev_counts, REGION_LEN+2 );
    destroy_mtrx( (void**)fwd_rates, REGION_LEN );
    destroy_mtrx( (void**)rev_rates, REGION_LEN );

    fprintf(stderr, "Done.\n");
    return 0;
}
