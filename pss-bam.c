#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <regex.h>

#include "fasta-genome-io.h"
#include "sam-parse.h"

#define TRUE 1
#define FALSE 0

static unsigned int REGION_LEN = 15;
static unsigned int MIN_READ_LEN = 0;
static unsigned long MAX_READ_LEN = 250000000;
static unsigned int MIN_MQ = 0;
static char UPSTR_BASE_CNTXT[5] = "ACTG";
static char DWNSTR_BASE_CNTXT[5] = "ACTG";
static int MERGED_ONLY = 0;


unsigned long** init_count_table() {
    /* extra 2 bases upstream/downstream of aln start */
    unsigned long** count_tab = malloc((REGION_LEN+2) * sizeof(unsigned long*));
    for (int i = 0; i < REGION_LEN+2; i++) {
        unsigned long* count_arr = malloc(16 * sizeof(unsigned long)); //16 total base substitutions
        for (int j = 0; j < 16; j++) {
             count_arr[j] = 0;
        }
        count_tab[i] = count_arr;
    }
    return count_tab;
}

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

int read_len_ok(int seq_len) {
    if ( (seq_len >= MIN_READ_LEN) && (seq_len <= MAX_READ_LEN) ) {
        return 1;
    }
    return 0;
}


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
    also check if 2 extra bases are present upstream/downstream of aln */
int context_base_ok(char fb_upstr, char sb_upstr, \
                    char fb_dwnstr, char sb_dwnstr) {
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

int context_count(unsigned long** count_table, const char fcb, const char scb) {
    char context_bases[3] = {scb, fcb}; //second context base is 1st row
    for (int i = 0; i < 2; i++) {
        switch (context_bases[i]) {
            case 'A':
                count_table[i][0] += 1;
                break;
            case 'C':
                count_table[i][5] += 1;
                break;
            case 'G':
                count_table[i][10] += 1;
                break;
            case 'T':
                count_table[i][15] += 1;
                break;
            default:
                return 0;
        }
    }
    return 0;
}

int fwd_substitution_count(unsigned long** fwd_counts, const char* ref_seq, \
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
    return 0;
}

int rev_substitution_count(unsigned long** rev_counts, const char* ref_seq, \
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
    return 0;
}


int fill_count_table(const char* fasta_fn, const char* bam_fn) {
    Genome* genome = init_genome(fasta_fn);
    unsigned long** fwd_counts = init_count_table();
    unsigned long** rev_counts = init_count_table();

    char cmd_buf[512];
    sprintf(cmd_buf, "samtools view %s", bam_fn);
    FILE* bp = popen(cmd_buf, "r"); //get file pointer to sam stdout
    if (bp == NULL) {
        fprintf(stderr, "Error: Unable to open %s with samtools view.\n", bam_fn);
        exit(1);
    }

    char saml_buf[MAX_LINE_LEN+1];
    Saml* sp = malloc(sizeof(Saml*));

    while (fgets(saml_buf, MAX_LINE_LEN+1, bp)) {
        int status = line2saml(saml_buf, sp);
        if (status) {
            fprintf(stderr, "Problem parsing alignment entry, continuing to next entry...\n");
            continue;
        }

        Seq* rseq = find_seq(genome, sp->rname);
        char* ref_seq = rseq->seq;

        unsigned long aln_start = (sp->pos)-1; //pos in sam file is 1-based
        unsigned long aln_end = aln_start+(sp->seq_len)-1;

        if ( (sp->mapq >= MIN_MQ) 
            && (read_len_ok(sp->seq_len))
            //&& (cigar_ok(sp->cigar)) 
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
                        int fc_count = context_count(fwd_counts, rfb_upstr, rsb_upstr); //count upstream context bases
                        int rc_count = context_count(rev_counts, rfb_dwnstr, rsb_dwnstr); //count downstream context bases
                        char* rev_read = do_reverse_complement(sp->seq);
                        int fs_count = fwd_substitution_count(fwd_counts, rev_ref, rev_read, aln_start);
                        int rs_count = rev_substitution_count(rev_counts, rev_ref, rev_read, sp->seq_len, aln_end);
                    }
                }

                char fb_upstr = ref_seq[aln_start-1];
                char sb_upstr = ref_seq[aln_start-2];
                char fb_dwnstr = ref_seq[aln_end+1];
                char sb_dwnstr = ref_seq[aln_end+2];

                if (context_base_ok(fb_upstr, sb_upstr, fb_dwnstr, sb_dwnstr)) {
                        int fc_count = context_count(fwd_counts, fb_upstr, sb_upstr);
                        int rc_count = context_count(rev_counts, fb_dwnstr, sb_dwnstr);
                        int fs_count = fwd_substitution_count(fwd_counts, ref_seq, sp->seq, aln_start);
                        int rs_count = rev_substitution_count(rev_counts, ref_seq, sp->seq, sp->seq_len, aln_end);
                }
        }
        else {
            printf("Read %s did not pass filters\n", sp->qname);
        }
    }
    printf("Forward tab:\n");
    for (int i = 0; i < REGION_LEN+2; i++) {
        for (int j = 0; j < 16; j++) {
            printf("%lu\t", fwd_counts[i][j]);
        }
        printf("\n");
    }
    
    return 0;
}

int main() {
    int a = fill_count_table("test.fa", "test.bam");
    return 0;
}

