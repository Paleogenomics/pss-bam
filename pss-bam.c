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
static char UPSTR_BASE_CONTXT[5] = "ACTG";
static char DWNSTR_BASE_CONTXT[5] = "ACTG";
static int MERGED_ONLY = 0;


unsigned long** init_output_tab() {
    /* extra 2 bases upstream/downstream of alignment start */
    unsigned long** output_tab = malloc((REGION_LEN+2) * sizeof(unsigned long*));
    return output_tab;
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


int fill_output_tab(const char* fasta_fn, const char* bam_fn) {
    Genome* genome = init_genome(fasta_fn);
    
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

        unsigned long pos = sp->pos;
        char* rname = sp->rname;
        Seq* rseq = find_seq(genome, rname);

        if ( (sp->mapq >= MIN_MQ) 
            && (read_len_ok(sp->seq_len))
            && (cigar_ok(sp->cigar)) 
            && (sp->read_paired == FALSE)  
            && (sp->read_unmapped == FALSE) 
            && (sp->not_primary == FALSE) 
            && (sp->not_passing_filters == FALSE)
            && (sp->is_duplicate == FALSE) 
            && (sp->supplementary_alignment == FALSE) ) {
            printf("Read %s passes filters\n", sp->qname);
        }
        else {
            printf("Read %s did not pass filters\n", sp->qname);
        }
        
    }
    return 0;
}

int main() {
    int a = fill_output_tab("test.fa", "test.bam");
    return 0;
}

