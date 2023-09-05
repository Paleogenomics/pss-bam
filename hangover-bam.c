#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sam-parse.h"

static int MAX_OH_LEN = 20;


FILE* bam_to_sorted_sam(const char* bam_fn) {
    char cmd_buf[MAX_LINE_LEN + 1];
    sprintf(cmd_buf, "samtools sort --output-fmt sam %s", bam_fn);
    FILE* sam_out = popen(cmd_buf, "r"); // get file pointer to sam stdout
    if (sam_out == NULL) {
        fprintf( stderr, "Error: Unable to run 'samtools sort' on %s.\n", bam_fn );
        exit(1);
    }
    return sam_out;
}


size_t get_right_pos(size_t left_pos, const char* cigar) {
    size_t aln_len = 0;
    int len;
    int cigar_len, idx = 0;
    cigar_len = strlen(cigar);
    char code;

    while (idx < cigar_len) {
    sscanf( &cigar[idx], "%d%c", &len, &code );
    if ( code == 'M' || code == 'I' ) {
        aln_len += len;
    }
    else if ( code == 'D' ) {
        aln_len -= len;
    }
    while (cigar[idx] != code) {
        idx++;
    }
    idx++;
  }
  return left_pos + aln_len;
}


long int get_LC_tag(Saml* sp) {
  char* token;
  long int pos;
  char buf[MAX_FIELD_WIDTH + 1];
  strcpy(buf, sp->tags);
  token = strtok(buf, "\t");
  while (token != NULL) {
    if ( (sscanf(token, "LC:i:%ld\n", &pos) == 1) ) {
      return pos;
    }
    else {
      token = strtok(NULL, "\t");
    }
  }
  return -1;
}


int main(int argc, char* argv[]) {

    int option;
    char* opts = ":i:o:m:";
    char in_fn[MAX_FN_LEN];
    char out_prefix[MAX_FN_LEN];
    memset(in_fn, '\0', MAX_FN_LEN);
    memset(out_prefix, '\0', MAX_FN_LEN);

    while ( (option = getopt(argc, argv, opts)) != -1 ) {
        switch (option) {
            case 'i':
                strcpy(in_fn, optarg);
                break;
            case 'o':
                strcpy(out_prefix, optarg);
                break;
            case 'm':
                MAX_OH_LEN = atoi(optarg);
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

    if (in_fn[0] == '\0' || out_prefix[0] == '\0') {
        fprintf( stderr, "hangover-bam: Program for characterizing 5' & 3' overhangs in BAM reads.\n" );
        fprintf( stderr, "The output is a BAM file consisting of read pairs with the overhang length\n" );
        fprintf( stderr, "recorded in the OH tag on the forward read of each pair. These read are\n" );
        fprintf( stderr, "potentially from the opposite strands of the same DNA template.\n" );
        fprintf( stderr, "-i <input BAM file, can be unsorted (required)>\n" );
        fprintf( stderr, "-o <output BAM prefix (required)>\n" );
        fprintf( stderr, "-m <max overhang length (default: 20)>\n" );
        exit(1);
    }


    if (MAX_OH_LEN < 0) {
        fprintf( stderr, "Error: Maximum overhang length must be >= 0.\n");
        exit(1);
    }

    FILE* in_fp = bam_to_sorted_sam(in_fn);
    // temp files for storing tagged (tmp1) &
    // untagged (tmp2) reads after 1st pass
    char tmp1_fn[] = "temp-XXXXXX";
    int fd1 = mkstemp(tmp1_fn);
    char tmp2_fn[] = "temp-XXXXXX";
    int fd2 = mkstemp(tmp2_fn);
    if (fd1 == -1 || fd2 == -1) {
        fprintf( stderr, "Error creating temp files. \n" );
        exit(1);
    }
    FILE* tmp1_fp = fdopen(fd1, "w");
    FILE* tmp2_fp = fdopen(fd2, "w");
    Saml* cur_sp = malloc(sizeof(Saml)); //pointer to current read
    Saml* next_sp = malloc(sizeof(Saml)); // pointer to next read
    char buf[MAX_LINE_LEN + 1];
    char cur_buf[MAX_LINE_LEN + 1];
    char* lp1 = fgets(buf, MAX_LINE_LEN + 1, in_fp);
    size_t cur_rc; //right coordinate of current read

    /* 1ST PASS STARTS */
    // skip SAM header
    while ( lp1[0] == '@' ) {
        fprintf(tmp1_fp, "%s", buf);
        fprintf(tmp2_fp, "%s", buf);
        lp1 = fgets(buf, MAX_LINE_LEN + 1, in_fp);
    }
    while ( lp1 ) {
        strcpy(cur_buf, buf);
        int res1 = line2saml(cur_buf, cur_sp);
        if (res1) {
            lp1 = fgets(buf, MAX_LINE_LEN + 1, in_fp);
            continue;
        }

        char* lp2 = fgets(buf, MAX_LINE_LEN + 1, in_fp);
        int res2 = line2saml(buf, next_sp);
        while ( lp2 && res2 ) {
            lp2 = fgets(buf, MAX_LINE_LEN + 1, in_fp);
            res2 = line2saml(buf, next_sp);
        }
        // skip this current read if already at end
        // of chromosome
        if (strcmp(next_sp->rname, cur_sp->rname) != 0) {
            lp1 = lp2;
            continue;
        }
        // overhang is always (+) because next read's pos
        // is always larger than currrent read's pos;
        // we only tag forward reads
        int oh = next_sp->pos - cur_sp->pos;
        if ( oh <= MAX_OH_LEN &&
             cur_sp->reverse &&
             !next_sp->reverse ) {
            fprintf( tmp1_fp, "%s", cur_buf );
            buf[strlen(buf)-1] = '\t';
            // so set OH to (-) if current read is reverse (3' overhang)
            fprintf( tmp1_fp, "%sOH:i:%d\n", buf, -1 * oh );
            lp1 = fgets(buf, MAX_LINE_LEN + 1, in_fp);
            continue;
        }

        else if ( oh <= MAX_OH_LEN &&
                  !cur_sp->reverse &&
                  next_sp->reverse ) {
            cur_buf[strlen(cur_buf)-1] = '\t';
            // OH stays (+) if current read is forward (5' overhang)
            fprintf( tmp1_fp, "%sOH:i:%d\n", cur_buf, oh );
            fprintf( tmp1_fp, "%s", buf );
            lp1 = fgets(buf, MAX_LINE_LEN + 1, in_fp);
            continue;
        }
        // if current read is not tagged, write it to the temp file,
        // replacing its POS field with the right coord and saving
        // the left coord in the LC tag
        cur_rc = get_right_pos(cur_sp->pos, cur_sp->cigar);
        if (cur_sp->qual[strlen(cur_sp->qual)-1] == '\n') {
            cur_sp->qual[strlen(cur_sp->qual)-1] = '\0'; //just to make sure no new line char here
        }
        if (strlen(cur_sp->tags) > 0) {
            cur_sp->tags[strlen(cur_sp->tags)-1] = '\0'; //same here
        }
        fprintf( tmp2_fp, "%s\t%u\t%s\t%lu\t%u\t%s\t%s\t%u\t%i\t%s\t%s\t%s\tLC:i:%lu\n",
                 cur_sp->qname,
                 cur_sp->flag,
	             cur_sp->rname,
	             cur_rc,
	             cur_sp->mapq,
	             cur_sp->cigar,
	             cur_sp->mrnm,
	             cur_sp->mpos,
	             cur_sp->isize,
	             cur_sp->seq,
	             cur_sp->qual,
	             cur_sp->tags,
                 cur_sp->pos );

        lp1 = lp2; // then set the next read as current read
    }
    fflush(tmp2_fp);

    /** 2ND PASS STARTS **/
    in_fp = bam_to_sorted_sam(tmp2_fn);
    lp1 = fgets(buf, MAX_LINE_LEN + 1, in_fp);

    // skip SAM header
    while ( lp1[0] == '@' ) {
        lp1 = fgets(buf, MAX_LINE_LEN + 1, in_fp);
    }
    while ( lp1 ) {
        strcpy(cur_buf, buf);
        int res1 = line2saml(cur_buf, cur_sp);
        if (res1) {
            lp1 = fgets(buf, MAX_LINE_LEN + 1, in_fp);
            continue;
        }

        char* lp2 = fgets(buf, MAX_LINE_LEN + 1, in_fp);
        int res2 = line2saml(buf, next_sp);
        while ( lp2 && res2 ) {
            lp2 = fgets(buf, MAX_LINE_LEN + 1, in_fp);
            res2 = line2saml(buf, next_sp);
        }
        // skip this current read if already at end
        // of chromosome
        if (strcmp(next_sp->rname, cur_sp->rname) != 0) {
            lp1 = lp2;
            continue;
        }

        int oh = next_sp->pos - cur_sp->pos;
        long int cur_lc, next_lc;
        // 3' overhang
        if ( oh <= MAX_OH_LEN &&
             cur_sp->reverse &&
             !next_sp->reverse ) {
            cur_lc = get_LC_tag(cur_sp);
            next_lc = get_LC_tag(next_sp);

            // if LC tag parsed ok for current & next reads,
            // write sam lines to output, replacing the POS
            // field with the actual left coords
            if (cur_lc != -1 && next_lc != -1) {
                fprintf( tmp1_fp, "%s\t%u\t%s\t%lu\t%u\t%s\t%s\t%u\t%i\t%s\t%s\t%s",
                         cur_sp->qname,
                         cur_sp->flag,
                         cur_sp->rname,
                         cur_lc,
                         cur_sp->mapq,
	                     cur_sp->cigar,
                         cur_sp->mrnm,
                         cur_sp->mpos,
                         cur_sp->isize,
                         cur_sp->seq,
                         cur_sp->qual,
                         cur_sp->tags );
                
                // next read is forward -> tag it
                next_sp->tags[strlen(next_sp->tags)-1] = '\0'; // remove new line char
                fprintf( tmp1_fp, "%s\t%u\t%s\t%lu\t%u\t%s\t%s\t%u\t%i\t%s\t%s\t%s\tOH:i:%d\n",
                         next_sp->qname,
                         next_sp->flag,
                         next_sp->rname,
                         next_lc,
                         next_sp->mapq,
	                     next_sp->cigar,
                         next_sp->mrnm,
                         next_sp->mpos,
                         next_sp->isize,
                         next_sp->seq,
                         next_sp->qual,
                         next_sp->tags,
                         -1 * oh );
                lp1 = fgets(buf, MAX_LINE_LEN + 1, in_fp);
                continue;
            }
        }

        // 5' overhang
        else if ( oh <= MAX_OH_LEN &&
                  !cur_sp->reverse &&
                  next_sp->reverse ) {
            cur_lc = get_LC_tag(cur_sp);
            next_lc = get_LC_tag(next_sp);
            if (cur_lc != -1 && next_lc != -1) {
                cur_sp->tags[strlen(cur_sp->tags)-1] = '\0'; //remove new line char
                //current read is forward -> tag it
                fprintf( tmp1_fp, "%s\t%u\t%s\t%lu\t%u\t%s\t%s\t%u\t%i\t%s\t%s\t%s\tOH:i:%d\n",
                         cur_sp->qname,
                         cur_sp->flag,
                         cur_sp->rname,
                         cur_lc,
                         cur_sp->mapq,
	                     cur_sp->cigar,
                         cur_sp->mrnm,
                         cur_sp->mpos,
                         cur_sp->isize,
                         cur_sp->seq,
                         cur_sp->qual,
                         cur_sp->tags,
                         oh );

                fprintf( tmp1_fp, "%s\t%u\t%s\t%lu\t%u\t%s\t%s\t%u\t%i\t%s\t%s\t%s",
                         next_sp->qname,
                         next_sp->flag,
                         next_sp->rname,
                         next_lc,
                         next_sp->mapq,
	                     next_sp->cigar,
                         next_sp->mrnm,
                         next_sp->mpos,
                         next_sp->isize,
                         next_sp->seq,
                         next_sp->qual,
                         next_sp->tags );
                lp1 = fgets(buf, MAX_LINE_LEN + 1, in_fp);
                continue;
            }
        }
        lp1 = lp2;
    }
    fflush(tmp1_fp); // write remaining buffered data to tmp1

    // now convert back to BAM, removing the LC tag
    char cmd_buf[MAX_LINE_LEN + 1];
    sprintf(cmd_buf, "samtools view -x LC -b %s -o %s.bam ", tmp1_fn, out_prefix);
    system(cmd_buf);
    unlink(tmp1_fn);
    unlink(tmp2_fn);
    free(cur_sp);
    free(next_sp);
    fclose(in_fp);
    return 0;
}