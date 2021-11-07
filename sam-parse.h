#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>

#define MAX_LINE_LEN (200000)
#define MAX_FN_LEN (2047)
#define MAX_FIELD_WIDTH (2047)
#define MAX_ID_LEN (512)
#define MATCH (1)
#define MISMATCH (4)
#define GAP_OPEN (6)
#define GAP_EXT (1)
/* (c) 2019 Astrea Forensics
   Ed Green
   Parser for sam lines */


typedef struct saml {
  char qname[MAX_FIELD_WIDTH + 1];
  unsigned int flag;
  char rname[ MAX_FIELD_WIDTH + 1];
  unsigned int pos;
  unsigned int mapq;
  char cigar[MAX_FIELD_WIDTH + 1];
  char mrnm[MAX_FIELD_WIDTH + 1];
  unsigned int mpos;
  int isize;
  int seq_len;
  char seq[MAX_FIELD_WIDTH + 1];
  char qual[MAX_FIELD_WIDTH + 1];
  char tags[MAX_FIELD_WIDTH + 1];
  char BC[MAX_FIELD_WIDTH + 1];
  char RG[MAX_FIELD_WIDTH + 1];
  int aln_seq_len;
  int NM;
  int AS; // alignment score
  int XM; // alignment mismatches
  int XO; // alignment gap opens
  int XG; // alignment gap extends
  // could add others
} Saml;

int line2saml( const char* line, Saml* sp );
int is_header( const char* line );
int aln_seq_len( const char* cigar );
int good_score( Saml* sp, float m, float b );
