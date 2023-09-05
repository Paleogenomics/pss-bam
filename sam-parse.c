#include "sam-parse.h"

/* (c) 2019 Astrea Forensics
   Ed Green
   Parser for sam lines */

/* Return 0 => Saml* sp is populations; copacetic
          1 => problem
*/
int line2saml( const char* line, Saml* sp ) {
  /* Note that not all TAGs are captured in the struct saml,
     only the ones I care about */
  /* 1 0x1 template having multiple segments in sequencing
     2 0x2 each segment properly aligned according to the aligner
     4 0x4 segment unmapped
     8 0x8 next segment in the template unmapped
     16 0x10 SEQ being reverse complemented
     32 0x20 SEQ of the next segment in 
             the template being reverse complemented
     64 0x40 the first segment in the template
     128 0x80 the last segment in the template
     256 0x100 secondary alignment
     512 0x200 not passing filters, such as 
               platform/vendor quality controls
     1024 0x400 PCR or optical duplicate
     2048 0x800 supplementary alignment */
  size_t pos, len;
  char tag[ 3 ];
  char type;
  int field = 0;
  int end_of_samline = 0; // set to true once we've parsed it all
  int got_AS = 0;
  int got_XM = 0;
  int got_XO = 0;
  int got_XG = 0;
  char value[ MAX_FIELD_WIDTH ];
  if ( sscanf( line,
	       "%s\t%u\t%s\t%lu\t%u\t%s\t%s\t%u\t%i\t%s\t%s\t%s",
	       sp->qname,
	       &sp->flag,
	       sp->rname,
	       &sp->pos,
	       &sp->mapq,
	       sp->cigar,
	       sp->mrnm,
	       &sp->mpos,
	       &sp->isize,
	       sp->seq,
	       sp->qual,
	       sp->tags ) >= 11 ) {

    if ( strlen( sp->seq ) == strlen( sp->qual ) ) {
      sp->seq_len = strlen( sp->seq );

      sp->paired = (sp->flag >> 0) & 1;
      sp->proper_pair = (sp->flag >> 1) & 1;
      sp->unmap = (sp->flag >> 2) & 1;
      sp->munmap = (sp->flag >> 3) & 1;
      sp->reverse = (sp->flag >> 4) & 1;
      sp->mreverse = (sp->flag >> 5) & 1;
      sp->read1 = (sp->flag >> 6) & 1;
      sp->read2 = (sp->flag >> 7) & 1;
      sp->secondary = (sp->flag >> 8) & 1;
      sp->qc_failed = (sp->flag >> 9) & 1;
      sp->duplicate = (sp->flag >> 10) & 1;
      sp->supplementary = (sp->flag >> 11) & 1;
      
      if ( sp->paired == 0 ) {
        sp->isize = strlen( sp->seq );
      }
      
      len = strlen( line );
      pos = 0;
      /* Advance pos to first position past mandatory fields
	 This is either the end of the line or the beginning
	 of the optional field */
      while( field < 11 ) {
        if ( line[pos] == '\t' ) {
          field++;
        }
        pos++;
      }

      if (pos < len) {
        strcpy(sp->tags, &line[pos]);
      }
      return 0;
    }
    return 1;
  }
  return 1;
}

/* Takes CIGAR string and returns the number of
   aligned bases from this read in the alignment 
   Args: const char* cigar: the cigar string
   Returns: int: the number of aligned bases
   Note: The aligned bases may be matches or mismatches
         The number returned is just the number of 
	 bases involved in the alignment
*/
int aln_seq_len( const char* cigar ) {
  int pos = 0;
  int aln_seq_len = 0;
  int cigar_len;
  char token;
  int len;

  cigar_len = strlen( cigar ); // length of cigar string
  while( pos < cigar_len ) {
    // get the next length and token
    sscanf( &cigar[pos], "%d%c", &len, &token );
    if ( token == 'M' ) {
      // if the token is M, add this length to the
      // aln_seq_len
      aln_seq_len += len;
    }
    while( cigar[pos] != token ) {
      // wind through to the token you just got
      pos++;
    }
    // and then one past it
    pos++;
  }
  return aln_seq_len;
}

/* Takes a sam line char* as argument
   Returns true if the line begins with the @
   character that signifies a header line */
int is_header( const char* line ) {
  if ( line[0] == '@' ) {
    return 1;
  }
  else {
    return 0;
  }
}

/* Args: (1) pointer to struct saml with alignment info
         (2) float m
	 (3) float b
   Returns: 1 (TRUE) if sp->AS >=  m * strlen(sp->seq) + b
                  or if there is no sp->AS
            0 (FALSE) iff sp->AS < m * strlen(sp->seq) + b

   NOTE: The score filtering is based on the length of the
         sequence. Some of it may not be aligned. There is no
	 examination of the cigar string to look for hard or 
	 soft masking.
         This may not be what you want!
*/ 
int good_score( Saml* sp, float m, float b ) {
  if ( sp->AS > 0 ) {
    if ((float)sp->AS >= (m * sp->seq_len) + b ) {
      return 1;
    }
    else {
      return 0;
    }
  }
  return 1;
}


