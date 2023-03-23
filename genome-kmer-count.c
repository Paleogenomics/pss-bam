#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include "fasta-genome-io.h"
#include "kmer.h"
#define K_DEF (4)

void count_kmers( const Seq* seq, KSP kcount );
int inx2kmer( size_t inx, char* kmer, size_t k );

void help( void ) {
  printf( "genome-kmer-count -f <fasta genome file>\n" );
  printf( "                  -k <kmer size; default = %u>\n",K_DEF );
  printf( "This program reports the number of observed number\n");
  printf( "of all possible kmers of the given length in the\n");
  printf( "input genome.\n" );
  exit( 0 );
}


int main( int argc, char* argv[] ) {
  extern char* optarg;
  int ich;
  size_t k = K_DEF;
  size_t i, max_k_inx, k_inx;
  char fa_in[ MAX_FN_LEN + 1 ]     = {'\0'};
  Genome* genome;
  Seq* seq;
  Fa_Src* fa_src;
  KSP kcount;
  char* kmer;
  
  while( (ich=getopt( argc, argv, "f:k:" ) ) != -1 ) {
    switch(ich) {
    case 'f' :
      strcpy( fa_in, optarg );
      break;
    case 'k' :
      k = (size_t)atoi( optarg );
      break;
    default :
      help();
    }
  }
  if ( strlen( fa_in ) == 0 ) {
    help();
  }

  genome = init_genome(fa_in);
  printf( "Parsed input genome. Found %lu sequences.\n",
	  genome->n_seqs );
  kcount = init_KSP(k);

  for (i = 0; i < genome->n_seqs; i++) {
    count_kmers( genome->seqs[i], kcount );
  }
  max_k_inx = 1<<(2*k);
  kmer = (char*)malloc(sizeof(char) * (k + 1));
  for (k_inx = 0; k_inx < max_k_inx; k_inx++) {
    inx2kmer( k_inx, kmer, k );
    printf( "%s\t%u\n", kmer, kmer2count( kmer, kcount ) );
  }
  exit( 0 );
}

void count_kmers( const Seq* seq, KSP kcount ) {
  size_t i;
  size_t num_kmer_pos;
  num_kmer_pos = seq->len - kcount->k + 1;
  for ( i = 0; i < num_kmer_pos ; i++ ) {
    if ( add_to_ksp( &seq->seq[i], kcount ) ) {
      //      fprintf( stderr, "Problem adding kmer to count\n" );
      ;
    }
  }
  return;
}


/* Converts an unsigned integer to a kmer using the
   conversion formula: 00=A, 01=C, 10=G, 11=T
   of requested length */
int inx2kmer( size_t inx, char* kmer, size_t k ) {
  // inx is bit string, 2 bits per base
  // kmer is pointer to make kmer in bases
  // k is length of kmer
  kmer[k] = '\0';
  size_t mask = 3; // 000011
  size_t base_code;
  char b;
  size_t i;
  for ( i = 0; i < k; i++ ) {
    base_code = (inx & mask);
    switch (base_code) {
    case 0 :
      kmer[k-i-1] = 'A';
      break;
    case 1 :
      kmer[k-i-1] = 'C';
      break;
    case 2 :
      kmer[k-i-1] = 'G';
      break;
    case 3 :
      kmer[k-i-1] = 'T';
      break;
    default :
      fprintf( stderr, "Can't decode %lu\n", base_code );
    }
    inx = inx >> 2;
  }
  return 0;
}
