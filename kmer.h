#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#define K_AR_SIZE (8) // default length of the array size of the kmer structure

typedef struct kmer_tree_node {
  struct kmer_tree_node* Ap;
  struct kmer_tree_node* Cp;
  struct kmer_tree_node* Gp;
  struct kmer_tree_node* Tp;
  unsigned int count;
} ktn;
typedef struct kmer_tree_node* ktnP;

typedef struct kmers {
  size_t k; // length of kmers
  size_t k_ar_size ; // length of kmer part that we'll handle in the array
                 // and not the tree
  ktnP* ka; // the array part;
} Kmers;
typedef struct kmers* KSP;

/* Function prototypes */
KSP init_KSP( int k );
int add_to_ksp( const char* kmer, KSP ks );
unsigned int kmer2count( const char* kmer, const KSP ks );
int kmer2inx( const char* kmer,
	      const size_t kmer_len,
	      size_t* inx );
int destroy_KSP(KSP ks);
