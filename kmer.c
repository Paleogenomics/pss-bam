#include "kmer.h"

KSP init_KSP( int k ) {
  KSP ks;
  size_t i, len;
  // shifting left by n bits = multiply by 2^n
  len = 1<<(K_AR_SIZE*2); // Length of the array for the array part
  ks = (KSP)malloc(sizeof(Kmers));
  ks->k_ar_size = K_AR_SIZE;
  ks->k = k;
  ks->ka = (ktnP*)malloc(sizeof(ktnP) * len);
  for( i = 0; i < len; i++ ) {
    ks->ka[i] = NULL;
  }
  return ks;
}


static inline ktnP init_ktn( void ) {
  ktnP new_ktn;
  new_ktn = (ktnP)malloc(sizeof(ktn));
  new_ktn->Ap = NULL;
  new_ktn->Cp = NULL;
  new_ktn->Gp = NULL;
  new_ktn->Tp = NULL;
  new_ktn->count = 0;
  return new_ktn;
}


/* This function increments the kmer count of a given kmer in the KSP ks
   Args: const char* kmer - pointer to the kmer
         KSP ks - pointer to the ks
   Returns: 0 - everything worked
           -1 - This is not a valid kmer
   Notes: Because the data structure (ks) has two parts, the function
          first finds the array index of the first ks->k_ar_size bases
	  in the kmer. That array position points to a suffix tree of
	  all kmers that extend from there.
	  If this kmer has not been seen, the function populates the tree
	  to the extent necessary to describe this kmer
*/
int add_to_ksp( const char* kmer, KSP ks ) {
  size_t inx;
  ktnP curr_node, next_node;
  size_t kmer_pos;
  char next_base;
  unsigned int init_k = ks->k;
  if (ks->k_ar_size < ks->k) {
    init_k = ks->k_ar_size;
  }

  /* Find what the index position is for the first K_AR_SIZE
     bases in this kmer */
  if ( kmer2inx( kmer, init_k , &inx ) ) {
    /* Start at the index position of the first ks->k_ar_size bases */
    curr_node = ks->ka[inx];
    
    /* If we've never seen that before, then initialize it */
    if ( curr_node == NULL ) {
      curr_node = init_ktn();
      ks->ka[inx] = curr_node;
    }
    
    kmer_pos = init_k;

    while( kmer_pos < ks->k ) {
      next_base = toupper(kmer[kmer_pos]);
      switch( next_base ) {
        case 'A' :
          if ( curr_node->Ap == NULL ) {
	          curr_node->Ap =  init_ktn();
	        }
	        next_node = curr_node->Ap;
	        break;
        case 'C' :
	        if ( curr_node->Cp == NULL ) {
	          curr_node->Cp = init_ktn();
	        }
	        next_node = curr_node->Cp;
	        break;
        case 'G' :
	        if ( curr_node->Gp == NULL ) {
	          curr_node->Gp = init_ktn();
	        }
	        next_node = curr_node->Gp;
	        break;
        case 'T' :
	        if ( curr_node->Tp == NULL ) {
	          curr_node->Tp = init_ktn();
	        }
	        next_node = curr_node->Tp;
	        break;
        default :
	      // not a good base, not a good kmer, we're done
	        return -1;
      }
      curr_node = next_node;
      kmer_pos++;
    }
  
    if ( curr_node->count < UINT_MAX ) {
      curr_node->count += 1;
    }
    return 0;
  }
  else {
    return -1;
  }
}


/* kmer2count
   Args: (1) a pointer to a character string starting at the kmer of interest
         (2) KSP ks (pointer to kmer structure)
   Returns: short int - count of the number of times we've seen this kmer
   Notes: This is an accesor function that is used to determine how many times
          a given kmer has been counted.
*/
unsigned int kmer2count( const char* kmer, const KSP ks ) {
  unsigned int count;
  size_t inx, kmer_pos;
  ktnP curr_node;
  char curr_base;
  unsigned int init_k = ks->k;
  if (ks->k_ar_size < ks->k) {
    init_k = ks->k_ar_size;
  }
  /* First, find the inx of the array part, i.e., the beginning */
  if ( kmer2inx( kmer, init_k, &inx ) ) {
    curr_node = ks->ka[inx];
    // never saw the kmer before (never init_ktn() at ks->ka[inx])
    if (curr_node == NULL) {
      return 0;
    }
  }
  // something wrong with the input kmer (contains 'N', etc.)
  else { 
    return 0;
  }
  /* Now, follow the pointers until you get to the end of the kmer
     or until we never saw this kmer before */
  for( kmer_pos = init_k; kmer_pos < ks->k; kmer_pos++ ) {
    curr_base = kmer[kmer_pos];
    switch( curr_base ) {
      case 'A' :
        curr_node = curr_node->Ap;
        break;
      case 'C' :
        curr_node = curr_node->Cp;
        break;
      case 'G' :
        curr_node = curr_node->Gp;
        break;
      case 'T' :
        curr_node = curr_node->Tp;
        break;
      default :
        curr_node = NULL;
    }
    if ( curr_node == NULL ) {
      return 0;
    }
  }
  return curr_node->count;
}


/* kmer2inx
   Args: (1) a pointer to a character string;
             the kmer to find the corresponding index of;
	     might not be null-terminated
	 (2) length of the kmer
	 (3) pointer to size_t to put the index
   Returns: TRUE if the index was set, FALSE if it could not
            be set because of some non A,C,G,T character
   Uses the formula A=>00, C=>01, G=>11, T=>11 to make a
   bit string for the kmer. Any other character is not allowed
   and will cause an error
   The bit string is constructed by reading the kmer from left
   to right. This bit-string is then interpreted as a variable
   of type size_t and is appropriate as an array index
*/
int kmer2inx( const char* kmer,
	      const size_t kmer_len,
	      size_t* inx ) {
  size_t l_inx = 0;
  int i = 0;
  char curr_char;

  while( i < kmer_len ) {
    l_inx = l_inx << 2;
    curr_char = toupper(kmer[i]); // Upper case it in case it is not
    switch( curr_char ) {
      case 'A' :
        l_inx += 0;
        break;
      case 'C' :
        l_inx += 1;
        break;
      case 'G' :
        l_inx += 2;
        break;
      case 'T' :
        l_inx += 3;
        break;
      default :
        return 0; // not valid!
    }
    i++;
  }
  *inx = l_inx;
  return 1; // valid!
}


/* destroy_KSP: Free memory allocated with init_KSP()
   Args: KSP ks - pointer to Kmers structure
   Returns: 0 after freeing */
int destroy_KSP(KSP ks) {
  if (!ks) {
    return 0;
  }
  // free 
  for (int i = 0; i < (1<<(ks->k*2)); i++) {
    free(ks->ka[i]);
  }
  free(ks->ka);
  free(ks);
  return 0;
}
