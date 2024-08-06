#include "../src/troy.h"

/*
Since the ReLU function doesn't satisty the homomorphic rule, 
we may want to simulate the function using a polynomial.
*/

/*
Follow are the parameters
*/

// These are the adjustable variables
const size_t poly_modulus_degree = 4096;
const size_t plain_modulus = 1024;


void poly_relu();
uint64_t ReLU();