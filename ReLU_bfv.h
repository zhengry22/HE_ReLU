#pragma once
#include "../src/troy.h"
#include "../Polynomial_Calc/polynomial.h"
#include "../Polynomial_Calc/SiLU.h"
using namespace troy;
using namespace std;
/*
Since the ReLU function doesn't satisty the homomorphic rule, 
we may want to simulate the function using a polynomial.
*/

/*
Follow are the parameters
*/

// These are the adjustable variables
const size_t poly_modulus_degree = 4096;
const size_t plain_modulus = 8192;


void poly_relu(const auto &encoder, const Evaluator &evaluator, const Ciphertext &x_encrypted, Ciphertext &encrypted_result);
void horner(const auto &encoder, const Evaluator &evaluator, const RelinKeys &relin_keys, const EncryptPolynomial &encpoly, const Ciphertext &x_encrypted, Ciphertext &encrypted_result);
uint64_t ReLU();