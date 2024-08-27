#pragma once
#include "../src/troy.h"
#include "../Polynomial_Calc/polynomial.h"
#include "../Polynomial_Calc/SiLU.h"
using namespace troy;
using namespace std;

double relu(double x);
void approx_with_fix(const CKKSEncoder &encoder, const Evaluator &evaluator, const RelinKeys &relin_keys, const double &scale,
                    const Ciphertext &x_encrypted, Ciphertext &encrypted_result);
void horner(const auto &encoder, const Evaluator &evaluator, const RelinKeys &relin_keys, const Polynomial<double> poly, const Ciphertext &x_encrypted, Ciphertext &encrypted_result);