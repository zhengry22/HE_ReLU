#pragma once
#include "ReLU_ckks.h"
#include "../examples/examples.h"
#include "../Polynomial_Calc/polynomial.h"
#include "../Polynomial_Calc/SiLU.h"
#include <iostream>

using namespace troy;
using namespace std;

double relu(double x) {
    return x > 0 ? x : 0;
}

/*
    Try to calculate the cipher with a fixed polynomial x^2 + 7x + 10
*/ 
void approx_with_fix(const CKKSEncoder &encoder, const Evaluator &evaluator, const RelinKeys &relin_keys, const double &scale,
                    const Ciphertext &x_encrypted, Ciphertext &encrypted_result) {
    

    // prepare the plaintext for 2 and 5
    Plaintext plain_2, plain_5;
    encoder.encode_float64_single(2.0, std::nullopt, scale, plain_2);
    encoder.encode_float64_single(5.0, std::nullopt, scale, plain_5);


    // x + 2
    Ciphertext x_plus_2 = x_encrypted;
    evaluator.add_plain_inplace(x_plus_2, plain_2);
    //evaluator.relinearize_inplace(x_plus_2, relin_keys);
    //evaluator.rescale_to_next_inplace(x_plus_2);

    // x + 5
    Ciphertext x_plus_5 = x_encrypted;
    evaluator.add_plain_inplace(x_plus_5, plain_5);
    //evaluator.relinearize_inplace(x_plus_5, relin_keys);
    //evaluator.rescale_to_next_inplace(x_plus_5);

    evaluator.multiply_inplace(x_plus_2, x_plus_5);
    evaluator.relinearize_inplace(x_plus_2, relin_keys);
    evaluator.rescale_to_next_inplace(x_plus_2);
    encrypted_result = x_plus_2;

}

int main() {

    EncryptionParameters parms(SchemeType::CKKS);

    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::create(poly_modulus_degree, { 60, 40, 40, 60 }));
    // Change this modulus chain if we need higher degree

    double scale = pow(2.0, 40);

    auto context = HeContext::create(parms, true, SecurityLevel::Classical128);
    print_parameters(*context);
    cout << endl;

    CKKSEncoder encoder(context);
    size_t slot_count = encoder.slot_count();
    cout << "Number of slots: " << slot_count << endl;

    context->to_device_inplace();
    encoder.to_device_inplace();

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key = keygen.create_public_key(false);
    RelinKeys relin_keys = keygen.create_relin_keys(false);
    GaloisKeys gal_keys = keygen.create_galois_keys(false);
    Encryptor encryptor(context); encryptor.set_public_key(public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);


    vector<complex<double>> input;
    input.reserve(slot_count);
    double curr_point = -2;
    //double step_size = 1.0 / (static_cast<double>(slot_count) - 1);
    double step_size = 0.25;
    for (size_t i = 0; i < 16; i++)
    {
        input.push_back(curr_point);
        curr_point += step_size;
    }
    cout << "Input vector: " << endl;
    print_vector(input, 16, 7);

    cout << "Evaluating polynomial x^2 + 7x + 10 ..." << endl;

    /*
        Try to calculate th
    */


    Plaintext x_plain;
    print_line(__LINE__);
    cout << "Encode input vectors." << endl;
    encoder.encode_complex64_simd(input, std::nullopt, scale, x_plain);
    Ciphertext x_encrypted;
    encryptor.encrypt_asymmetric(x_plain, x_encrypted);
    Ciphertext encrypted_result;
    approx_with_fix(encoder, evaluator, relin_keys, scale, x_encrypted, encrypted_result);

    /*
    Decrypt, decode, and print the result.
    */
    Plaintext plain_result;
    decryptor.decrypt(encrypted_result, plain_result);
    print_line(__LINE__);
    vector<complex<double>> result;
    encoder.decode_complex64_simd(plain_result, result);
    cout << "    + Computed result ...... Correct." << endl;
    cout << "    size of the result: " << result.size() << endl;
    //print_vector(result, 16, 7);
    for (int i = 0; i < 16; i++) {
        cout << "Relu(x) / calc: " << relu(input[i].real()) << " " << result[i].real() / 15 << endl;
    }
    cout << endl;
    return 0;

}
