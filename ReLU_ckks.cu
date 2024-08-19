#pragma once
#include "ReLU_ckks.h"
#include "../examples/examples.h"
#include "../Polynomial_Calc/polynomial.h"
#include "../Polynomial_Calc/SiLU.h"
#include <iostream>

using namespace troy;
using namespace std;

int main() {

    EncryptionParameters parms(SchemeType::CKKS);

    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::create(poly_modulus_degree, { 60, 40, 40, 60 }));

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
    double curr_point = 0;
    //double step_size = 1.0 / (static_cast<double>(slot_count) - 1);
    double step_size = 0.25;
    for (size_t i = 0; i < 10; i++)
    {
        input.push_back(curr_point);
        curr_point += step_size;
    }
    cout << "Input vector: " << endl;
    print_vector(input, 10, 7);

    cout << "Evaluating polynomial PI*x^3 + 0.4x + 1 ..." << endl;

    /*
    We create plaintexts for PI, 0.4, and 1 using an overload of CKKSEncoder::encode
    that encodes the given floating-point value to every slot in the vector.
    */
    Plaintext plain_coeff2, plain_coeff1, plain_coeff0;
    encoder.encode_float64_single(3.0, std::nullopt, scale, plain_coeff2);
    encoder.encode_float64_single(2.0, std::nullopt, scale, plain_coeff1);
    encoder.encode_float64_single(1.0, std::nullopt, scale, plain_coeff0);

    Plaintext x_plain;
    print_line(__LINE__);
    cout << "Encode input vectors." << endl;
    encoder.encode_complex64_simd(input, std::nullopt, scale, x_plain);
    Ciphertext x1_encrypted;
    encryptor.encrypt_asymmetric(x_plain, x1_encrypted);

    Ciphertext x3_encrypted;
    print_line(__LINE__);
    cout << "Compute x^2 and relinearize:" << endl;
    evaluator.square(x1_encrypted, x3_encrypted);
    evaluator.relinearize_inplace(x3_encrypted, relin_keys);
    cout << "    + Scale of x^2 before rescale: " << log2(x3_encrypted.scale()) << " bits" << endl;

    print_line(__LINE__);
    cout << "Rescale x^2." << endl;
    evaluator.rescale_to_next_inplace(x3_encrypted);
    cout << "    + Scale of x^2 after rescale: " << log2(x3_encrypted.scale()) << " bits" << endl;


    // Ensure the plain_coeff2 has the same parms_id as the rescaled x3_encrypted
    Plaintext plain_coeff2_adjusted;
    encoder.encode_float64_single(3.0, x3_encrypted.parms_id(), x3_encrypted.scale(), plain_coeff2_adjusted);

    // Then multiply the rescaled ciphertext with the adjusted plaintext
    Ciphertext x1_encrypted_coeff3;
    evaluator.multiply_plain(x3_encrypted, plain_coeff2_adjusted, x1_encrypted_coeff3);

    cout << "    + Scale of 3x^2 before rescale: " << log2(x1_encrypted_coeff3.scale()) << " bits" << endl;
    evaluator.rescale_to_next_inplace(x1_encrypted_coeff3);
    cout << "    + Scale of 3x^2 after rescale: " << log2(x1_encrypted_coeff3.scale()) << " bits" << endl;
    print_line(__LINE__);
    /*
    Decrypt, decode, and print the result.
    */
    Plaintext plain_result;
    decryptor.decrypt(x1_encrypted_coeff3, plain_result);
    print_line(__LINE__);
    vector<complex<double>> result;
    encoder.decode_complex64_simd(plain_result, result);
    cout << "    + Computed result ...... Correct." << endl;
    cout << "    size of the result: " << result.size() << endl;
    print_vector(result, 10, 7);

    return 0;

}
