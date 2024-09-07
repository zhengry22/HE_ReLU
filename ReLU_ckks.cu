#include "ReLU_ckks.h"
#include "../examples/examples.h"
#include "../Polynomial_Calc/polynomial.h"
#include "../Polynomial_Calc/SiLU.h"
#include <iostream>
#include <chrono>
#include <ctime>
#define BATCH_SIZE 40
using namespace troy;
using namespace std;

/*
    Try to calculate the cipher with a fixed polynomial 
*/ 
void approx_with_fix(const CKKSEncoder &encoder, const Evaluator &evaluator, const RelinKeys &relin_keys, const double &scale,
                    const Ciphertext &x_encrypted, Ciphertext &encrypted_result) {
    
    // Prepare the plaintext for 2, 5, and 1
    Plaintext plain_2, plain_5, plain_1;
    encoder.encode_float64_single(2.0, std::nullopt, scale, plain_2);
    encoder.encode_float64_single(5.0, std::nullopt, scale, plain_5);
    encoder.encode_float64_single(1.0, std::nullopt, scale, plain_1);

    // Compute x + 2
    Ciphertext x_plus_2 = x_encrypted;
    evaluator.add_plain_inplace(x_plus_2, plain_2);

    // Compute x + 5
    Ciphertext x_plus_5 = x_encrypted;
    evaluator.add_plain_inplace(x_plus_5, plain_5);

    // Compute x + 1
    Ciphertext x_plus_1 = x_encrypted;
    evaluator.add_plain_inplace(x_plus_1, plain_1);

    // Compute (x + 2) * (x + 5)
    evaluator.multiply_inplace(x_plus_2, x_plus_5);
    evaluator.relinearize_inplace(x_plus_2, relin_keys);
    evaluator.rescale_to_next_inplace(x_plus_2);

    // Compute (x + 2) * (x + 5) * (x + 1)
    evaluator.mod_switch_to_inplace(x_plus_1, x_plus_2.parms_id());
    x_plus_1.scale() = x_plus_2.scale();
    evaluator.multiply_inplace(x_plus_2, x_plus_1);
    evaluator.relinearize_inplace(x_plus_2, relin_keys);
    evaluator.rescale_to_next_inplace(x_plus_2);

    // Compute (x + 2) * (x + 5) * (x + 1) * (x + 1)
    evaluator.mod_switch_to_inplace(x_plus_1, x_plus_2.parms_id());
    x_plus_1.scale() = x_plus_2.scale();
    evaluator.multiply_inplace(x_plus_2, x_plus_1);
    evaluator.relinearize_inplace(x_plus_2, relin_keys);
    evaluator.rescale_to_next_inplace(x_plus_2);

    evaluator.mod_switch_to_inplace(x_plus_1, x_plus_2.parms_id());
    x_plus_1.scale() = x_plus_2.scale();
    evaluator.multiply_inplace(x_plus_2, x_plus_1);
    evaluator.relinearize_inplace(x_plus_2, relin_keys);
    evaluator.rescale_to_next_inplace(x_plus_2);

    evaluator.mod_switch_to_inplace(x_plus_1, x_plus_2.parms_id());
    x_plus_1.scale() = x_plus_2.scale();
    evaluator.multiply_inplace(x_plus_2, x_plus_1);
    evaluator.relinearize_inplace(x_plus_2, relin_keys);
    evaluator.rescale_to_next_inplace(x_plus_2);

    // Store the result
    encrypted_result = x_plus_2;
}


void horner(const CKKSEncoder &encoder, const Evaluator &evaluator, const RelinKeys &relin_keys, const double &scale,
            const Polynomial<double> poly, const Ciphertext &x_encrypted, Ciphertext &encrypted_result) {

    /*
        Conduct the horner's algorithm in linear complexity
    */
    Ciphertext my_cipher = x_encrypted; // set the initial value to be x
    size_t poly_deg = poly.get_degree();
    Plaintext plain_coeff;
    double largest_coeff = poly.get_coeff_by_rank(poly_deg);
    //double largest_coeff = 1.0;
    encoder.encode_float64_single(largest_coeff, std::nullopt, scale, plain_coeff);
    evaluator.multiply_plain_inplace(my_cipher, plain_coeff);
    evaluator.rescale_to_next_inplace(my_cipher); // a_n*x

    for (int i = poly_deg - 1; i >= 0; i--) {
        if (i != poly_deg - 1) {
            // multiply x
            Ciphertext this_x = x_encrypted;
            // First change scale 
            evaluator.mod_switch_to_inplace(this_x, my_cipher.parms_id());
            this_x.scale() = my_cipher.scale();
            // Second multiply
            evaluator.multiply_inplace(my_cipher, this_x);
            evaluator.relinearize_inplace(my_cipher, relin_keys);
            evaluator.rescale_to_next_inplace(my_cipher);
        }
        // Add const
        Plaintext this_coeff;
        encoder.encode_float64_single(poly.get_coeff_by_rank((size_t)(i)), std::nullopt, my_cipher.scale(), this_coeff);
        //encoder.encode_float64_single(0.0, std::nullopt, my_cipher.scale(), this_coeff);
        evaluator.add_plain_inplace(my_cipher, this_coeff);
    }

    encrypted_result = my_cipher;
}


int main() {

    int deg;
    cout << "Input deg: " << endl;
    cin >> deg;

    vector<size_t> mod_chain;
    for (int i = 0; i < deg + 2; i++) {
        if (i == 0) {
            mod_chain.push_back(60);
        }
        else if (i == deg + 1) {
            mod_chain.push_back(60);
        }
        else {
            mod_chain.push_back(40);
        }
    }

    EncryptionParameters parms(SchemeType::CKKS);

    Remez<double, double> my_p(deg, gelu_and_sqplus);
    Polynomial<double> poly = my_p.generate_approx(deg, 0);
    poly.prune();
    poly.check();

    size_t poly_modulus_degree = 32768;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::create(poly_modulus_degree, mod_chain));
    // Change this modulus chain if we need higher degree


    double scale = pow(2.0, 40);

    auto context = HeContext::create(parms, true, SecurityLevel::Classical128);
    print_parameters(*context);
    cout << endl;

    CKKSEncoder encoder(context);
    size_t slot_count = encoder.slot_count();
    //cout << "Number of slots: " << slot_count << endl;

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
    double curr_point = -5;
    //double step_size = 1.0 / (static_cast<double>(slot_count) - 1);
    double step_size = (double)10 / (double)BATCH_SIZE;
    for (size_t i = 0; i < BATCH_SIZE; i++)
    {
        input.push_back(curr_point);
        curr_point += step_size;
    }
    cout << "Input vector: " << endl;
    print_vector(input, BATCH_SIZE, 7);

    /*
        Try to calculate th
    */


    Plaintext x_plain;
    cout << "Encode input vectors." << endl;
    encoder.encode_complex64_simd(input, std::nullopt, scale, x_plain);
    Ciphertext x_encrypted;
    encryptor.encrypt_asymmetric(x_plain, x_encrypted);
    Ciphertext encrypted_result;
    //approx_with_fix(encoder, evaluator, relin_keys, scale, x_encrypted, encrypted_result);
    auto start = std::chrono::high_resolution_clock::now();
    horner(encoder, evaluator, relin_keys, scale, poly, x_encrypted, encrypted_result);
    
    /*
    Decrypt, decode, and print the result.
    */
    Plaintext plain_result;
    decryptor.decrypt(encrypted_result, plain_result);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    vector<complex<double>> result;
    encoder.decode_complex64_simd(plain_result, result);
    //print_vector(result, 16, 7);
    for (int i = 0; i < BATCH_SIZE; i++) {
        //cout << "x: " << input[i].real() << "       Relu(x): " << relu(input[i].real()) << "      SiLU(x): " << silu(input[i].real()) << "       计算结果：" << result[i].real() << endl;
        //cout << "x: " << input[i].real() << "       Relu(x): " << relu(input[i].real()) << "    abs(x): " << abs_test(input[i].real()) << "       计算结果：" << result[i].real() << endl;
        cout << "x: " << input[i].real() << "       Relu(x): " << relu(input[i].real()) << "       计算结果：" << result[i].real() << endl;
    }
    std::cout << "运行时间: " << (double)duration.count() / (double)1000 << " ms" << std::endl;
    cout << endl;
    return 0;

}
