#pragma once
#include "ReLU_ckks.h"
#include "../examples/examples.h"
#include "../Polynomial_Calc/polynomial.h"
#include "../Polynomial_Calc/SiLU.h"
#include <iostream>
#include <chrono>
#include <ctime>
#include <codecvt>
#include <locale>
//#include "./matplotlib-cpp/matplotlibcpp.h"
//#define MYDEBUG
#define BATCH_SIZE 40
#define PARTS 2
using namespace troy;
using namespace std;
//namespace plt = matplotlibcpp;
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
            const Polynomial<double> &poly, const Ciphertext &x_encrypted, Ciphertext &encrypted_result) {

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


void multiply_x(const CKKSEncoder &encoder, const Evaluator &evaluator, const RelinKeys &relin_keys, 
    const double &scale, const int &degree, const Ciphertext &x_encrypted, Ciphertext &encrypted_result) {
    Ciphertext my_cipher = x_encrypted; // set the initial value to be x

    int deg_left = degree - 1; // Its x ^ 1 at the beginning
    while (deg_left > 0) {
        if (deg_left % 2 == 1) {
            // multiply x
            Ciphertext this_x = x_encrypted;
            // First change scale 
            evaluator.mod_switch_to_inplace(this_x, my_cipher.parms_id());
            this_x.scale() = my_cipher.scale();
            // Second multiply
            evaluator.multiply_inplace(my_cipher, this_x);
            evaluator.relinearize_inplace(my_cipher, relin_keys);
            evaluator.rescale_to_next_inplace(my_cipher);
            deg_left -= 1;
        }
        else {
            Ciphertext tmp;
            evaluator.square(my_cipher, tmp);
            evaluator.relinearize_inplace(tmp, relin_keys);
            evaluator.rescale_to_next_inplace(tmp);
            my_cipher = tmp;
            deg_left >>= 1;  
        }    
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

    /*
        Create multiple threads
    */
    vector<Polynomial<double>> fragments;
    for (int i = 0; i < PARTS; i++) {
        Polynomial<double> tmp = poly.slice((poly.get_degree() + 1) * i / PARTS, ((poly.get_degree() + 1) * (i + 1) / PARTS) - 1);
        fragments.push_back(tmp);
        tmp.check();
    }

    auto start = std::chrono::high_resolution_clock::now();
    vector<Ciphertext> ciphers(2 * PARTS - 1);

    vector<thread> threads;
    mutex mtx;

    for (int i = 0; i < (2 * PARTS - 1); i++) {
        if (i % 2 == 0) {
            threads.push_back(thread(horner, cref(encoder), cref(evaluator), cref(relin_keys), 
                                    cref(scale), cref(fragments[i / 2]), cref(x_encrypted), ref(ciphers[i])));
        } else {
            int degree = poly.get_degree() + 1;
            int part_degree = degree * (i + 1) / (PARTS * 2);
            threads.push_back(thread(multiply_x, cref(encoder), cref(evaluator), cref(relin_keys), 
                                    cref(scale), cref(part_degree), cref(x_encrypted), ref(ciphers[i])));
        }
    }

    for (auto &t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }

#ifdef MYDEBUG
    cout << "Debugging..." << endl;
    for (int i = 0; i < (2 * PARTS - 1); i++) {
        cout << "i = " << i << endl;
        Plaintext plain_result;
        decryptor.decrypt(ciphers[i], plain_result);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        vector<complex<double>> result;
        encoder.decode_complex64_simd(plain_result, result);
        for (int j = 0; j < BATCH_SIZE; j++) {
            //cout << "x: " << input[i].real() << "       Relu(x): " << relu(input[i].real()) << "      SiLU(x): " << silu(input[i].real()) << "       计算结果：" << result[i].real() << endl;
            cout << "x: " << input[j].real() << "       Relu(x): " << relu(input[j].real()) << "       计算结果：" << result[j].real();
            if (i % 2 == 0) {
                cout << "     True value: " << fragments[i / 2].get_poly_value(input[j].real());
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << "Debugging end..." << endl;
#endif

    vector<thread> mult(PARTS - 1);

    for (int i = 1; i < (2 * PARTS - 1); i += 2) {
        // multiply x
        mult.push_back(thread([&evaluator, &relin_keys, &ciphers, i]() mutable {
            Ciphertext this_x = ciphers[i];
            // First change scale 
            evaluator.mod_switch_to_inplace(this_x, ciphers[i + 1].parms_id());
            this_x.scale() = ciphers[i + 1].scale();
            // Second multiply
            evaluator.multiply_inplace(ciphers[i + 1], this_x);
            evaluator.relinearize_inplace(ciphers[i + 1], relin_keys);
            evaluator.rescale_to_next_inplace(ciphers[i + 1]);
        }));
    }

    for (auto &t : mult) {
        if (t.joinable()) {
            t.join();
        }
    }

    Ciphertext enc;

    

    for (int i = 2; i < 2 * PARTS - 1; i += 2) {
        Ciphertext this_x = ciphers[i - 2];
        evaluator.mod_switch_to_inplace(this_x, ciphers[i].parms_id());
        this_x.scale() = ciphers[i].scale();
        evaluator.add_inplace(ciphers[i], this_x);
        if ((i + 2) > (2 * PARTS - 1)) {
            enc = ciphers[i];
        }
    }

    //approx_with_fix(encoder, evaluator, relin_keys, scale, x_encrypted, encrypted_result);
    
    // horner(encoder, evaluator, relin_keys, scale, poly, x_encrypted, encrypted_result);
    
    /*
    Decrypt, decode, and print the result.
    */
    Plaintext plain_result;
    decryptor.decrypt(enc, plain_result);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    vector<complex<double>> result;
    encoder.decode_complex64_simd(plain_result, result);
    //print_vector(result, 16, 7);

    /*
        Ploting the curve: 
    */

    vector<double> x_coord, relu_y, poly_y;
    for (int i = 0; i < BATCH_SIZE; i++) {
        //cout << "x: " << input[i].real() << "       Relu(x): " << relu(input[i].real()) << "      SiLU(x): " << silu(input[i].real()) << "       计算结果：" << result[i].real() << endl;
        cout << "x: " << input[i].real() << "       Relu(x): " << relu(input[i].real()) << "       计算结果：" << result[i].real() << endl;
        x_coord.push_back(input[i].real());
        relu_y.push_back(relu(input[i].real()));
        poly_y.push_back(result[i].real());
    }
    std::cout << "运行时间: " << (double)duration.count() / (double)1000 << " ms" << std::endl;
    // Py_Initialize();
    
    // // 获取 Python 解释器路径
    // const wchar_t* python_executable_w = Py_GetProgramName();
    
    // // 将 wchar_t* 转换为 std::string
    // std::wstring_convert<std::codecvt_utf8<wchar_t>> converter;
    // std::string python_executable = converter.to_bytes(python_executable_w);

    // std::cout << "Python executable path: " << python_executable << std::endl;
    
    // Py_Finalize();

    // plt::plot(x_coord, relu_y);
    // plt::plot(x_coord, poly_y);  
    // plt::show();
    return 0;

}
