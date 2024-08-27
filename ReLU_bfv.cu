#pragma once
#include "ReLU_bfv.h"
#include "../Polynomial_Calc/polynomial.h"
#include "../Polynomial_Calc/SiLU.h"
#include <vector>
#include <cmath>
using namespace troy;
using namespace std;
#define MYDEBUG

extern const size_t poly_modulus_degree;
// extern const vector<Modulus> coeff_modulus;
extern const size_t plain_modulus;


void poly_relu(const auto &encoder, const Evaluator &evaluator, const Ciphertext &x_encrypted, Ciphertext &encrypted_result) {
    /*
        This is the initial version of the poly_relu function.
        In the first version, we use x^2 + 7x + 10 to simulate the pt's polynomial * 15.
        In order to decrypt, we should first obtain the ciphertext and then divide the result by 15.
    */

    // First encode x + 2
    Ciphertext x_plus_two = x_encrypted; 
    Plaintext plain_two;
    encoder.encode_polynomial({2}, plain_two);
    evaluator.add_plain_inplace(x_plus_two, plain_two);

    // // Next encode x + 5
    Ciphertext x_plus_five = x_encrypted;
    Plaintext plain_five;
    encoder.encode_polynomial({5}, plain_five);
    evaluator.add_plain_inplace(x_plus_five, plain_five);

    // At last, multiply the two polynomials together

    evaluator.multiply(x_plus_two, x_plus_five, encrypted_result);
} 

void test_4degree_horner(const auto &encoder, const Evaluator &evaluator, const RelinKeys &relin_keys, const Ciphertext &x_encrypted, Ciphertext &encrypted_result) {
    
    /*
        This function is used to test whether using 4 degree polynomial has a problem
    */
    cout << "Testing horner's method using a 4 degree polynomial: " << endl;
    Plaintext zero;
    Ciphertext my_cipher = x_encrypted; // x
    encoder.encode_polynomial({0}, zero);

    evaluator.multiply_plain_inplace(my_cipher, zero); // x * 0
    evaluator.relinearize_inplace(my_cipher, relin_keys);
    Plaintext this_pt;
    encoder.encode_polynomial({1}, this_pt);
    evaluator.add_plain_inplace(my_cipher, this_pt); // x * 0 + 1

    Ciphertext mid;
    evaluator.multiply(my_cipher, x_encrypted, mid); // 1 * x
    evaluator.relinearize_inplace(mid, relin_keys); 
    my_cipher = mid;

    encoder.encode_polynomial({0}, this_pt); // 1 * x + 0
    evaluator.add_plain_inplace(my_cipher, this_pt);

    evaluator.multiply(my_cipher, x_encrypted, mid); // x^2
    evaluator.relinearize_inplace(mid, relin_keys); 
    my_cipher = mid;

    evaluator.multiply(my_cipher, x_encrypted, mid); // x^2
    evaluator.relinearize_inplace(mid, relin_keys); 
    my_cipher = mid;

    //     evaluator.multiply(my_cipher, x_encrypted, mid); // x^2
    // evaluator.relinearize_inplace(mid, relin_keys); 
    // my_cipher = mid;
    // encoder.encode_polynomial({8181}, this_pt); // x*2 - 11
    // evaluator.add_plain_inplace(my_cipher, this_pt);

    // evaluator.multiply(my_cipher, x_encrypted, mid); // x^3 - 11x
    // evaluator.relinearize_inplace(mid, relin_keys); 
    // my_cipher = mid;

    // encoder.encode_polynomial({8168}, this_pt); // x^3 - 11x - 24
    // evaluator.add_plain_inplace(my_cipher, this_pt);

    // evaluator.multiply(my_cipher, x_encrypted, mid); // x^4 - 11x^2 - 24x
    // evaluator.relinearize_inplace(mid, relin_keys); 
    // my_cipher = mid;

    // encoder.encode_polynomial({0}, this_pt); // x^4 - 11x^2 - 24x
    // evaluator.add_plain_inplace(my_cipher, this_pt);

    encrypted_result = my_cipher;
}

void horner(const auto &encoder, const Evaluator &evaluator, const RelinKeys &relin_keys, const EncryptPolynomial &encpoly, const Ciphertext &x_encrypted, Ciphertext &encrypted_result) {
    /*
        This function has bug in it! or the way we decrypt!
    */
    
    cout << "Using horners' method to calculate: " << endl;
    
    /*
        In order to calculate the final result, we may use the linear horner's algorithm. 
        Note that this is not likely the most efficient way, but is convenient for testing
    */

    cout << "Generate relinearization keys." << endl;

    int deg = encpoly.poly.get_degree();
    Plaintext zero;
    Ciphertext my_cipher = x_encrypted;
    encoder.encode_polynomial({0}, zero);
    for (int i = deg; i >= 0; i--) {
        // First multiple x and then add coeff[i]
        if (i == deg) {
            evaluator.multiply_plain_inplace(my_cipher, zero);
            evaluator.relinearize_inplace(my_cipher, relin_keys);
        }   
        else {
            Ciphertext mid;
            evaluator.multiply(my_cipher, x_encrypted, mid);
            evaluator.relinearize_inplace(mid, relin_keys);
            my_cipher = mid;
        }
        uint64_t this_coeff = encpoly.poly.get_coeff_by_rank(i);
        Plaintext this_pt;
        encoder.encode_polynomial({this_coeff}, this_pt);
        evaluator.add_plain_inplace(my_cipher, this_pt);
    }

    encrypted_result = my_cipher;
}


int main() {
    int deg;
    cout << "Input deg: " << endl;
    cin >> deg;

    // Generate the polynomial for approximation
    Taylor<double, double> taylor(deg, silu);
    Polynomial<double> poly = taylor.generate_approx(deg, 0);
    poly.prune();
    poly.check();
    EncryptPolynomial encpoly = round_polynomial(poly);
    encpoly.show();


    // These are the parameters used for encryption: ReLU(x_) and ReLU(y_)
    uint64_t x_ = 3;

    // Set the encryption scheme to be BFV
    EncryptionParameters parms(SchemeType::BFV);

    // Set the modulus of the polynomial modulus: x^d + 1, where d is a power of 2
    parms.set_poly_modulus_degree(poly_modulus_degree);

    // Set the coeff_modulus
    parms.set_coeff_modulus(CoeffModulus::bfv_default(poly_modulus_degree, SecurityLevel::Classical128));

    // Set the plain_modulus
    parms.set_plain_modulus(plain_modulus);

    auto context = HeContext::create(parms, true, SecurityLevel::Classical128);
    auto encoder = BatchEncoder(context);
    context->to_device_inplace();
    encoder.to_device_inplace();

    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key = keygen.create_public_key(false);

    Encryptor encryptor(context);
    encryptor.set_public_key(public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    /*
        First, encode x to poly_relu as poly_relu(x), 
        where poly_relu is the polynomial we use to simulate ReLU.
        This procedure should be done via these steps: 
            1. Encode the number x to 'x as a polynomial'
    */
    Plaintext x_plain;
    encoder.encode_polynomial({x_}, x_plain);
    /*
            2. Encrypt x
    */
    Ciphertext x_encrypted;
    encryptor.encrypt_asymmetric(x_plain, x_encrypted);
    /*
            3. Now we get the ciphertext c, calculate d, so that Relu(Dec(c)) = Dec(d)
    */
    Ciphertext encrypted_result;

    // poly_relu(encoder, evaluator, x_encrypted, encrypted_result);

    RelinKeys relin_keys = keygen.create_relin_keys(false);
    //horner(encoder, evaluator, relin_keys, encpoly, x_encrypted, encrypted_result);
    test_4degree_horner(encoder, evaluator, relin_keys, x_encrypted, encrypted_result);
    /*
        Decrypt `encrypted result` and divide the result by 15
    */

    Plaintext decrypted_result;
    decryptor.decrypt(encrypted_result, decrypted_result);
    std::string my_answer = decrypted_result.to_string();
    int int_value = std::stoi(my_answer, nullptr, 16);

    //int_value = ((plain_modulus >> 1) > int_value) ? int_value : int_value - plain_modulus;

    cout << int_value << " " << (double)int_value / (double)(encpoly.k) << endl;

    long long result = std::lround((double)int_value / (double)(encpoly.k));

    // 输出结果
    std::cout << "Relu(x) = " << result << endl;
    

    return 0;
}