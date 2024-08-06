#pragma once
#include "ReLU_bfv.h"
#include <vector>

using namespace troy;
using namespace std;

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

int main() {

    // These are the parameters used for encryption: ReLU(x_) and ReLU(y_)
    uint64_t x_ = 0;

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

    cout << "We are here" << endl;
    poly_relu(encoder, evaluator, x_encrypted, encrypted_result);
    // /*
    //     Decrypt `encrypted result` and divide the result by 15
    // */

    Plaintext decrypted_result;
    decryptor.decrypt(encrypted_result, decrypted_result);
    std::string my_answer = decrypted_result.to_string();
    int int_value = std::stoi(my_answer, nullptr, 16);

    // 整数除以 15
    int result = int_value / 15;

    // 输出结果
    std::cout << "Relu(x) = " << result << endl;
    

    return 0;
}