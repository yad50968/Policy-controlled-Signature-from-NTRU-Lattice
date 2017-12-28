#include "sha512.h"
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <assert.h>
#include <complex.h>
#include <gmp.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <vector>

#include "Algebra.h"
#include "FFT.h"
#include "Random.h"
#include "Sampling.h"
#include "params.h"

using namespace std;
using namespace NTL;

const ZZX phi = Cyclo();

long* string_to_long_array(long* la, std::string s) {

    for (int i = 0; i < N0; i++) {
        la[i] = s[i] - '0';
    }
    return la;
}

std::string long_array_to_string(long* la) {

    std::string s;
    for (int i = 0; i < N0; i++) {
        s += la[i] + '0';
    }
    return s;
}

std::string binary_xor(std::string s1, std::string s2) {
    std::string output;
    for (int i = 0; i < s2.length(); i++)
        output += ((s1[i] - '0') ^ (s2[i] - '0')) + '0';
    return output;
}

const char* hex_char_to_bin(char c) {
    // TODO handle default / error
    switch (toupper(c)) {
    case '0':
        return "0000";
    case '1':
        return "0001";
    case '2':
        return "0010";
    case '3':
        return "0011";
    case '4':
        return "0100";
    case '5':
        return "0101";
    case '6':
        return "0110";
    case '7':
        return "0111";
    case '8':
        return "1000";
    case '9':
        return "1001";
    case 'A':
        return "1010";
    case 'B':
        return "1011";
    case 'C':
        return "1100";
    case 'D':
        return "1101";
    case 'E':
        return "1110";
    case 'F':
        return "1111";
    }
}

std::string hex_str_to_bin_str(const std::string& hex) {
    // TODO use a loop from <algorithm> or smth
    std::string bin;
    for (unsigned i = 0; i != hex.length(); ++i)
        bin += hex_char_to_bin(hex[i]);
    return bin;
}
//==============================================================================
// Generates from parameters N and q :
// - a public key : polynomial h
// - a private key : polynomials f,g,F,G
//==============================================================================
void Keygen(ZZ_pX& PublicKey, ZZX* PrivateKey) {
    ZZ SqNorm;
    ZZX f, g, F, G;

    SqNorm = conv<ZZ>(1.36 * q0 / 2);

    GenerateBasis(f, g, F, G, SqNorm);
    PrivateKey[0] = f;
    PrivateKey[1] = g;
    PrivateKey[2] = F;
    PrivateKey[3] = G;

    for (unsigned int i = 0; i < 4; i++) {
        PrivateKey[i].SetLength(N0);
    }

    PublicKey = Quotient(f, g);
}

//==============================================================================
// Computes the private basis B from private key PrivateKey and parameter N
//==============================================================================
void CompletePrivateKey(mat_ZZ& B, const ZZX* const PrivateKey) {
    ZZX f, g, F, G;
    f = PrivateKey[0];
    g = PrivateKey[1];
    F = PrivateKey[2];
    G = PrivateKey[3];

    f = -f;
    F = -F;

    B = BasisFromPolynomials(g, f, G, F);
}

void GPV(RR_t* v, const RR_t* const c, const RR_t s,
         unique_ptr<MSK_Data> &MSKD) {

    int i;
    unsigned j;
    RR_t ci[2 * N0], zi, cip, sip, aux;

    for (j = 0; j < 2 * N0; j++) {
        ci[j] = c[j];
    }

    for (j = 0; j < 2 * N0; j++) {
    }

    for (i = 2 * N0 - 1; i >= 0; i--) {
        aux = (MSKD->GS_Norms)[i];
        cip = DotProduct(ci, MSKD->Bstar[i]) / (aux * aux);
        sip = s / aux;
        zi = Sample4(cip, sip * PiPrime);

        for (j = 0; j < 2 * N0; j++) {
            ci[j] -= zi * (MSKD->B)[i][j];
        }
    }

    for (j = 0; j < 2 * N0; j++) {
        v[j] = c[j] - ci[j];
    }
}

//==============================================================================
//==============================================================================
//                            MAIN PROGRAMS
//==============================================================================
//==============================================================================

void CompleteMSK(unique_ptr<MSK_Data> &MSKD, ZZX* MSK) {
    unsigned int i, j;
    mat_ZZ B0;

    for (i = 0; i < 4; i++) {
        MSKD->PrK[i] = MSK[i];
        ZZXToFFT(MSKD->PrK_fft[i], MSK[i]);
    }

    CompletePrivateKey(B0, MSK);

    for (i = 0; i < 2 * N0; i++) {
        for (j = 0; j < 2 * N0; j++) {
            MSKD->B[i][j] = ((RR_t)conv<double>(B0[i][j]));
        }
    }

    for (i = 0; i < 1; i++) {
        FastMGS(MSKD->Bstar, MSKD->B);
    }

    for (i = 0; i < 2 * N0; i++) {
        MSKD->GS_Norms[i] = sqrt(DotProduct(MSKD->Bstar[i], MSKD->Bstar[i]));
    }

    MSKD->sigma = 2 * MSKD->GS_Norms[0];
}

void CompleteMPK(unique_ptr<MPK_Data> &MPKD, ZZ_pX MPK) {
    MPKD->h = MPK;
    ZZXToFFT(MPKD->h_FFT, conv<ZZX>(MPK));
}

void IBE_Extract(ZZX SK_id[2], vec_ZZ id, unique_ptr<MSK_Data> &MSKD) {
    unsigned int i;
    RR_t c[2 * N0], sk[2 * N0], sigma;
    ZZX f, g, aux;

    f = MSKD->PrK[0];
    g = MSKD->PrK[1];
    sigma = MSKD->sigma;
    SK_id[0].SetLength(N0);
    SK_id[1].SetLength(N0);

    for (i = 0; i < N0; i++) {
        c[i] = ((RR_t)conv<double>(id[i]));
        c[i + N0] = 0;
    }

    GPV(sk, c, sigma, MSKD);

    for (i = 0; i < N0; i++) {
        sk[i] = c[i] - sk[i];
        sk[i + N0] = -sk[i + N0];
    }

    for (i = 0; i < N0; i++) {
        SK_id[0][i] = sk[i];
        SK_id[1][i] = sk[i + N0];
    }
}


void IBE_Encrypt(long C[2][N0], const long m[N0], const long id0[N0],
                 unique_ptr<MPK_Data> &MPKD) {

    unsigned long i;
    long r[N0], e1[N0], e2[N0];
    CC_t r_FFT[N0], t_FFT[N0], aux1_FFT[N0], aux2_FFT[N0];

    for (i = 0; i < N0; i++) {
        e1[i] = (rand() % 3) - 1;
        e2[i] = (rand() % 3) - 1;
        r[i] = (rand() % 3) - 1;
    }

    MyIntFFT(r_FFT, r);
    MyIntFFT(t_FFT, id0);

    for (i = 0; i < N0; i++) {
        aux1_FFT[i] = r_FFT[i] * ((MPKD->h_FFT)[i]);
        aux2_FFT[i] = r_FFT[i] * t_FFT[i];
    }

    MyIntReverseFFT(C[0], aux1_FFT);
    MyIntReverseFFT(C[1], aux2_FFT);

    for (i = 0; i < N0; i++) {
        C[0][i] = (C[0][i] + e1[i] + q0 / 2) % q0 - (q0 / 2);
        C[1][i] = (C[1][i] + e2[i] + (q0 / 2) * m[i] + q0 / 2) % q0 - (q0 / 2);
    }
}

void IBE_Decrypt(long message[N0], const long C[2][N0],
                 const CC_t* const SKid_FFT) {
    unsigned int i;
    CC_t c0_FFT[N0], aux_FFT[N0];

    MyIntFFT(c0_FFT, C[0]);

    for (i = 0; i < N0; i++) {
        aux_FFT[i] = c0_FFT[i] * SKid_FFT[i];
    }

    MyIntReverseFFT(message, aux_FFT);

    for (i = 0; i < N0; i++) {
        message[i] = C[1][i] - message[i];
        message[i] = ((unsigned long)(message[i])) % q0;
        message[i] = (message[i] + (q0 >> 2)) / (q0 >> 1);
        message[i] %= 2;
    }
}

void CreGen(std::string policy, ZZX SK_policy[2], unique_ptr<MSK_Data> &MSKD) {
    vec_ZZ policy_id = StringToVecZZ(policy);
    IBE_Extract(SK_policy, policy_id, MSKD);
}

void Signature(std::vector<std::string> policy, std::string message,
               unique_ptr<MPK_Data> &MPKD, unique_ptr<MSK_Data> &SSKD,
               ZZX* signature, std::string& message_p_hash,
               long int C_i[][2][N0], std::vector<std::string>& R_i) {

    vec_ZZ t = RandomVector();
    std::string m_i = Gen_Random_Str(N0);
    std::string message_p = message + VecZZToString(t);
    string message_p_p = message_p;

    int policy_size = policy.size();
    for (int i = 0; i < policy_size; i++) {
        vec_ZZ policy_id = StringToVecZZ(policy[i]);
        long int identity[N0];
        for (int i = 0; i < N0; i++) {
            identity[i] = conv<long int>(policy_id[i]);
        }

        long long_array_m_i[512];
        string_to_long_array(long_array_m_i, m_i);
        IBE_Encrypt(C_i[i], long_array_m_i, identity, MPKD);
        message_p_hash = sha512(message_p);
        std::string tmp_R_i =
            binary_xor(VecZZToString(t),
                       hex_str_to_bin_str(sha512(message_p_hash + "1" + m_i)));
        message_p_p += tmp_R_i;
        R_i.push_back(tmp_R_i);
    }

    vec_ZZ message_p_p_zz = StringToVecZZ(sha512(message_p_p));
    IBE_Extract(signature, message_p_p_zz, SSKD);
}

void Verify(std::vector<std::string> policy, ZZX policy_SK[2],
            std::string message, unique_ptr<MPK_Data> &SPKD, ZZX* signature,
            std::string& message_p_hash, long int C_i[][2][N0],
            std::vector<std::string> R_i, int policy_key_pos) {

    int policy_size = policy.size();
    CC_t policy_SK_FFT[N0];
    ZZXToFFT(policy_SK_FFT, policy_SK[1]);

    long int m_i_hat[N0];
    IBE_Decrypt(m_i_hat, C_i[policy_key_pos], policy_SK_FFT);

    std::string t_hat =
        binary_xor(R_i[policy_key_pos],
                   hex_str_to_bin_str(sha512(message_p_hash + "1" +
                                             long_array_to_string(m_i_hat))));

    std::string M_hat = message + t_hat;
    for (int i = 0; i < policy_size; i++) {
        M_hat += R_i[i];
    }

    ZZX tt = conv<ZZX>(StringToVecZZ(sha512(M_hat)));

    ZZX o;
    o.SetLength(N0);
    o = (signature[0] + signature[1] * (conv<ZZX>(SPKD->h))) % phi;

    for (int i = 0; i < N0; i++) {
        o[i] %= q1;
    }

    if (IsZero(o - tt) == 0) {
        cout << "Success" << endl;
    } else {
        cout << "Fail" << endl;
    }
}
