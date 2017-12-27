#ifndef LIBE_SCHEME_H
#define LIBE_SCHEME_H

#include "Algebra.h"
#include "FFT.h"
#include "Random.h"
#include "Sampling.h"
#include "params.h"
#include <vector>
using namespace std;

long *string_to_long_array(long *la, std::string s);
std::string long_array_to_string(long *la);
std::string binary_xor(std::string s1, std::string s2);
const char *hex_char_to_bin(char c);
std::string hex_str_to_bin_str(const std::string &hex);

void Keygen(ZZ_pX &PublicKey, ZZX *PrivateKey);
void CompletePrivateKey(mat_ZZ &B, const ZZX *const PrivateKey);
void GPV(RR_t *v, const RR_t *const c, const RR_t s,
         const MSK_Data *const MSKD);
void CompleteMSK(MSK_Data *MSKD, ZZX *MSK);
void CompleteMPK(MPK_Data *MPKD, ZZ_pX MPK);
void IBE_Extract(ZZX SK_id[2], vec_ZZ id, const MSK_Data *const MSKD);
unsigned long IBE_Verify_Key(const ZZX SK_id[2], const vec_ZZ id,
                             const MSK_Data *const MSKD);
void IBE_Encrypt(long C[2][N0], const long m[N0], const long id0[N0],
                 const MPK_Data *const MPKD);
void IBE_Decrypt(long message[N0], const long C[2][N0],
                 const CC_t *const SKid_FFT);
void CreGen(std::string policy, ZZX SK_policy[2], const MSK_Data *const MSKD);
void Signature(std::vector<std::string> policy, std::string message,
               const MPK_Data *const MPKD, const MSK_Data *const SSKD, ZZX *,
               std::string &, long int C_i[][2][N0],
               std::vector<std::string> &R_i);
void Verify(std::vector<std::string> policy_size, ZZX policy_SK[2],
            std::string message, const MPK_Data *const SPKD, ZZX *signature,
            std::string &message_p_hash, long int C_i[][2][N0],
            std::vector<std::string> R_i, int policy_key_pos);
#endif
