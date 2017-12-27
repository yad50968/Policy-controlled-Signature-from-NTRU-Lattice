/*

Copyright or © or Copr. Liu Zi Yuan.

yad50986@gmail.com

This software is a computer program which purpose is to provide to the
research community a proof-of-concept implementation of the Policy-controlled
signature from NTRU lattices, described in the NTMS2018 paper
"Policy-controlled signature from NTRU lattices", of
Zi-Yuan Liu, Jen-Chieh Hsu, Raylin Tso

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

*/

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
#include "Scheme.h"
#include "io.h"
#include "params.h"

using namespace std;
using namespace NTL;

//==============================================================================
//==============================================================================
//                                  MAIN
//==============================================================================
//==============================================================================

int main() {

    ZZX MSK[4], SSK[4];
    ZZ_pX phiq, MPK, SPK;
    unsigned int i;
    MSK_Data* MSKD = new MSK_Data;
    MPK_Data* MPKD = new MPK_Data;

    MSK_Data* SSKD = new MSK_Data;
    MPK_Data* SPKD = new MPK_Data;

    const ZZX phi = Cyclo();

    srand(rdtsc()); // initialisation of rand

    ZZ_p::init(q1);
    zz_p::init(q0);

    phiq = conv<ZZ_pX>(phi);
    ZZ_pXModulus PHI(phiq);

    // Generate TA/Signer's key
    for (i = 0; i < 1; i++) {
        Keygen(MPK, MSK);
        Keygen(SPK, SSK);
    }

    CompleteMSK(MSKD, MSK);
    CompleteMPK(MPKD, MPK);

    CompleteMSK(SSKD, SSK);
    CompleteMPK(SPKD, SPK);

    // Sign
    int policy_size = 2;
    std::vector<std::string> policy;
    policy.push_back("Alice and Bob");
    policy.push_back("Tom");
    ZZX policy_SK[2];

    std::string message = "message";
    long int C_i[policy_size][2][N0];
    std::vector<std::string> R_i;
    std::string message_p_hash;
    ZZX signature[2];

    Signature(policy, message, MPKD, SSKD, signature, message_p_hash, C_i, R_i);

    // Verify apply the policy key
    int policy_key_pos = 1; // which policy he want to apply
    CreGen(policy[policy_key_pos], policy_SK, MSKD);

    // Verify
    Verify(policy, policy_SK, message, SPKD, signature, message_p_hash, C_i,
           R_i, policy_key_pos);

    return 0;
}
