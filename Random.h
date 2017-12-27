#ifndef LIBE_RANDOM_H
#define LIBE_RANDOM_H

#include "Sampling.h"
#include "params.h"

vec_ZZ RandomVector();
ZZX RandomPoly(const unsigned int degree);
ZZX RandomPolyFixedSqNorm(const ZZ& SqNorm, const unsigned int degree);
std::string Gen_Random_Str(const int);
#endif
