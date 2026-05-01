#ifndef MATRIX_ENSEMBLES_H
#define MATRIX_ENSEMBLES_H

#include "bitmatrix.hpp"

//parity check matrix from Gallager's ensemble
bitmatrix gallager(size_t k, size_t D, size_t bsize, std::mt19937_64 &eng);

//F_2 matrix from erdos_renyi ensemble. p=0.5 is uniformly random matrix.
bitmatrix erdos_renyi(size_t n, size_t m, double p, std::mt19937_64 &eng);

#endif