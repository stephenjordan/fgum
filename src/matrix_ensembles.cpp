#include "matrix_ensembles.hpp"
#include "bitmatrix.hpp"

//parity check matrix from Gallager's ensemble
bitmatrix gallager(size_t k, size_t D, size_t bsize, std::mt19937_64 &eng) {
    size_t m = D*bsize;
    size_t n = k*bsize;
    bitmatrix B(n,m);
    std::vector<int> perm(m);
    for(size_t i = 0; i < m; i++) perm[i] = i;
    for(size_t a = 0; a < k; a++) {
        std::shuffle(perm.begin(), perm.end(), eng);
        for(size_t b = 0; b < bsize; b++) {
            size_t row_index = a*bsize+b;
            for(size_t entry_index = 0; entry_index < D; entry_index++) {
                B.set(row_index, perm[bsize*entry_index + row_index%bsize], 1);
            }
        }
    }
    return B;
}

//F_2 matrix from erdos_renyi ensemble. p=0.5 is uniformly random matrix.
bitmatrix erdos_renyi(size_t n, size_t m, double p, std::mt19937_64 &eng) {
    bitmatrix H(n,m);
    std::bernoulli_distribution coin(p);
    for(size_t i = 0; i < n; i++) {
        for(size_t j = 0; j < m; j++) {
            H.set(i,j,coin(eng));
        }
    }
    return H;
}