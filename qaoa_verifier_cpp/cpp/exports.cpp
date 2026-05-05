/**
 * C-linkage exports for Python ctypes binding (forward evaluation only).
 */

#include "contract.h"

extern "C" {

double qaoa_contract(const double* gammas, const double* betas,
                     int p, int D, int k,
                     int use_dd, int verbose) {
    return contract_symmetric_tree(gammas, betas, p, D, k,
                                   use_dd != 0, verbose != 0);
}

}  // extern "C"
