#pragma once
/**
 * Public API for QAOA symmetric tree contraction (forward evaluation only).
 */

#include <complex>
#include <vector>

using C128 = std::complex<double>;

/**
 * Exact <Z^{otimes k}> for depth-p QAOA on a D-regular k-uniform tree.
 *
 * @param gammas  Phase separator angles, length p.
 * @param betas   Mixer angles, length p.
 * @param p       QAOA depth.
 * @param D       Vertex degree.
 * @param k       Hyperedge size.
 * @param use_dd  If true, use double-double precision.
 * @param verbose If true, print timing to stderr.
 * @return The expectation value <Z^{otimes k}>.
 */
double contract_symmetric_tree(const double* gammas, const double* betas,
                               int p, int D, int k,
                               bool use_dd = false, bool verbose = false);
