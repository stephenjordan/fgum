#include <iostream>
#include <fstream>
#include <random>
#include "bitmatrix.hpp"

int save_max_xorsat_instance(const bitmatrix &H, const bitvector &v, std::string filename) {
    std::ofstream outfile;
    outfile.open(filename);
    if(!outfile.is_open()) {
        std::cout << "Error: unable to create " << std::endl;
        return 0;
    }
    outfile << H.num_rows() << "\t" << H.num_cols() << "\n";
    std::vector <std::vector <size_t> > scols = H.to_sparse_cols();
    for(size_t j = 0; j < H.num_cols(); j++) {
        double coeff = 0.5;
        if(v.get(j)) coeff = -0.5;
        outfile << coeff;
        for(size_t index : scols[j]) outfile << "\t" << index;
        outfile << "\n";
    }
    outfile.close();
    return 1;
}

int main(int argc, char *argv[]) {
    int k,D,bsize;
    if(argc != 4) {
        std::cout << "Usage: random_regular k D bsize" << std::endl;
        return 0;
    }
    k = std::stoi(argv[1]);
    D = std::stoi(argv[2]);
    bsize = std::stoi(argv[3]);
    if(k < 2) {
        std::cout << "k must be at least 2." << std::endl;
        return 0;
    }
    if(D <= k) {
        std::cout << "D must exceed k." << std::endl;
        return 0;
    }
    std::random_device rd;
    uint64_t seed = rd();
    std::mt19937_64 eng(seed);
    std::cout << "seed = " << seed << std::endl;
    std::string filename = "instance_" + std::to_string(k) + "_" + std::to_string(D) + "_" + std::to_string(bsize) + ".tsv";
    bitmatrix H = gallager(k,D,bsize,eng);
    bitvector v(H.num_cols());
    v.random(eng);
    save_max_xorsat_instance(H, v, filename);
    return 0;
}