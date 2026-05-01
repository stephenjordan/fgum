#include <iostream>
#include <random>
#include "xoropt.hpp"

int main(int argc, char *argv[]) {
    if(argc != 2 && argc != 3) {
        std::cout << "Usage: full_anneal instance.tsv <sweeps>" << std::endl;
        return 0;
    }
    xorsat_instance I;
    if(!I.load(argv[1])) return 0;
    std::random_device rd;
    uint64_t seed = rd();
    std::cout << "seed: " << seed << std::endl;
    std::mt19937_64 eng(seed);
    walker w(I);
    w.randomize(eng);
    int sweeps = 5000;
    if(argc == 3) {
        sweeps = atoi(argv[2]);
        if(sweeps < 1 || sweeps > 1E8) {
            std::cout << "Error: invalid sweep count " << sweeps << " reverting to 5000" << std::endl;
            sweeps = 5000;
        }
    }
    anneal(w, sweeps, 0.0, 3.0, eng);
    return 0;
}