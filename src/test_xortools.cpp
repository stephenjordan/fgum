#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "bitmatrix.hpp"
#include "utils.hpp"

//bitmatrix tests:---------------------------------------------------------------------------------

int non_invertible_test(std::string filename) {
    bool success, invertible;
    bitmatrix A,I;
    int failed = 0;
    success = (int)A.load_dense(filename);
    if(!success) { notify("noninvertiblefail1: " + filename); failed++; }
    else {
        I = A.inverse(invertible);
        if(invertible) { notify("noninvertiblefail2: " + filename); failed++; }
    }
    return failed;
}

int invertible_test(std::string testname, std::string inversename) {
    bool success, invertible;
    bitmatrix A, I, IR;
    int failed = 0;
    success = A.load_dense(testname);
    if(!success) { notify("invertiblefail1: " + testname); failed++; }
    else {
        I = A.inverse(invertible);
        if(!invertible) { notify("invertiblefail2: " + testname); failed++; }
        success = IR.load_dense(inversename);
        if(!success) { notify("invertiblefail3: " + testname); failed++; }
        else {
            if(IR != I) { notify("invertiblefail4: " + testname); failed++; }
        }
    }
    return failed;
}

int test_rank(std::string filename, size_t rows, size_t cols, int rank) {
    int failed = 0;
    bool success;
    bitmatrix M;
    success = M.load_dense(filename);
    if(!success) { notify("rankfail1:" + filename); failed++; }
    else {
        if(rows != M.num_rows()) { notify("rankfail2: " + filename); failed++; }
        if(cols != M.num_cols()) { notify("rankfail3: " + filename); failed++; }
        if(M.rank() != rank) {
            notify("rankfail4: " + filename);
            failed++;
        }
    }
    return failed;
}

int test_determinant(std::string filename, bool val) {
    int failed = 0;
    bool success;
    bitmatrix M;
    success = M.load_dense(filename);
    if(!success) { notify("determinantfail1: " + filename); failed++; }
    else {
        if(M.determinant() != val) { notify("determinantfail2 " + filename); failed++; }
    }
    return failed;
}

int inversion_tests() {
    int failed = 0;
    failed += non_invertible_test("data/inv_test0.txt");                  //inv_test0 is not invertible
    failed += invertible_test("data/inv_test1.txt", "data/inverse1.txt"); //inv_test1 is invertible
    failed += non_invertible_test("data/inv_test2.txt");                  //inv_test2 is not invertible
    failed += non_invertible_test("data/inv_test3.txt");                  //inv_test3 is not invertible
    failed += non_invertible_test("data/inv_test4.txt");                  //inv_test4 is not invertible
    failed += non_invertible_test("data/inv_test5.txt");                  //inv_test5 is not invertible
    failed += invertible_test("data/inv_test6.txt", "data/inverse6.txt"); //inv_test6 is invertible
    return failed;
}

int rank_tests() {
    int failed = 0;
    failed += test_rank("data/testB0.txt",4,5,3);
    failed += test_rank("data/testB1.txt",6,10,6);
    failed += test_rank("data/testB2.txt",70,60,60);
    failed += test_rank("data/testB3.txt",30,20,20);
    failed += test_rank("data/testB4.txt",131,131,129);
    return failed;
}

int determinant_tests() {
    int failed = 0;
    failed += test_determinant("data/testC0.txt", 1);
    failed += test_determinant("data/testC1.txt", 0);
    failed += test_determinant("data/testC2.txt", 1);
    failed += test_determinant("data/testC3.txt", 0);
    failed += test_determinant("data/testC4.txt", 0);
    failed += test_determinant("data/testC5.txt", 0);
    failed += test_determinant("data/testC6.txt", 1);
    return failed;
}

int multiply_tests() {
    int failed = 0;
    bool loaded;
    bitmatrix A,B,C,D;
    for(int num = 0; num <= 3; num++) {
        std::string filenameA = "data/testP" + std::to_string(num) + "A.txt";
        std::string filenameB = "data/testP" + std::to_string(num) + "B.txt";
        std::string filenameC = "data/testP" + std::to_string(num) + "C.txt";
        loaded = A.load_dense(filenameA);
        if(!loaded) {
            notify("multiplyfail1");
            failed++;
            continue;
        }
        loaded = B.load_dense(filenameB);
        if(!loaded) {
            notify("multiplyfail2");
            failed++;
            continue;
        }
        loaded = C.load_dense(filenameC);
        if(!loaded) {
            notify("multiplyfail3");
            failed++;
            continue;
        }
        D = A*B;
        if(D != C) { notify("multiplyfail4"); failed++; }
    }
    return failed;
}

int io_tests(std::mt19937_64 &eng) {
    bitmatrix B,C;
    bool loaded;
    int failed = 0;
    for(int i = 1; i < 5; i++) {
        bitmatrix A(15*i,22*i);
        A.random(eng);
        A.save_dense("data/densetest.txt");
        A.save_sparse("data/sparsetest.txt");
        loaded = B.load_sparse("data/sparsetest.txt");
        if(!loaded) { notify("iofail1"); failed++; }
        else if(A != B) { notify("iofail2"); failed++; }
        loaded = C.load_dense("data/densetest.txt");
        if(!loaded) { notify("iofail3"); failed++; }
        else if(A != C) { notify("iofail4"); failed++; }
        remove("data/densetest.txt");
        remove("data/sparsetest.txt");
    }
    return failed;
}

int sparse_tests(std::mt19937_64 &eng) {
    int failed = 0;
    std::uniform_int_distribution<uint64_t> distr;
    for(int bsize = 1; bsize < 5; bsize++) {
        bitmatrix A;
        A.random(32*bsize,28*bsize,eng);
        size_t ones = A.count_ones();
        std::vector <std::vector <size_t> > sparse_rows = A.to_sparse_rows();
        std::vector <std::vector <size_t> > sparse_cols = A.to_sparse_cols();
        bitmatrix C,D;
        C.from_sparse_rows(sparse_rows, A.num_rows(), A.num_cols());
        D.from_sparse_cols(sparse_cols, A.num_rows(), A.num_cols());
        if(C != A) { notify("sparsefail1"); failed++; }
        if(D != A) { notify("sparsefail2"); failed++; }
        size_t count = 0;
        for(size_t i = 0; i < sparse_rows.size(); i++) count += sparse_rows[i].size();
        if(count != ones) { notify("sparsefail3"); failed++; }
        count = 0;
        for(size_t i = 0; i < sparse_cols.size(); i++) count += sparse_cols[i].size();
        if(count != ones) { notify("sparsefail4"); failed++; }
        for(int rep = 0; rep < 20; rep++) { //swaps should not alter number of ones
            size_t r1 = distr(eng)%A.num_rows();
            size_t r2 = distr(eng)%A.num_rows();
            A.row_swap(r1,r2);
            r1 = distr(eng)%A.num_cols();
            r2 = distr(eng)%A.num_cols();
            A.col_swap(r1,r2);
        }
        if(A.count_ones() != ones) { notify("sparsefail5"); failed++; }
    }
    return failed;
}

int entrywise_tests(std::mt19937_64 &eng) {
    int failed = 0;
    for(size_t asize = 1; asize < 3; asize++) {
        for(size_t bsize = 1; bsize < 3; bsize++) {
            bitmatrix Rand,Zero,Ones;
            Rand.random(25*asize, 25*bsize, eng);
            Zero.zeros(25*asize,25*bsize);
            Ones = ~Zero;
            bitmatrix RandAndZero = Rand & Zero;
            bitmatrix RandOrZero = Rand | Zero;
            bitmatrix RandPlusZero = Rand + Zero;
            bitmatrix RandAndOnes = Rand & Ones;
            bitmatrix RandOrOnes = Rand | Ones;
            bitmatrix RandPlusOnes = Rand + Ones;
            bitmatrix ZeroAndOnes = Zero & Ones;
            bitmatrix ZeroOrOnes = Zero | Ones;
            bitmatrix ZeroPlusOnes = Zero + Ones;
            bitmatrix RandPlusRand = Rand + Rand;
            if(RandAndZero != Zero) { notify("entrywisefail1"); failed++; }
            if(RandOrZero != Rand) { notify("entrywisefail2"); failed++; }
            if(RandPlusZero != Rand) { notify("entrywisefail3"); failed++; }
            if(RandAndOnes != Rand) { notify("entrywisefail4"); failed++; }
            if(RandOrOnes != Ones) { notify("entrywisefail5"); failed++; }
            if(RandPlusOnes != ~Rand) { notify("entrywisefail6"); failed++; }
            if(ZeroAndOnes != Zero) { notify("entrywisefail7"); failed++; }
            if(ZeroOrOnes != Ones) { notify("entrywisefail8"); failed++; }
            if(RandPlusRand != Zero) { notify("entrywisefail9"); failed++; }
        }
    }
    return failed;
}

int singles_tests(std::mt19937_64 &eng) {
    int failed = 0;
    for(size_t asize = 1; asize < 3; asize++) {
        for(size_t bsize = 1; bsize < 3; bsize++) {
            bitmatrix rand;
            rand.random(25*asize, 25*bsize, eng);
            std::vector <bitmatrix> rows;
            std::vector <bitmatrix> cols;
            std::vector <bitmatrix> coltrans;
            for(size_t i = 0; i < rand.num_rows(); i++) rows.push_back(rand.row(i));
            for(size_t j = 0; j < rand.num_cols(); j++) {
                cols.push_back(rand.col(j));
                coltrans.push_back(rand.col_transpose(j));
            }
            bitmatrix from_rows, from_cols, from_coltrans;
            from_rows = rows[0];
            for(size_t i = 1; i < rows.size(); i++) from_rows.append_below(rows[i]);
            from_cols = cols[0];
            for(size_t j = 1; j < cols.size(); j++) from_cols.append_right(cols[j]);
            from_coltrans = coltrans[0];
            for(size_t j = 1; j < coltrans.size(); j++) from_coltrans.append_below(coltrans[j]);
            from_coltrans = from_coltrans.transpose();
            if(from_cols != rand) { notify("singlesfail1"); failed++; }
            if(from_rows != rand) { notify("singlesfail2"); failed++; }
            if(from_coltrans != rand) { notify("singlesfail3"); failed++; }
        }
    }
    return failed;
}

int op_tests(std::mt19937_64 &eng) {
    int failed = 0;
    std::uniform_int_distribution<uint64_t> distr;
    for(size_t asize = 1; asize < 3; asize++) {
        for(size_t bsize = 1; bsize < 3; bsize++) {
            bitmatrix rand;
            rand.random(25*asize, 25*bsize, eng);
            bitmatrix tmp;
            std::vector <size_t> i = random_subset<size_t>(rand.num_rows(), 3, eng);
            std::vector <size_t> j = random_subset<size_t>(rand.num_cols(), 3, eng);
            tmp = rand;
            tmp.row_swap(i[0],i[1]);
            tmp.row_swap(i[2], i[2]); //should do nothing
            if(tmp.row(i[0]) != rand.row(i[1])) { notify("opfail1"); failed++; }
            if(tmp.row(i[1]) != rand.row(i[0])) { notify("opfail2"); failed++; }
            tmp = rand;
            tmp.col_swap(j[0],j[1]);
            tmp.col_swap(j[2],j[2]); //should do nothing
            if(tmp.col(j[0]) != rand.col(j[1])) { notify("opfail3"); failed++; }
            if(tmp.col(j[1]) != rand.col(j[0])) { notify("opfail4"); failed++; }
            tmp = rand;
            tmp.row_add(i[0],i[1]);
            tmp.row_add(i[0],i[1]); //should cancel
            if(tmp != rand) { notify("opfail5"); failed++; }
            tmp = rand;
            tmp.col_add(j[0],j[1]);
            tmp.col_add(j[0],j[1]); //should cancel
            if(tmp != rand) { notify("opfail6"); failed++; }
            tmp = rand;
            tmp.row_add(i[0],i[1]);
            tmp.row_add(i[1],i[0]);
            tmp.row_add(i[0],i[1]);
            tmp.row_swap(i[0],i[1]); //3 cnots cancels a swap
            if(tmp != rand) { notify("opfail7"); failed++; }
            tmp = rand;
            tmp.col_add(j[0],j[1]);
            tmp.col_add(j[1],j[0]);
            tmp.col_add(j[0],j[1]);
            tmp.col_swap(j[0],j[1]); //3 cnots cancels a swap
            if(tmp != rand) { notify("opfail8"); failed++; }
            tmp = rand;
            tmp.row_add(i[0],i[0]); //adding a row to itself should yield all zeros
            tmp.col_add(j[0],j[0]); //adding a col to itself should yield all zeros
            bitmatrix R = tmp.row(i[0]);
            bitmatrix C = tmp.col(j[0]);
            if(R.count_ones() != 0) { notify("opfail9"); failed++; }
            if(C.count_ones() != 0) { notify("opfail10"); failed++; }
        }
    }
    return failed;
}

int dot_tests() {
    int failed = 0;
    bitmatrix I;
    int Z;
    bool Z2;
    for(size_t asize = 1; asize < 3; asize++) {
        I.identity(asize*30);
        for(size_t i = 0; i < asize*30; i++) {
            for(size_t j = 0; j < asize*30; j++) {
                Z = Z_dot(I.col(i), I.col(j));
                Z2 = Z2_dot(I.col(i), I.col(j));
                if(Z2 != (i == j)) { notify("dotfail1"); failed++; }
                if(Z != (i == j)) { notify("dotfail2"); failed++; }
                Z = Z_dot(I.row(i), I.row(j));
                Z2 = Z2_dot(I.row(i), I.row(j));
                if(Z2 != (i == j)) { notify("dotfail3"); failed++; }
                if(Z != (i == j)) { notify("dotfail4"); failed++; }
                Z = Z_dot(I.row(i), I.col(j));
                Z2 = Z2_dot(I.row(i), I.col(j));
                if(Z2 != (i == j)) { notify("dotfail5"); failed++; }
                if(Z != (i == j)) { notify("dotfail6"); failed++; }
                Z = Z_dot(I.col(i), I.row(j));
                Z2 = Z2_dot(I.col(i), I.row(j));
                if(Z2 != (i == j)) { notify("dotfail7"); failed++; }
                if(Z != (i == j)) { notify("dotfail8"); failed++; }
                Z = Z_dot(I.col_transpose(i), I.col_transpose(j));
                Z2 = Z2_dot(I.col_transpose(i), I.col_transpose(j));
                if(Z2 != (i == j)) { notify("dotfail9"); failed++; }
                if(Z != (i == j)) { notify("dotfail10"); failed++; }
            }
        }
        I.zeros(asize*31, asize*21);
        I = ~I;
        for(size_t i = 0; i < I.num_rows(); i++) {
            for(size_t j = 0; j < I.num_rows(); j++) {
                if(Z_dot(I.row(i), I.row(j)) != I.num_cols()) { notify("dotfail11"); failed++; }
                if(Z2_dot(I.row(i), I.row(j)) != I.num_cols()%2) { notify("dotfail12"); failed++; }
            }
        }
        for(size_t i = 0; i < I.num_cols(); i++) {
            for(size_t j = 0; j < I.num_cols(); j++) {
                if(Z_dot(I.col(i), I.col(j)) != I.num_rows()) { notify("dotfail13"); failed++; }
                if(Z2_dot(I.col(i), I.col(j)) != I.num_rows()%2) { notify("dotfail14"); failed++; }
            }
        }
    }
    return failed;
}

int get_set_tests(std::mt19937_64 &eng) {
    int failed = 0;
    bitmatrix M,N1,N2;
    for(size_t asize = 1; asize < 3; asize++) {
        for(size_t bsize = 1; bsize < 3; bsize++) {
            M.random(asize*25, bsize*25, eng);
            N1.zeros(asize*25, bsize*25);
            N2.zeros(asize*25, bsize*25);
            for(size_t i = 0; i < M.num_rows(); i++) {
                for(size_t j = 0; j < M.num_cols(); j++) {
                    N1.set(i, j, M.get(i,j));
                    if(M.get(i,j)) N2.flip(i,j);
                }
            }
            if(M != N1) { notify("getsetfail1"); failed++; }
            if(M != N2) { notify("getsetfail2"); failed++; }
        }
    }
    return failed;
}

int symmetry_tests(std::mt19937_64 &eng) {
    int failed = 0;
    bitmatrix I,R;
    for(size_t asize = 1; asize < 3; asize++) {
        I.identity(asize*25);
        R.random(asize*25, asize*25, eng);
        R += R.transpose();
        if(I.transpose() != I) { notify("symmetryfail1"); failed++; }
        if(R.transpose() != R) { notify("symmetryfail2"); failed++; }
    }
    return failed;
}

int slice_tests(std::mt19937_64 &eng) {
    int failed = 0;
    bitmatrix R1,R2,piece1,piece2;
    std::uniform_int_distribution<uint64_t> distr;
    for(size_t asize = 1; asize < 3; asize++) {
        for(size_t bsize = 1; bsize < 3; bsize++) {
            R1.random(25*asize, 25*bsize, eng);
            size_t row_slicepoint = distr(eng)%(R1.num_rows()-2)+1; //we don't want an empty row-slice
            size_t col_slicepoint = distr(eng)%(R1.num_cols()-2)+1; //we don't want an empty col-slice
            piece1 = R1.row_range(0, row_slicepoint);
            piece2 = R1.row_range(row_slicepoint, R1.num_rows());
            R2 = vstack(piece1, piece2);
            if(R2 != R1) { notify("slicefail1"); failed++; }
            R2 = piece1;
            R2.append_below(piece2);
            if(R2 != R1) { notify("slicefail2"); failed++; }
            piece1 = R1.col_range(0, col_slicepoint);
            piece2 = R1.col_range(col_slicepoint, R1.num_cols());
            R2 = hstack(piece1, piece2);
            if(R2 != R1) { notify("slicefail3"); failed++; }
            R2 = piece1;
            R2.append_right(piece2);
            if(R2 != R1) { notify("slicefail4"); failed++; }
        }
    }
    return failed;
}

int tripleslice_tests(std::mt19937_64 &eng) {
    int failed = 0;
    bitmatrix R1,R2,piece1,piece2,piece3;
    for(size_t asize = 1; asize < 3; asize++) {
        for(size_t bsize = 1; bsize < 3; bsize++) {
            R1.random(25*asize, 25*bsize, eng);
            std::vector <size_t>row_slicepoints = random_subset<size_t>(R1.num_rows()-1, 2, eng, 1);
            std::vector <size_t>col_slicepoints = random_subset<size_t>(R1.num_cols()-1, 2, eng, 1);
            piece1 = R1.row_range(0, row_slicepoints[0]);
            piece2 = R1.row_range(row_slicepoints[0], row_slicepoints[1]);
            piece3 = R1.row_range(row_slicepoints[1], R1.num_rows());
            R2 = vstack(piece1, piece2);
            R2.append_below(piece3);
            if(R2 != R1) { notify("tripleslicefail1"); failed++; }
            piece1 = R1.col_range(0, col_slicepoints[0]);
            piece2 = R1.col_range(col_slicepoints[0], col_slicepoints[1]);
            piece3 = R1.col_range(col_slicepoints[1], R1.num_cols());
            R2 = hstack(piece1, piece2);
            R2.append_right(piece3);
            if(R2 != R1) { notify("tripleslicefail2"); failed++; }
        }
    }
    return failed;
}

//returns the index of the column containing the first one in row i, if it exists
//if the row is all zeros the returns the number of columns
size_t leading_entry(const bitmatrix &M, size_t i) {
    if(i >= M.num_rows()) {
        std::cout << "Error:index out of range in leading_entry" << std::endl;
        return 0;
    }
    size_t j;
    for(j = 0; j < M.num_cols(); j++) {
        if(M.get(i,j)) break;
    }
    return j;
}

bool row_echelon(const bitmatrix &M) {
    if(M.num_rows() == 1) return true;
    std::vector <size_t> leading_entries(M.num_rows());
    for(size_t row_index = 0; row_index < M.num_rows(); row_index++) leading_entries[row_index] = leading_entry(M,row_index);
    for(size_t row_index = 1; row_index < M.num_rows(); row_index++) {
        if(leading_entries[row_index] <= leading_entries[row_index-1] && leading_entries[row_index] != M.num_cols()) return false;
    }
    return true;
}

bool reduced_row_echelon(const bitmatrix &M) {
    if(M.num_rows() == 1) return true;
    std::vector <size_t> leading_entries(M.num_rows());
    for(size_t row_index = 0; row_index < M.num_rows(); row_index++) leading_entries[row_index] = leading_entry(M,row_index);
    for(size_t row_index = 1; row_index < M.num_rows(); row_index++) {
        if(leading_entries[row_index] <= leading_entries[row_index-1] && leading_entries[row_index] != M.num_cols()) return false;
    }
    for(size_t row_index = 1; row_index < M.num_rows(); row_index++) {
        if(leading_entries[row_index] < M.num_cols()) {
            bitmatrix col = M.col(leading_entries[row_index]);
            size_t num_ones = col.count_ones();
            if(num_ones != 1) return false;
        }
    }
    return true;
}

int echelon_tests(std::mt19937_64 &eng) {
    int failed = 0;
    bitmatrix M,E,RE;
    for(size_t asize = 1; asize < 3; asize++) {
        for(size_t bsize = 1; bsize < 3; bsize++) {
            M.random(25*asize,25*bsize,eng);
            E  = M.row_echelon_form();
            RE = M.reduced_row_echelon_form();
            if(!row_echelon(E)) { notify("echelonfail1"); failed++; }
            if(!reduced_row_echelon(RE)) { notify("echelonfail2"); failed++; }
        }
    }
    return failed;
}

int decomp_tests(std::mt19937_64 &eng) {
    int failed = 0;
    bitmatrix M,E,RE,G,R,RG,RR;
    for(size_t asize = 1; asize < 3; asize++) {
        for(size_t bsize = 1; bsize < 3; bsize++) {
            M.random(25*asize,25*bsize,eng);
            E  = M.row_echelon_form();
            RE = M.reduced_row_echelon_form();
            M.row_echelon_decomp(G,R);
            M.reduced_row_echelon_decomp(RG,RR);
            if(R != E) { notify("decompfail1"); failed++; }
            if(RR != RE) { notify("decompfail2"); failed++; }
            if(G * M != R) { notify("decompfail3"); failed++; }
            if(RG * M != RR) { notify("decompfail4"); failed++; }
        }
    }
    return failed;
}

int string_tests(std::mt19937_64 &eng) {
    int failed = 0;
    bitmatrix R1,R2,R3;
    std::string s;
    for(size_t asize = 1; asize < 3; asize++) {
        for(size_t bsize = 1; bsize < 3; bsize++) {
            R1.random(25*asize,25*bsize,eng);
            s = R1.to_str();
            R2.from_str(s, '\n');
            std::ostringstream ss;
            ss << R1;
            R3.from_str(ss.str(), '\n');         
            if(R1 != R2) { notify("stringfail1"); failed++; }
            if(R1 != R3) { notify("stringfail2"); failed++; }
        }
    }
    return failed;
}

int low_rank_tests(std::mt19937_64 &eng) {
    int failed = 0;
    bitmatrix R1,R2,R3;
    int rank;
    for(size_t asize = 1; asize < 3; asize++) {
        for(size_t bsize = 1; bsize < 3; bsize++) {
            R1.random(12*asize,25*bsize,eng);
            R2.zeros(12*asize, 25*bsize);
            rank = R1.rank();
            R3 = vstack(R1,R2);
            for(size_t row_index = 2; row_index < 12*asize; row_index++) {
                //fill in some of the bottom rows with sums of pairs of the top rows
                R3.row_add(row_index-1,row_index-1+12*asize);
                R3.row_add(row_index, row_index+12*asize);
                if(R3.rank() != rank) { notify("low_rank_fail"); failed++; }
            }
        }
    }
    return failed;
}

//failure = crash
void bounds_tests(std::mt19937_64 &eng) {
    std::ostringstream msg;
    notify("REDIRECT", msg);
    std::string errormsg = "Error: dimension mismatch in bitmatrix::+\n";
    errormsg += "Error: dimension mismatch in bitmatrix::*\n";
    errormsg += "Error: dimension mismatch in bitmatrix::&\n";
    errormsg += "Error: dimension mismatch in bitmatrix::|\n";
    errormsg += "Error: dimension mismatch in bitmatrix::append_right\n";
    errormsg += "Error: dimension mismatch in bitmatrix::append_below\n";
    errormsg += "Error: index out of bounds in bitmatrix::get\n";
    errormsg += "Error: index out of bounds in bitmatrix::set\n";
    errormsg += "Error: index out of bounds in bitmatrix::flip\n";
    errormsg += "Error: index out of bounds in bitmatrix::col_add\n";
    errormsg += "Error: index out of bounds in bitmatrix::row_add\n";
    errormsg += "Error: index out of bounds in bitmatrix::col_swap\n";
    errormsg += "Error: index out of bounds in bitmatrix::row_swap\n";
    errormsg += "Error: index out of bounds in bitmatrix::col_zero\n";
    errormsg += "Error: index out of bounds in bitmatrix::row_zero\n";
    errormsg += "Error: index out of range in bitmatrix::row\n";
    errormsg += "Error: index out of range in bitmatrix::col\n";
    bitmatrix R1,R2,R3;
    for(size_t asize = 1; asize < 3; asize++) {
        for(size_t bsize = 1; bsize < 3; bsize++) {
            R1.random(11*asize,12*bsize,eng);
            R2.zeros(13*asize, 14*bsize);
            R3 = R1 + R2;
            R3 = R1 * R2;
            R3 = R1 & R2;
            R3 = R1 | R2;
            R3 = hstack(R1,R2);
            R3 = vstack(R1,R2);
            R1.get(900,900);
            R1.get(900,0);
            R1.get(0,900);
            R2.set(900,900,1);
            R2.set(900,0,1);
            R2.set(0,900,1);
            R3.flip(800,800);
            R3.flip(800,0);
            R3.flip(0,800);
            R1.col_add(0,800);
            R1.col_add(800,0);
            R1.row_add(0,800);
            R1.row_add(800,0);
            R1.col_swap(0,800);
            R1.col_swap(800,0);
            R1.row_swap(0,800);
            R1.row_swap(800,0);
            R1.col_zero(800);
            R1.row_zero(800);
            R3 = R1.row(800);
            R3 = R1.col(800);
        }
    }
    notify("REDIRECT"); //back to std::cout
    if(msg.str() != errormsg) {
        std::cout << "Error message not as expected: " << std::endl;
        std::cout << msg.str() << std::endl;
        std::cout << " vs " << std::endl;
        std::cout << errormsg << std::endl;
    }
}

int zero_tests(std::mt19937_64 &eng) {
    int failed = 0;
    std::uniform_int_distribution<uint64_t> distr;
    bitmatrix R1,R2,R3;
    size_t i,j;
    for(size_t asize = 1; asize < 3; asize++) {
        for(size_t bsize = 1; bsize < 3; bsize++) {
            R1.random(35*asize,25*bsize,eng);
            R2 = R1;
            R3.zeros(1,25*bsize);
            i = distr(eng)%R2.num_rows();
            R2.row_zero(i);
            for(j = 0; j < R2.num_rows(); j++) {
                if(j != i && R2.row(j) != R1.row(j)) { notify("zero_fail1"); failed++; }
                if(j == i && R2.row(j) != R3) { notify("zero_fail2"); failed++; }
            }
            R2 = R1;
            R3.zeros(35*asize,1);
            i = distr(eng)%R2.num_cols();
            R2.col_zero(i);
            for(j = 0; j < R2.num_cols(); j++) {
                if(j != i && R2.col(j) != R1.col(j)) { notify("zero_fail3"); failed++; }
                if(j == i && R2.col(j) != R3) { notify("zero_fail4"); failed++; }
            }
        }
    }
    return failed;
}

void timing_tests(std::mt19937_64 &eng) {
    bitmatrix R1,R2,R3,R4,R5,R6,R7,R8;
    std::cout << "Gathering timing data. (All times in milliseconds.)" << std::endl;
    for(size_t asize = 1; asize < 5; asize++) {
        auto t1 = std::chrono::high_resolution_clock::now();
        R1.random(1000*asize, 1000*asize, eng);
        R2.random(1000*asize, 1000*asize, eng);
        auto t2 = std::chrono::high_resolution_clock::now();
        R1.determinant();
        R2.determinant();
        auto t3 = std::chrono::high_resolution_clock::now();
        R3 = R1 * R2;
        R4 = R2 * R1;
        auto t4 = std::chrono::high_resolution_clock::now();
        R5 = R1.simple_transpose();
        R6 = R2.simple_transpose();
        auto t5 = std::chrono::high_resolution_clock::now();
        R7 = R1.transpose();
        R8 = R2.transpose();
        auto t6 = std::chrono::high_resolution_clock::now();
        auto gentime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        auto dettime = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
        auto prodtime = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3);
        auto trantime = std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4);
        auto ftrantime = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5);
        std::cout << "n = " << 1000*asize << std::endl;
        std::cout << "two random generations: " << gentime.count() << std::endl;
        std::cout << "two determinants: " << dettime.count() << std::endl;
        std::cout << "two matrix-matrix products: " << prodtime.count() << std::endl;
        std::cout << "two transposes: " << trantime.count() << std::endl;
        std::cout << "two fast transposes: " << ftrantime.count() << std::endl;
    }
}

bool is_identity(const std::vector<size_t> &p) {
    for(size_t i = 0; i < p.size(); i++) if(p[i] != i) return false;
    return true;
}

int permutation_tests(std::mt19937_64 &eng) {
    int failed = 0;
    std::uniform_int_distribution<> dice(5,100);
    for(int trial = 0; trial < 50; trial++) {
        size_t n = dice(eng);
        bitmatrix B,C,D;
        B.identity(n);
        std::vector<size_t> perm(n);
        std::vector<size_t> perm_inverse(n);
        for(size_t i = 0; i < n; i++) perm[i] = i;
        while(is_identity(perm)) std::shuffle(perm.begin(), perm.end(), eng);
        for(size_t i = 0; i < n; i++) perm_inverse[perm[i]] = i;
        C = permute_rows(B, perm);
        D = permute_rows(C, perm_inverse);
        if(C == B) failed++;
        if(D != B) failed++;
        C = permute_columns(B, perm);
        D = permute_columns(C, perm_inverse);
        if(C == B) failed++;
        if(D != B) failed++;
    }
    return failed;
}

int weightk_tests(std::mt19937_64 &eng) {
    int failed = 0;
    std::vector<std::string> all;
    for(int trial = 0; trial < 100; trial++) {
        bitvector v = rand_weight_k(eng, 1000, 100);
        if(v.count() != 100) failed++;
        all.push_back(v.to_str());
    }
    if(!all_distinct(all)) failed++; //at n=1000, k=100 the likelihood is 10^-136
    return failed;
}

int submatrix_tests(std::mt19937_64 &eng) {
    int failed = 0;
    std::uniform_int_distribution<> dice(0,99);
    for(int trial = 0; trial < 100; trial++) {
        size_t start = 0;
        size_t end = 0;
        while(start == end) {
            start = dice(eng);
            end = dice(eng);
        }
        if(end < start) { //swap them
            size_t tmp = end;
            end = start;
            start = tmp;
        }
        bitmatrix A(100,100);
        A.random(eng);
        bitmatrix B = A.col_range(start, end);
        bitmatrix C = B.row_range(start, end);
        bitvector v(100);
        for(size_t i = start; i < end; i++) v.set(i,1);
        bitmatrix D = submatrix(A, v);
        if(D != C) failed++;
    }
    return failed;
}

int directsum_tests(std::mt19937_64 &eng) {
    int failed = 0;
    std::uniform_int_distribution<> dice(1,70);
    std::uniform_int_distribution<> blockcount(3,10);
    for(int trial = 0; trial < 50; trial++) {
        bitmatrix A(dice(eng),dice(eng));
        bitmatrix B(dice(eng),dice(eng));
        A.random(eng);
        B.random(eng);
        bitmatrix C = direct_sum(A,B);
        bitmatrix AR = block(C, 0, A.num_rows(), 0, A.num_cols());
        bitmatrix BR = block(C, A.num_rows(), A.num_rows() + B.num_rows(), A.num_cols(), A.num_cols() + B.num_cols());
        if(AR != A) failed++;
        if(BR != B) failed++;
    }
    for(int trial = 0; trial < 50; trial++) {
        int num_blocks = blockcount(eng);
        std::vector<bitmatrix> blocks;
        for(int i = 0; i < num_blocks; i++) {
            bitmatrix tmp(dice(eng),dice(eng));
            tmp.random(eng);
            blocks.push_back(tmp);
        }
        bitmatrix S = direct_sum(blocks);
        std::vector<bitmatrix> decomp;
        size_t row_start = 0;
        size_t col_start = 0;
        for(size_t i = 0; i < blocks.size(); i++) {
            bitmatrix bl = block(S, row_start, row_start+blocks[i].num_rows(), col_start, col_start+blocks[i].num_cols());
            decomp.push_back(bl);
            row_start += blocks[i].num_rows();
            col_start += blocks[i].num_cols();
        }
        for(size_t i = 0; i < blocks.size(); i++) {
            if(decomp[i] != blocks[i]) failed++;
        }
    }
    return failed;
}

int transpose_test(size_t rows, size_t cols, std::mt19937_64 &eng) {
    bitmatrix A;
    A.random(rows, cols, eng);
    bitmatrix AT = A.simple_transpose();
    bitmatrix ATF = A.transpose();
    if(AT != ATF) {
        notify("Error: fast transpose test failed.");
        return 1; 
    }
    return 0;
}

int transpose_tests(std::mt19937_64 &eng) {
    int failed = 0;
    failed += transpose_test(90, 100, eng);
    failed += transpose_test(5,5,eng);
    failed += transpose_test(1,90, eng);
    failed += transpose_test(90,1,eng);
    failed += transpose_test(5,1,eng);
    failed += transpose_test(1,5,eng);
    failed += transpose_test(1000,1000,eng);
    failed += transpose_test(64,64,eng);
    return failed;
}

int bitmatrix_tests(std::mt19937_64 &eng) {
    int failed;
    int total_failed = 0;
    failed = inversion_tests();
    total_failed += failed;
    if(failed == 0) std::cout << "All inversion tests passed." << std::endl;
    else std::cout << failed << " inversion tests failed." << std::endl;
    failed = rank_tests();
    total_failed += failed;
    if(failed == 0) std::cout << "All rank tests passed." << std::endl;
    else std::cout << failed << " rank tests failed." << std::endl;
    failed = determinant_tests();
    total_failed += failed;
    if(failed == 0) std::cout << "All determinant tests passed." << std::endl;
    else std::cout << failed << " determinant tests failed." << std::endl;
    failed = multiply_tests();
    total_failed += failed;
    if(failed == 0) std::cout << "All multiply tests passed." << std::endl;
    else std::cout << failed << " multiply tests failed." << std::endl;
    failed = io_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All io tests passed." << std::endl;
    else std::cout << failed << " io tests failed." << std::endl;
    failed = sparse_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All sparse tests passed." << std::endl;
    else std::cout << failed << " sparse tests failed." << std::endl;
    failed = entrywise_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All entrywise tests passed." << std::endl;
    else std::cout << failed << " entrywise tests failed." << std::endl;
    failed = singles_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All singles tests passed." << std::endl;
    else std::cout << failed << " singles tests failed." << std::endl;
    failed = op_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All op tests passed." << std::endl;
    else std::cout << failed << " op tests failed." << std::endl;
    failed = dot_tests();
    total_failed += failed;
    if(failed == 0) std::cout << "All dot tests passed." << std::endl;
    else std::cout << failed << " dot tests failed." << std::endl;
    failed = get_set_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All getset tests passed." << std::endl;
    else std::cout << failed << " getset tests failed." << std::endl;
    failed = symmetry_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All symmetry tests passed." << std::endl;
    else std::cout << failed << " symmetry tests failed." << std::endl;
    failed = slice_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All slice tests passed." << std::endl;
    else std::cout << failed << " slice tests failed." << std::endl;
    failed = tripleslice_tests(eng);
    if(failed == 0) std::cout << "All tripleslice tests passed." << std::endl;
    else std::cout << failed << " tripleslice tests failed." << std::endl;
    failed = echelon_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All echelon tests passed." << std::endl;
    else std::cout << failed << " echelon tests failed." << std::endl;
    failed = decomp_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All decomp tests passed." << std::endl;
    else std::cout << failed << " decomp tests failed." << std::endl;
    failed = string_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All string tests passed." << std::endl;
    else std::cout << failed << " string tests failed." << std::endl;
    failed = low_rank_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All low rank tests passed." << std::endl;
    else std::cout << failed << " low rank tests failed." << std::endl;
    failed = zero_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All zero tests passed." << std::endl;
    else std::cout << failed << " zero tests failed." << std::endl;
    failed = permutation_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All permutation tests passed." << std::endl;
    else std::cout << failed << " permutation tests failed." << std::endl;
    failed = weightk_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All weightk tests passed." << std::endl;
    else std::cout << failed << " weightk tests failed." << std::endl;
    failed = submatrix_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All submatrix tests passed." << std::endl;
    else std::cout << failed << " submatrix tests failed." << std::endl;
    failed = directsum_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All direct sum tests passed." << std::endl;
    else std::cout << failed << " direct sum tests failed." << std::endl;
    failed = transpose_tests(eng);
    if(failed == 0) std::cout << "All transpose tests passed." << std::endl;
    else std::cout << failed << " direct transpose tests failed." << std::endl;
    total_failed += failed;
    std::cout << "Starting bounds tests...";
    bounds_tests(eng);
    std::cout << " passed." << std::endl;
    timing_tests(eng);
    return total_failed;
}

//bitvector tests:-------------------------------------------------------------------------------------

int vec_get_set_tests(std::mt19937_64 &eng) {
    int failed = 0;
    bitvector V,N1,N2;
    for(size_t asize = 1; asize < 3; asize++) {
        V.random(asize*25, eng);
        N1.zeros(asize*25);
        N2.zeros(asize*25);
        for(size_t i = 0; i < V.size(); i++) {
            N1.set(i, V.get(i));
            if(V.get(i)) N2.flip(i);
        }
        if(V != N1) { notify("vec_getsetfail1"); failed++; }
        if(V != N2) { notify("vec_getsetfail2"); failed++; }
    }
    return failed;
}

int vec_io_tests(std::mt19937_64 &eng) {
    bitvector B,C;
    bool loaded;
    int failed = 0;
    for(int i = 1; i < 5; i++) {
        bitvector A(22*i);
        A.random(eng);
        A.save_dense("data/vec_densetest.txt");
        A.save_sparse("data/vec_sparsetest.txt");
        loaded = B.load_sparse("data/vec_sparsetest.txt");
        if(!loaded) { notify("vec_iofail1"); failed++; }
        else if(A != B) { notify("vec_iofail2"); failed++; }
        loaded = C.load_dense("data/vec_densetest.txt");
        if(!loaded) { notify("vec_iofail3"); failed++; }
        else if(A != C) { notify("vec_iofail4"); failed++; }
        remove("data/vec_densetest.txt");
        remove("data/vec_sparsetest.txt");
    }
    return failed;
}

int vec_sparse_tests(std::mt19937_64 &eng) {
    int failed = 0;
    std::uniform_int_distribution<uint64_t> distr;
    for(int bsize = 1; bsize < 5; bsize++) {
        bitvector A;
        A.random(32*bsize,eng);
        size_t ones = A.count();
        std::vector <size_t> one_list = A.to_sparse();
        bitvector B;
        B.from_sparse(one_list, A.size());
        if(B != A) { notify("vec_sparsefail1"); failed++; }
        if(one_list.size() != ones) { notify("vec_sparsefail2"); failed++; }
    }
    return failed;
}

int vec_entrywise_tests(std::mt19937_64 &eng) {
    int failed = 0;
    for(size_t asize = 1; asize < 3; asize++) {
        bitvector Rand,Zero,Ones;
        Rand.random(25*asize, eng);
        Zero.zeros(25*asize);
        Ones = ~Zero;
        bitvector RandAndZero = Rand & Zero;
        bitvector RandOrZero = Rand | Zero;
        bitvector RandPlusZero = Rand + Zero;
        bitvector RandAndOnes = Rand & Ones;
        bitvector RandOrOnes = Rand | Ones;
        bitvector RandPlusOnes = Rand + Ones;
        bitvector ZeroAndOnes = Zero & Ones;
        bitvector ZeroOrOnes = Zero | Ones;
        bitvector ZeroPlusOnes = Zero + Ones;
        bitvector RandPlusRand = Rand + Rand;
        if(RandAndZero != Zero) { notify("vec_entrywisefail1"); failed++; }
        if(RandOrZero != Rand) { notify("vec_entrywisefail2"); failed++; }
        if(RandPlusZero != Rand) { notify("vec_entrywisefail3"); failed++; }
        if(RandAndOnes != Rand) { notify("vec_entrywisefail4"); failed++; }
        if(RandOrOnes != Ones) { notify("vec_entrywisefail5"); failed++; }
        if(RandPlusOnes != ~Rand) { notify("vec_entrywisefail6"); failed++; }
        if(ZeroAndOnes != Zero) { notify("vec_entrywisefail7"); failed++; }
        if(ZeroOrOnes != Ones) { notify("vec_entrywisefail8"); failed++; }
        if(RandPlusRand != Zero) { notify("vec_entrywisefail9"); failed++; }
    }
    return failed;
}

int vec_dot_tests() {
    int failed = 0;
    int Z;
    bool Z2;
    for(size_t asize = 1; asize < 3; asize++) {
        std::vector<bitvector> I(asize*30);
        for(size_t i = 0; i < asize*30; i++) {
            I[i].zeros(asize*30);
            I[i].set(i,1);
        }
        for(size_t i = 0; i < asize*30; i++) {
            for(size_t j = 0; j < asize*30; j++) {
                Z = Z_dot(I[i], I[j]);
                Z2 = Z2_dot(I[i], I[j]);
                if(Z2 != (i == j)) { notify("vec_dotfail1"); failed++; }
                if(Z != (i == j)) { notify("vec_dotfail2"); failed++; }
            }
        }
    }
    for(size_t asize = 1; asize < 3; asize++) {
        bitvector J,K;
        J.zeros(asize*31);
        K = ~J;
        if(Z_dot(J,J) != 0) { notify("vec_dotfail3"); failed++; }
        if(Z2_dot(J,J) != 0) { notify("vec_dotfail4"); failed++; }
        if(Z_dot(K,J) != 0) { notify("vec_dotfail5"); failed++; }
        if(Z2_dot(K,J) != 0) { notify("vec_dotfail6"); failed++; }
        if(Z_dot(J,K) != 0) { notify("vec_dotfail7"); failed++; }
        if(Z2_dot(J,K) != 0) { notify("vec_dotfail8"); failed++; }
        if(Z_dot(K,K) != K.size()) { notify("vec_dotfail9"); failed++; }
        if(Z2_dot(K,K) != K.size()%2) { notify("vec_dotfail10"); failed++; }
    }
    return failed;
}

int vec_slice_tests(std::mt19937_64 &eng) {
    int failed = 0;
    bitvector R1,R2,piece1,piece2,piece3;
    std::uniform_int_distribution<uint64_t> distr;
    //slice into two pieces and put back together
    for(size_t asize = 1; asize < 3; asize++) {
        R1.random(25*asize, eng);
        size_t row_slicepoint = distr(eng)%(R1.size()-2)+1; //we don't want an empty row-slice
        piece1 = R1.range(0, row_slicepoint);
        piece2 = R1.range(row_slicepoint, R1.size());
        R2 = hstack(piece1, piece2);
        if(R2 != R1) { notify("vec_slicefail1"); failed++; }
    }
    //slice into three pieces and put back together
    for(size_t asize = 1; asize < 3; asize++) {
        R1.random(25*asize, eng);
        std::vector <size_t>row_slicepoints = random_subset<size_t>(R1.size()-1, 2, eng, 1);
        piece1 = R1.range(0, row_slicepoints[0]);
        piece2 = R1.range(row_slicepoints[0], row_slicepoints[1]);
        piece3 = R1.range(row_slicepoints[1], R1.size());
        R2 = hstack(piece1, piece2);
        R2.append(piece3);
        if(R2 != R1) { notify("vec_slicefail2"); failed++; }
    }
    return failed;
}

int vec_string_tests(std::mt19937_64 &eng) {
    int failed = 0;
    bitvector R1,R2,R3;
    std::string s;
    for(size_t asize = 1; asize < 3; asize++) {
        R1.random(25*asize,eng);
        s = R1.to_str();
        R2.from_str(s);
        std::ostringstream ss;
        ss << R1;
        R3.from_str(ss.str());         
        if(R1 != R2) { notify("vec_stringfail1"); failed++; }
        if(R1 != R3) { notify("vec_stringfail2"); failed++; }
    }
    return failed;
}

int vec_pushback_tests(std::mt19937_64 &eng) {
    int failed = 0;
    bitvector R1;
    std::string s;
    for(size_t asize = 1; asize < 3; asize++) {
        R1.random(25*asize,eng);
        bitvector R2;
        for(size_t i = 0; i < R1.size(); i++) R2.push_back(R1.get(i));
        if(R2 != R1) { notify("vec_pushbackfail1"); failed++; }
        size_t count = R2.count();
        for(size_t i = 0; i < 25*asize; i++) R2.push_back(1);
        for(size_t i = 0; i < 25*asize; i++) R2.push_back(0);
        for(size_t i = 0; i < 25*asize; i++) R2.push_back(true);
        for(size_t i = 0; i < 25*asize; i++) R2.push_back(false);
        if(R2.count() != count + 50*asize) {notify("vec_pushbackfail2"); failed++; }
    }
    return failed;
}

int mixed_tests(std::mt19937_64 &eng) {
    int failed = 0;
    bitvector V1,V2,V3;
    bitmatrix M;
    for(size_t asize = 1; asize < 3; asize++) {
        V1.random(25*asize,eng);
        M.identity(25*asize);
        V2 = V1 * M;
        if(V2 != V1) { notify("mixed_fail1"); failed++; }
        V2 = M * V1;
        if(V2 != V1) { notify("mixed_fail2"); failed++; }
        M.zeros(20*asize, 25*asize);
        V2 = M * V1;
        V3.zeros(20*asize);
        if(V3 != V2) { notify("mixed_fail3"); failed++; }
        M.zeros(30*asize, 25*asize);
        V2 = M * V1;
        V3.zeros(30*asize);
        if(V3 != V2) { notify("mixed_fail4"); failed++; }
        M.zeros(25*asize, 20*asize);
        V2 = V1 * M;
        V3.zeros(20*asize);
        if(V3 != V2) { notify("mixed_fail5"); failed++; }
        M.zeros(25*asize, 30*asize);
        V2 = V1 * M;
        V3.zeros(30*asize);
        if(V3 != V2) { notify("mixed_fail6"); failed++; }
        M.zeros(25*asize, 25*asize);
        M = ~M; //all ones matrix
        V2 = M * V1; //V2 should be all ones (if V1 has odd Hamming weight) or all zeros (if V1 has even Hamming weight)
        if(V2.count() != 0 && V2.count() != V2.size()) { notify("mixed_fail5"); failed++; }
    }
    return failed;
}

int bigmix_tests(std::mt19937_64 &eng) {
    int failed = 0;
    bitvector V1, V2, V3;
    bitmatrix M1,Z,R,C;
    for(size_t asize = 1; asize < 3; asize++) {
        V1.random(25*asize,eng);
        R = V1.to_row_matrix();
        C = V1.to_col_matrix();
        V2 = R.row_vec(0);
        V3 = C.col_vec(0);
        if(V1 != V2) { notify("bigmix_fail1"); failed++; }
        if(V1 != V3) { notify("bigmix_fail2"); failed++; }
        M1.random(25*asize, 12*asize, eng);
        Z.zeros(25*asize, 12*asize);
        if(M1 * true != M1) { notify("bigmix_fail3"); failed++; }
        if(true * M1 != M1) { notify("bigmix_fail4"); failed++; }
        if(M1 * false != Z) { notify("bigmix_fail5"); failed++; }
        if(false * M1 != Z) { notify("bigmix_fail6"); failed++; }
        V2.zeros(25*asize);
        if(V1 * true != V1) { notify("bigmix_fail7"); failed++; }
        if(true * V1 != V1) { notify("bigmix_fail8"); failed++; }
        if(V1 * false != V2) { notify("bigmix_fail9"); failed++; }
        if(false * V1 != V2) { notify("bigmix_fail10"); failed++; }
    }
    return failed;
}

//failure = crash
void vec_bounds_tests(std::mt19937_64 &eng) {
    std::ostringstream msg;
    notify("REDIRECT", msg);
    std::string errormsg = "Error: dimension mismatch in bitvector::+\n";
    errormsg += "Error: dimension mismatch in bitvector::&\n";
    errormsg += "Error: dimension mismatch in bitvector::|\n";
    errormsg += "Error: index out of bounds in bitvector::get\n";
    errormsg += "Error: index out of bounds in bitvector::set\n";
    errormsg += "Error: index out of bounds in bitvector::flip\n";
    errormsg += "Error: dimension mismatch in bitvector dot product\n";
    errormsg += "Error: index out of range in bitmatrix::row_vec\n";
    errormsg += "Error: index out of range in bitmatrix::col_vec\n";
    errormsg += "Error: dimension mismatch in matrix times vector\n";
    errormsg += "Error: dimension mismatch in vector times matrix\n";
    bitvector V1,V2,V3;
    bitmatrix M;
    for(size_t asize = 1; asize < 3; asize++) {
        M.identity(35*asize);
        V1.random(11*asize, eng);
        V2.zeros(13*asize);
        V3 = V1 + V2;
        V3 = V1 & V2;
        V3 = V1 | V2;
        V1.get(900);
        V2.set(900,1);
        V3.flip(800);
        Z_dot(V1,V2);
        Z2_dot(V1,V2);
        M.row_vec(800);
        M.col_vec(800);
        V3 = M * V1;
        V3 = V1 * M;
    }
    notify("REDIRECT"); //back to std::cout
    if(msg.str() != errormsg) {
        std::cout << "Error message not as expected: " << std::endl;
        std::cout << msg.str() << std::endl;
        std::cout << " vs " << std::endl;
        std::cout << errormsg << std::endl;
    }
}

int append_tests(std::mt19937_64 &eng) {
    int failed = 0;
    bitmatrix M1, M2, M3, M4;
    bitvector V1,V2, V3;
    V1.random(300, eng);
    V2.random(300, eng);
    V3.random(300, eng);
    M1.append_below(V1);
    M1.append_below(V2);
    M1.append_below(V3);
    M2.append_right(V1);
    M2.append_right(V2);
    M2.append_right(V3);
    M3.append_below(M1);
    M3.append_below(M1);
    M4.append_right(M2);
    M4.append_right(M2);
    if(M1.row_vec(0) != V1) { notify("append_tests fail1"); failed++; }
    if(M1.row_vec(1) != V2) { notify("append_tests fail2"); failed++; }
    if(M1.row_vec(2) != V3) { notify("append_tests fail3"); failed++; }
    if(M2.col_vec(0) != V1) { notify("append_tests fail4"); failed++; }
    if(M2.col_vec(1) != V2) { notify("append_tests fail5"); failed++; }
    if(M2.col_vec(2) != V3) { notify("append_tests fail6"); failed++; }
    if(M3.row_vec(0) != V1) { notify("append_tests fail7"); failed++; }
    if(M3.row_vec(1) != V2) { notify("append_tests fail8"); failed++; }
    if(M3.row_vec(2) != V3) { notify("append_tests fail9"); failed++; }
    if(M3.row_vec(3) != V1) { notify("append_tests fail10"); failed++; }
    if(M3.row_vec(4) != V2) { notify("append_tests fail11"); failed++; }
    if(M3.row_vec(5) != V3) { notify("append_tests fail12"); failed++; }
    if(M4.col_vec(0) != V1) { notify("append_tests fail13"); failed++; }
    if(M4.col_vec(1) != V2) { notify("append_tests fail14"); failed++; }
    if(M4.col_vec(2) != V3) { notify("append_tests fail15"); failed++; }
    if(M4.col_vec(3) != V1) { notify("append_tests fail16"); failed++; }
    if(M4.col_vec(4) != V2) { notify("append_tests fail17"); failed++; }
    if(M4.col_vec(5) != V3) { notify("append_tests fail18"); failed++; }
    return failed;
}

int complement_tests(std::mt19937_64 &eng) {
    int failed = 0;
    std::uniform_int_distribution<> dice(5,100);
    std::bernoulli_distribution coin(0.5);
    for(int trial = 0; trial < 100; trial++) {
        size_t n = dice(eng);
        bitvector v1(n);
        bitvector v2(n);
        for(size_t i = 0; i < n; i++) {
            bool val = coin(eng);
            v1.set(i,val);
            v2.set(i,!val);
        }
        bitvector v3 = v1.complement();
        if(v3 != v2) failed++;
    }
    return failed;
}

size_t simplecount(const bitvector &v, bool val) {
    for(size_t i = 0; i < v.size(); i++) if(v.get(i) != val) return i;
    return v.size();
}

int lsb_tests(std::mt19937_64 &eng) {
    int failed = 0;
    std::uniform_int_distribution<> dice_n(5,200);
    for(int trial = 0; trial < 20; trial++) {
        int n = dice_n(eng);
        std::uniform_int_distribution<> dice_k(1,n-1);
        size_t k = dice_k(eng);
        bitvector u,v;
        u.random(n, eng);
        v = u;
        v.set_lsb_zeros(n);
        if(v.count() != 0) { notify("lsb_test fail1"); failed++; }
        v.set_lsb_ones(n);
        if(v.count() != v.size()) { notify("lsb_test fail2"); std::cout << v.to_str() << std::endl; failed++; }
        v = u;
        v.set_lsb_zeros(0);
        if(v != u) { notify("lsb_test_fail3"); failed++; }
        v = u;
        v.set_lsb_ones(0);
        if(v != u) { notify("lsb_test fail4"); failed++; }
        v = u;
        v.set_lsb_zeros(k);
        bool good = true;
        for(size_t i = 0; i < v.size(); i++) {
            if(i < k && v.get(i) != 0) good = false;
            if(i >= k && v.get(i) != u.get(i)) good = false;
        }
        if(!good) { notify("lsb_test fail5"); failed++; }
        v = u;
        v.set_lsb_ones(k);
        good = true;
        for(size_t i = 0; i < v.size(); i++) {
            if(i < k && v.get(i) != 1) good = false;
            if(i >= k && v.get(i) != u.get(i)) good = false;
        }
        if(!good) { notify("lsb_test fail6"); failed++; }
        v = u;
        if(v.get_lsb_zeros() != simplecount(v, 0)) { notify("lsb_test fail7"); failed++; }
        v.set_lsb_zeros(n);
        if(v.get_lsb_zeros() != simplecount(v, 0)) { notify("lsb_test fail8"); failed++; }
        v.set_lsb_ones(n);
        if(v.get_lsb_zeros() != simplecount(v, 0)) { notify("lsb_test fail9"); failed++; }
    }
    return failed;
}

int bitvector_tests(std::mt19937_64 &eng) {
    int failed;
    int total_failed = 0;
    failed = vec_get_set_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All vector getset tests passed." << std::endl;
    else std::cout << failed << " vector getset tests failed." << std::endl;
    failed = vec_io_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All vector io tests passed." << std::endl;
    else std::cout << failed << " vector io tests failed." << std::endl;
    failed = vec_sparse_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All vector sparse tests passed." << std::endl;
    else std::cout << failed << " vector sparse tests failed." << std::endl;
    failed = vec_entrywise_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All vector entrywise tests passed." << std::endl;
    else std::cout << failed << " vector entrywise tests failed." << std::endl;
    failed = vec_dot_tests();
    total_failed += failed;
    if(failed == 0) std::cout << "All vector dot tests passed." << std::endl;
    else std::cout << failed << " vector dot tests failed." << std::endl;
    failed = vec_slice_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All vector slice tests passed." << std::endl;
    else std::cout << failed << " vector slice tests failed." << std::endl;
    failed = vec_string_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All vector string tests passed." << std::endl;
    else std::cout << failed << " vector string tests failed." << std::endl;
    failed = vec_pushback_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All vector pushback tests passed." << std::endl;
    else std::cout << failed << " vector pushback tests failed." << std::endl;
    failed = mixed_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All mixed tests passed." << std::endl;
    else std::cout << failed << " mixed tests failed." << std::endl;
    failed = bigmix_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All bigmix tests passed." << std::endl;
    else std::cout << failed << " bigmix tests failed." << std::endl;
    failed = append_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All append tests passed." << std::endl;
    else std::cout << failed << " append tests failed." << std::endl;
    failed = complement_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All complement tests passed." << std::endl;
    else std::cout << failed << " complement tests failed." << std::endl;
    failed = lsb_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All lsb tests passed." << std::endl;
    else std::cout << failed << " lsb tests failed." << std::endl;
    std::cout << "Starting vector bounds tests...";
    vec_bounds_tests(eng);
    std::cout << " passed." << std::endl;
    return total_failed;
}

//utils tests:-------------------------------------------------------------------------------------

int utils_random_subset_tests(std::mt19937_64 &eng) {
    int failed = 0;
    std::uniform_int_distribution<uint64_t> distr;
    for(size_t asize = 1; asize < 3; asize++) {
        size_t k = distr(eng)%10;
        size_t n = 50*asize;
        std::vector <size_t> v = random_subset<size_t>(n, k, eng);
        if(v.size() != k) { notify("random_subset_fail1"); failed++; }
        for(size_t i = 0; i < v.size(); i++) if(v[i] < 0 || v[i] >= n) { notify("random_subset_fail2"); failed++; }
        std::sort(v.begin(), v.end());
        for(size_t i = 1; i < v.size(); i++) if(v[i] == v[i-1]) { notify("random_subset_fail2"); failed++; }
    }
    return failed;
}

std::string detokenize(std::vector <std::string> tokens, char separator) {
    bool first = true;
    std::string returnval;
    for(size_t i = 0; i < tokens.size(); i++) {
        if(!first) returnval += separator;
        returnval += tokens[i];
        first = false;
    }
    return returnval;
}

std::vector<std::string> pop_all(std::string line, char separator) {
    std::vector<std::string> returnval;
    std::string token;
    std::string s = line;
    while(s.length() > 0) {
        token = poptoken(s, separator);
        if(token.length() > 0) returnval.push_back(token);
    }
    return returnval;
}

std::vector<std::string> pop_all(std::string line, std::string separators) {
    std::vector<std::string> returnval;
    std::string token;
    std::string s = line;
    while(s.length() > 0) {
        token = poptoken(s, separators);
        if(token.length() > 0) returnval.push_back(token);
    }
    return returnval;
}

int utils_tokenize_tests() {
    int failed = 0;
    std::string s1 = "This is a string that tests the capabilities of the tokenizer from utils.";
    std::string s2 = "  This is    a more    challenging string for the tokenizer.  ";
    std::string s3 = "This is a more challenging string for the tokenizer.";
    std::string s4 = " This__is another challenge_ _ for the tokenizer.";
    std::string s5 = "This is another challenge for the tokenizer.";
    std::string reconstituted;
    std::vector <std::string> tokens;
    tokens = tokenize(s1, " ");
    reconstituted = detokenize(tokens, ' ');
    if(reconstituted != s1) { notify("tokenize_fail1"); failed++; }
    tokens = pop_all(s1, " ");
    reconstituted = detokenize(tokens, ' ');
    if(reconstituted != s1) { notify("poptoken_fail1"); failed++; }
    tokens = tokenize(s2, " ");
    reconstituted = detokenize(tokens, ' ');
    if(reconstituted != s3) { notify("tokenize_fail2"); failed++; }
    tokens = pop_all(s2, " ");
    reconstituted = detokenize(tokens, ' ');
    if(reconstituted != s3) { notify("poptoken_fail2"); failed++; }
    tokens = tokenize(s1, ' ');
    reconstituted = detokenize(tokens, ' ');
    if(reconstituted != s1) { notify("tokenize_fail3"); failed++; }
    tokens = pop_all(s1, ' ');
    reconstituted = detokenize(tokens, ' ');
    if(reconstituted != s1) { notify("poptoken_fail3"); failed++; }
    tokens = tokenize(s2, ' ');
    reconstituted = detokenize(tokens, ' ');
    if(reconstituted != s3) { notify("tokenize_fail4"); failed++; }
    tokens = pop_all(s2, ' ');
    reconstituted = detokenize(tokens, ' ');
    if(reconstituted != s3) { notify("poptoken_fail4"); failed++; }
    tokens = tokenize(s4, " _");
    reconstituted = detokenize(tokens, ' ');
    if(reconstituted != s5) { notify("tokenize_fail5"); failed++; }
    tokens = pop_all(s4, " _");
    reconstituted = detokenize(tokens, ' ');
    if(reconstituted != s5) { notify("poptoken_fail5"); failed++; }
    return failed;
}

int utils_approx_tests(std::mt19937_64 &eng) {
    int failed = 0;
    std::normal_distribution<double> big_normal(0.0,1.0E6);
    std::normal_distribution<double> small_normal(0.0,0.01);
    std::uniform_int_distribution<size_t> index_uniform;
    double target;
    for(size_t asize = 1; asize < 3; asize++) {
        size_t n = 50*asize;
        for(size_t rep = 0; rep < 10; rep++) {
            std::vector <double> list;
            for(size_t i = 0; i < n; i++) list.push_back(big_normal(eng));
            size_t rand_index = index_uniform(eng)%list.size();
            target = list[rand_index] + small_normal(eng);
            bool result = approx_present(list, target, 0.1);
            int location = approx_find(list, target, 0.1);
            if(!result) {
                std::cout << "approxfail1: target " << target << " not found amongst:" << std::endl;
                for(size_t i = 0; i < list.size(); i++) std::cout << list[i] << std::endl;
                failed++;
            }
            if(fabs(target - list[location]) > 0.1) {
                std::cout << "approxfail2: target " << target << " found at list[" << location << "] = " << list[location];
                std::cout << "instead of list[" << rand_index << "] = " << list[rand_index] << std::endl;
                for(size_t i = 0; i < list.size(); i++) std::cout << list[i] << std::endl;
                failed++;
            }
            target = big_normal(eng) + 1.0E7;
            result = approx_present(list, target, 0.1);
            if(result) {
                std::cout << "approxfail3: target " << target << " found amongst:" << std::endl;
                for(size_t i = 0; i < list.size(); i++) std::cout << list[i] << std::endl;
                failed++;
            }
            std::ostringstream msg;
            notify("CLEAR");
            notify("REDIRECT", msg);
            location = approx_find(list, target, 0.1);
            if(msg.str() != "Error: approximate value not found.\n") {
                std::cout << "msg = " << msg.str() << std::endl;
                std::cout << "approxfail4: target " << target << " found at location " << location << " amongst:" << std::endl;
                for(size_t i = 0; i < list.size(); i++) std::cout << i << ": " << list[i] << std::endl;
                failed++;
            }
            notify("REDIRECT");
        }
    }
    return failed;
}

int utils_io_tests(std::mt19937_64 &eng) {
    int failed = 0;
    std::normal_distribution<double> big_normal(0.0,1.0E6);
    std::normal_distribution<float> fbig_normal(0.0,1.0E6);
    std::uniform_int_distribution<long int> uniform;
    std::uniform_int_distribution<unsigned int> uuniform;
    std::string s;
    double x, xr, tolerance;
    float fx, fxr;
    long int y, yr;
    unsigned int uy, uyr;
    for(size_t num_digits = 2; num_digits < 6; num_digits++) {
        tolerance = pow(10,1-num_digits);
        for(size_t trial = 0; trial < 100; trial++) {
            x = big_normal(eng);
            s = to_digits(x, num_digits);
            xr = stod(s);
            if(fabs(x-xr) > tolerance) { notify("utils_io_fail1"); failed++; }
            s = fwf(x, num_digits);
            xr = stod(s);
            if(error_frac(x, xr) > 0.01) { notify("utils_io_fail2"); }
            fx = fbig_normal(eng);
            s = to_digits(fx, num_digits);
            fxr = stof(s);
            if(fabs(fx-fxr) > tolerance) { notify("utils_io_fail3"); failed++; }
            s = fwf(fx, num_digits);
            fxr = stod(s);
            if(error_frac(fx, fxr) > 0.01) { notify("utils_io_fail4"); }
        }
    }
    if(trailing_zeros("001.0227000") != 3) { notify("utils_io_fail5"); failed++; }
    if(trailing_zeros("-7.0227") != 0) { notify("utils_io_fail6"); failed++; }
    if(trailing_zeros("-100000.000000000") != 9) { notify("utils_io_fail7"); failed++; }
    if(trailing_zeros("3.12300007") != 0) { notify("utils_io_fail8"); failed++; }
    if(trim("001.0227000") != "001.0227000") { notify("utils_io_fail9"); failed++; }
    if(trim("-7.0227") != "-7.0227") { notify("utils_io_fail10"); failed++; }
    if(trim("-100000.000000000") != "-100000.000") { notify("utils_io_fail11"); failed++; }
    if(trim("3.12300007") != "3.12300007") { notify("utils_io_fail12"); failed++; }
    if(trim("-3.1230000") != "-3.123000") { notify("utils_io_fail13"); failed++; }
    for(size_t trial = 0; trial < 100; trial++) {
        y = uniform(eng);
        s = fwi(y);
        yr = stol(s);
        if(y != yr) { notify("utils_io_fail14"); failed++; }
        uy = uuniform(eng);
        s = fwi(uy);
        uyr = (unsigned int)stol(s);
        if(uy != uyr) { notify("utils_io_fail15"); failed++; }
    }
    if(fws("dkj", 5) != "dkj  ") { notify("utils_io_fail16"); failed++; }
    if(fws("  ()$HG(JHHHFdkj", 2) != "  ()$HG(JHHHFdkj") { notify("utils_io_fail17"); failed++; }
    if(fws("  ()$HG(JHHHFdkj", 22) != "  ()$HG(JHHHFdkj      ") { notify("utils_io_fail18"); failed++; }
    return failed;
}

int utils_math_tests() {
    int failed = 0;
    if(!approx_equal(logfactorial(5), 4.7874917)) { notify("utils_math_fail1"); failed++; }
    if(!approx_equal(logfactorial(10), 15.1044126)) { notify("utils_math_fail2"); failed++; }
    if(!approx_equal(logfactorial(256), 1167.2572786)) { notify("utils_math_fail3"); failed++; }
    if(error_frac(logfactorial(1272628), 1.6616196E7) > 1.0E-6) { notify("utils_math_fail4"); failed++; }
    if(!approx_equal(log_binomial(5,3), 2.3025851)) { notify("utils_math_fail5"); failed++; }
    if(!approx_equal(log_binomial(100,1), 4.6051702)) { notify("utils_math_fail6"); failed++; }
    if(!approx_equal(log_binomial(100,99), 4.6051702)) { notify("utils_math_fail7"); failed++; }
    if(!approx_equal(log_binomial(100,50), 66.7838417)) { notify("utils_math_fail8"); failed++; }
    if(!approx_equal(log_binomial(123450,50000), 83322.1312564)) { notify("utils_math_fail9"); failed++; }
    return failed;
}

int utils_vec_tests(std::mt19937_64 &eng) {
    int failed = 0;
    std::vector <std::string> all;
    std::vector <int> moduli;
    std::uniform_int_distribution<size_t> distr;
    std::vector <int> f;
    for(int length = 1; length <= 5; length++) {
        moduli.clear();
        all.clear();
        f.resize(length);
        for(int i = 0; i < length; i++) moduli.push_back(distr(eng)%5 + 1);
        for(vec_zero(f); vec_lessthan(f,moduli); vec_increment(f,moduli)) all.push_back(vec_string(f, ""));
        size_t size = (size_t)moduli[0];
        for(size_t i = 1; i < moduli.size(); i++) size *= moduli[i];
        if(all.size() != size) { notify("utils_vec_fail1"); failed++; }
        std::sort(all.begin(), all.end());
        bool unique = true;
        for(size_t i = 1; i < all.size(); i++) if(all[i] == all[i-1]) unique = false;
        if(!unique) { notify("utils_vec_fail2"); failed++; }
    }
    return failed;
}

bool valid_subset(std::vector<size_t>s, size_t n, size_t k) {
    if(s.size() != k) return false;
    for(size_t index : s) if(index > n) return false;
    if(!all_distinct(s)) return false;
    return true;
}

int utils_subset_tests(std::mt19937_64 &eng) {
    int failed = 0;
    for(int trial = 0; trial < 30; trial++) {
        std::uniform_int_distribution<> ndice(2,20);
        int n = ndice(eng);
        std::uniform_int_distribution<> kdice(1,n);
        int k = kdice(eng);
        std::vector<std::string> all_subset_strings;
        for(subset s(n,k); s.notdone(); s.increment()) {
            std::vector<size_t> indices = s.indices();
            if(!valid_subset(indices, n, k)) failed++;
            std::sort(indices.begin(), indices.end());
            all_subset_strings.push_back(vec_string(indices,","));

        }
        if(!all_distinct(all_subset_strings)) failed++;
        if(all_subset_strings.size() != (size_t)std::round(binomial(n,k))) failed++;
    }
    return failed;
}

int utils_tests(std::mt19937_64 &eng) {
    int failed;
    int total_failed = 0;
    failed = utils_random_subset_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All random subset tests passed." << std::endl;
    else std::cout << failed << " random subset tests failed." << std::endl;
    failed = utils_tokenize_tests();
    total_failed += failed;
    if(failed == 0) std::cout << "All tokenize tests passed." << std::endl;
    else std::cout << failed << " tokenize tests failed." << std::endl;
    failed = utils_approx_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All approx tests passed." << std::endl;
    else std::cout << failed << " approx tests failed." << std::endl;
    failed = utils_io_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All utils_io tests passed." << std::endl;
    else std::cout << failed << " utils_io tests failed." << std::endl;
    failed = utils_vec_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All utils_vec tests passed." << std::endl;
    else std::cout << failed << " utils_vec tests failed." << std::endl;
    failed = utils_math_tests();
    total_failed += failed;
    if(failed == 0) std::cout << "All utils_math tests passed." << std::endl;
    else std::cout << failed << " utils_math tests failed." << std::endl;
    failed = utils_subset_tests(eng);
    total_failed += failed;
    if(failed == 0) std::cout << "All utils_subset tests passed." << std::endl;
    else std::cout << failed << " utils_subset tests failed." << std::endl;
    return total_failed;
}


//run the tests:-------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
    std::random_device rd;
    uint64_t seed = rd();
    if(argc == 2) seed = std::stoull(argv[1]);
    std::mt19937_64 eng(seed);
    std::cout << "seed = " << seed << std::endl;

    int total_failed = bitmatrix_tests(eng);
    total_failed += bitvector_tests(eng);
    total_failed += utils_tests(eng);

    if(!total_failed) std::cout << "ALL TESTS PASSED" << std::endl;
    else std::cout << total_failed << " TESTS FAILED" << std::endl;
    return total_failed == 0 ? 0 : 1;
}
