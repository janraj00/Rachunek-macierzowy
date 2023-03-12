#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <functional>

using namespace std;

#define SMALL_DIM 4

typedef vector<vector<int>> MInt;

MInt generate_matrix(size_t n, bool random = true){
    srand(time(nullptr));
    MInt matrix(n);
    for (int j=0; j < n; j++) {
        matrix[j] = vector<int>(n);
    }

    if (random) {
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                matrix[j][i] = rand() % 10;
            }
        }
    }
    return matrix;
}

bool eq_matrices(MInt &m, MInt &other) {
    size_t n = m.size();
    for (int j=0; j<n; j++) {
        for (int i=0; i<n; i++) {
            if (m[j][i] != other[j][i]) {
                return false;
            }
        }
    }
    return true;
}

void transpose(MInt &src, MInt &dst) {
    size_t n = src.size();
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            dst[j][i] = src[i][j];
        }
    }
}

void classic_multiply_matrices(MInt &first, MInt &second, MInt &res_matrix) {
    size_t n = first.size();
    for (int j=0; j<n; j++){
        for(int i=0; i<n; i++){
            int sum = 0;
            for(int k=0; k<n; k++){
                sum += first[j][k] * second[k][i];
            }
            res_matrix[j][i] = sum;
        }
    }
}

void classic_transposed_multiply_matrices(MInt &first, MInt &second, MInt &res_matrix) {
    size_t n = first.size();
    MInt transposed_second(n);
    for (int j=0; j < n; j++) {
        transposed_second[j] = vector<int>(n);
    }
    transpose(second, transposed_second);

    for (int j=0; j<n; j++){
        for(int i=0; i<n; i++){
            int sum = 0;
            for(int k=0; k<n; k++){
                sum += first[j][k] * transposed_second[i][k];
            }
            res_matrix[j][i] = sum;
        }
    }
}


void inner_block_recursive_multiply_matrices(MInt &fst, size_t r_fst, size_t c_fst,
                                             MInt &sec, size_t r_sec, size_t c_sec,
                                             MInt &res, size_t r_res, size_t c_res,
                                             size_t dim) {
    if (dim > SMALL_DIM) {
        size_t rd = dim / 2;
        // C 11 = A11*B11 + A12*B21
        inner_block_recursive_multiply_matrices(fst, r_fst, c_fst,
                                                sec, r_sec, c_sec,
                                                res, r_res, c_res, rd);
        inner_block_recursive_multiply_matrices(fst, r_fst, c_fst + rd,
                                                sec, r_sec + rd, c_sec,
                                                res, r_res, c_res, rd);
        // C 12 = A11*B12 + A12*B22
        inner_block_recursive_multiply_matrices(fst, r_fst, c_fst,
                                                sec, r_sec, c_sec + rd,
                                                res, r_res, c_res + rd, rd);
        inner_block_recursive_multiply_matrices(fst, r_fst, c_fst + rd,
                                                sec, r_sec + rd, c_sec + rd,
                                                res, r_res, c_res + rd, rd);
        // C 21 = A21*B11 + A22*B21
        inner_block_recursive_multiply_matrices(fst, r_fst + rd, c_fst,
                                                sec, r_sec, c_sec,
                                                res, r_res + rd, c_res, rd);
        inner_block_recursive_multiply_matrices(fst, r_fst + rd, c_fst + rd,
                                                sec, r_sec + rd, c_sec,
                                                res, r_res + rd, c_res, rd);

        // C 22 = A21*B12 + A22*B22
        inner_block_recursive_multiply_matrices(fst, r_fst + rd, c_fst,
                                                sec, r_sec, c_sec + rd,
                                                res, r_res + rd, c_res + rd, rd);
        inner_block_recursive_multiply_matrices(fst, r_fst + rd, c_fst + rd,
                                                sec, r_sec + rd, c_sec + rd,
                                                res, r_res + rd, c_res + rd, rd);


    } else {
        for (int j=0; j<dim; j++) {
            for (int i=0; i<dim; i++) {
                int sum = 0;
                for (int k=0; k<dim; k++) {
                    sum += fst[r_fst + j][c_fst + k] * sec[r_sec + k][c_sec + i];
                }
                res[r_res + j][c_res + i] += sum;
            }
        }
    }
}

void inner_transposed_block_recursive_multiply_matrices(MInt &fst, size_t r_fst, size_t c_fst,
                                                        MInt &sec, size_t r_sec, size_t c_sec,
                                                        MInt &res, size_t r_res, size_t c_res,
                                                        size_t dim) {
    if (dim > SMALL_DIM) {
        size_t rd = dim / 2;
        // C 11 = A11*B11 + A12*B21
        inner_transposed_block_recursive_multiply_matrices(fst, r_fst, c_fst,
                                                sec, r_sec, c_sec,
                                                res, r_res, c_res, rd);
        inner_transposed_block_recursive_multiply_matrices(fst, r_fst, c_fst + rd,
                                                sec, r_sec + rd, c_sec,
                                                res, r_res, c_res, rd);
        // C 12 = A11*B12 + A12*B22
        inner_transposed_block_recursive_multiply_matrices(fst, r_fst, c_fst,
                                                sec, r_sec, c_sec + rd,
                                                res, r_res, c_res + rd, rd);
        inner_transposed_block_recursive_multiply_matrices(fst, r_fst, c_fst + rd,
                                                sec, r_sec + rd, c_sec + rd,
                                                res, r_res, c_res + rd, rd);
        // C 21 = A21*B11 + A22*B21
        inner_transposed_block_recursive_multiply_matrices(fst, r_fst + rd, c_fst,
                                                sec, r_sec, c_sec,
                                                res, r_res + rd, c_res, rd);
        inner_transposed_block_recursive_multiply_matrices(fst, r_fst + rd, c_fst + rd,
                                                sec, r_sec + rd, c_sec,
                                                res, r_res + rd, c_res, rd);

        // C 22 = A21*B12 + A22*B22
        inner_transposed_block_recursive_multiply_matrices(fst, r_fst + rd, c_fst,
                                                sec, r_sec, c_sec + rd,
                                                res, r_res + rd, c_res + rd, rd);
        inner_transposed_block_recursive_multiply_matrices(fst, r_fst + rd, c_fst + rd,
                                                sec, r_sec + rd, c_sec + rd,
                                                res, r_res + rd, c_res + rd, rd);


    } else {
        for (int j=0; j<dim; j++) {
            for (int i=0; i<dim; i++) {
                int sum = 0;
                for (int k=0; k<dim; k++) {
                    sum += fst[r_fst + j][c_fst + k] * sec[c_sec + i][r_sec + k];
                }
                res[r_res + j][c_res + i] += sum;
            }
        }
    }
}

void block_recursive_multiply_matrices(MInt &fst, MInt &sec, MInt &res) {
    size_t dim = fst.size();
    inner_block_recursive_multiply_matrices(fst, 0, 0,
                                            sec, 0, 0,
                                            res, 0, 0,
                                            dim);
}

void block_transposed_recursive_multiply_matrices(MInt &fst, MInt &sec, MInt &res) {
    size_t dim = fst.size();
    MInt transposed_second(dim);
    for (int j = 0; j < dim; j++) {
        transposed_second[j] = vector<int>(dim);
    }
    transpose(sec, transposed_second);
    inner_transposed_block_recursive_multiply_matrices(fst, 0, 0,
                                                       transposed_second, 0, 0,
                                                       res, 0, 0,
                                                       dim);
}





void test_algo_correctness() {
    size_t size = 128;
    MInt fst = generate_matrix(size);
    MInt sec = generate_matrix(size);
    MInt res1 = generate_matrix(size, false);
    MInt res2 = generate_matrix(size, false);
    MInt res3 = generate_matrix(size, false);
    MInt res4 = generate_matrix(size, false);

    //
    classic_multiply_matrices(fst, sec, res1);
    classic_transposed_multiply_matrices(fst, sec, res2);
    block_recursive_multiply_matrices(fst, sec, res3);
    block_transposed_recursive_multiply_matrices(fst, sec, res4);

    cout << boolalpha << eq_matrices(res1, res2) << endl; // true if correct
    cout << boolalpha << eq_matrices(res1, res3) << endl;
    cout << boolalpha << eq_matrices(res1, res4) << endl;
}

void test_function(const function<void(MInt&, MInt&, MInt&)>& function_being_tested, int exp, string &name, unsigned attempts = 1) {
    printf("<|%s|-[ex:*att=%d*:avg]>:\n", name.c_str(), attempts);
        clock_t t, total_t = 0;
        int size = pow(2, exp);
        MInt first = generate_matrix(size);
        MInt second = generate_matrix(size);
        MInt res = generate_matrix(size, false);

        cout << exp << " ";
        for (int att=0; att<attempts; att++) {
            t = clock();
            function_being_tested(first, second, res);
            t = clock() - t;
            cout << ((float)t) / CLOCKS_PER_SEC << " ";
            total_t += t;
        }
        cout << ((float)total_t / attempts) / CLOCKS_PER_SEC << endl;
}

void measure(string &name, int attempts) {

}

int main() {
    int max_exp = 10;
    string name = "classic";
    test_function(classic_transposed_multiply_matrices, 10, name, 2);
}
