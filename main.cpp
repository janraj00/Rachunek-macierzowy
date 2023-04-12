#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <functional>
#include <random>
#include <iomanip>

using namespace std;

#define SMALL_DIM 256

typedef vector<vector<float>> Mfloat;
typedef pair<int, int> Icoors;

Mfloat generate_matrix(size_t n, bool random = false) {
    srand(time(nullptr));
    Mfloat matrix(n);
    for (int j=0; j < n; j++) {
        matrix[j] = vector<float>(n, 0);
    }

    if (random) {
        uniform_real_distribution<float> unif(1, 9);
        default_random_engine re;
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                matrix[j][i] = unif(re);
            }
        }
    }
    return matrix;
}

bool eq_matrices(Mfloat &m, Mfloat &other) {
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

void transpose(Mfloat &src, Mfloat &dst) {
    size_t n = src.size();
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            dst[j][i] = src[i][j];
        }
    }
}

void inner_block_recursive_multiply_matrices(Mfloat &fst, size_t r_fst, size_t c_fst,
                                             Mfloat &sec, size_t r_sec, size_t c_sec,
                                             Mfloat &res, size_t r_res, size_t c_res,
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
                float sum = 0;
                for (int k=0; k<dim; k++) {
                    sum += fst[r_fst + j][c_fst + k] * sec[r_sec + k][c_sec + i];
                }
                res[r_res + j][c_res + i] += sum;
            }
        }
    }
}


void block_recursive_multiply_matrices(Mfloat &fst, Mfloat &sec, Mfloat &res) {
    size_t dim = fst.size();
    inner_block_recursive_multiply_matrices(fst, 0, 0,
                                            sec, 0, 0,
                                            res, 0, 0,
                                            dim);
}


void copy_matrix(Mfloat &fr, int r_fr, int c_fr,
                 Mfloat &to, int r_to, int c_to, size_t dim) {
    for (int j=0; j < dim; j++) {
        for (int i=0; i < dim; i++) {
            to[j + r_to][i + c_to] = fr[j + r_fr][i + c_fr];
        }
    }
}

void multiply_matrices(vector<Mfloat> &matrices, vector<Icoors> rc, Mfloat &res, Icoors rc_res, int dim) {

    inner_block_recursive_multiply_matrices(matrices[0], rc[0].first, rc[0].second,
                                                       matrices[1], rc[1].first, rc[1].second,
                                                       res, rc_res.first, rc_res.second, dim);

    for (int i=2; i < matrices.size(); i++) {
        Mfloat tmp_res = generate_matrix(dim);
        inner_block_recursive_multiply_matrices(res, rc_res.first, rc_res.second,
                                                           matrices[i], rc[i].first, rc[i].second,
                                                           tmp_res, 0, 0, dim);
        copy_matrix(tmp_res, 0, 0, res, rc_res.first, rc_res.second, dim);
    }
}

void subtruct_matrices(Mfloat &fst, int r_fst, int c_fst,
                       Mfloat &sec, int r_sec, int c_sec,
                       Mfloat &res, int r_res, int c_res,
                       int dim) {

    for (int j=0; j < dim; j++) {
        for (int i=0; i < dim; i++) {
            res[j + r_res][i + c_res] = fst[j + r_fst][i + c_fst] - sec[j + r_sec][i + c_sec];
        }
    }
}


void add_diag(Mfloat &M, int r, int c, int dim, float val = 1) {
    for (int j=0; j < dim; j++) {
        M[j + r][j + c] += val;
    }
}

void minus_matrix(Mfloat &M, int r, int c, int dim) {
    for (int j=0; j < dim; j++) {
        for (int i=0; i < dim; i++) {
            M[j + r][i + c] = -M[j + r][i + c];
        }
    }
}


Mfloat inner_inverse_matrix(Mfloat &A, int r, int c, int dim) {
    if (dim == 1) {
        return {{1 / A[r][c]}};
    }
    Mfloat B = generate_matrix(dim);

    int hdim = dim / 2;
    // A11 inverse
    Mfloat inv_A11 = inner_inverse_matrix(A, r, c, hdim);
    // S22
    Mfloat S22 = generate_matrix(hdim);
    Mfloat tmp_A_mull = generate_matrix(hdim);
    vector<Mfloat> to_multiply = {A, inv_A11, A};
    multiply_matrices(to_multiply, {{r + hdim, c}, {0, 0}, {r, c + hdim}},
                      tmp_A_mull, {0, 0}, hdim);
    /* S22 = */subtruct_matrices(A, r + hdim, c + hdim,
                      tmp_A_mull, 0, 0,
                      S22, 0, 0, hdim);
    Mfloat inv_S22 = inner_inverse_matrix(S22, 0, 0, hdim);


    // B11
    Mfloat tmp_b11_mull = generate_matrix(hdim);
    to_multiply = {A, inv_S22, A, inv_A11};
    multiply_matrices(to_multiply,
                      {{r, c + hdim}, {0, 0}, {r + hdim, c}, {0, 0}},
                      tmp_b11_mull, {0, 0}, hdim);
    add_diag(tmp_b11_mull, 0, 0, hdim);
    /* B11= */inner_block_recursive_multiply_matrices(inv_A11, 0, 0,
                                                                 tmp_b11_mull, 0, 0,
                                                                 B, r, c, hdim);


    // B12
    to_multiply = {inv_A11, A, inv_S22};
    multiply_matrices(to_multiply, {{0, 0}, {r, c + hdim}, {0, 0}},
                      B, {r, c + hdim}, hdim);
    minus_matrix(B, r, c + hdim, hdim);



    // B21
    to_multiply = {inv_S22, A, inv_A11};
    multiply_matrices(to_multiply, {{0, 0}, {r + hdim, c}, {0, 0}},
                      B, {r + hdim, c}, hdim);
    minus_matrix(B, r + hdim, c, hdim);



    // B22
    copy_matrix(inv_S22, 0, 0, B, r + hdim, c + hdim, hdim);

    return B;
}

Mfloat inverse_matrix(Mfloat &M) {
    return inner_inverse_matrix(M, 0, 0, M.size());
}

void print_mat(Mfloat &M, int width=10, int prec = 6) {
    size_t n = M.size();
    for (int j=0; j < n; j++) {
        for (int i=0; i < n; i++) {
            cout << fixed << setw(width) << setprecision(prec) << setfill(' ') << M[j][i];
        }
        cout << '\n';
    }
}


// testing

void unit_test_inverse() {
    int w = 13;
    size_t dim = 8;
    Mfloat M = generate_matrix(dim, true);
    cout << "M:\n";
    print_mat(M, w);
    cout << "\n\ninv_M:\n";
    Mfloat inv_M = inner_inverse_matrix(M, 0, 0, M.size());
    print_mat(inv_M, w);

    Mfloat M_inv_M = generate_matrix(M.size());
    inner_block_recursive_multiply_matrices(M, 0, 0,
                                            inv_M, 0, 0,
                                            M_inv_M, 0, 0, M.size());
    cout << "\n M @ inv_M:\n";
    print_mat(M_inv_M, w);
}

pair<Mfloat, Mfloat> lu_decomp_matrix(Mfloat &A, int r, int c, int dim) {
    if (dim == 1) {
        Mfloat L1 = generate_matrix(dim);
        L1[0][0] = 1;
        Mfloat U1 = generate_matrix(dim);
        U1[0][0] = A[r][c];
        return pair<Mfloat, Mfloat>(L1, U1);
    }

    Mfloat L = generate_matrix(dim);
    Mfloat U = generate_matrix(dim);

    int hdim = dim / 2;
    Mfloat zeros = generate_matrix(hdim);

    pair<Mfloat, Mfloat> LU_11 = lu_decomp_matrix(A, r, c, hdim);

    //Oblicz rekurencyjnie [L11,U11] = LU(A11)
    Mfloat L11 = LU_11.first;
    Mfloat U11 = LU_11.second;

    //Oblicz rekurencyjne U_11_-1 = inverse(U11)
    Mfloat U11_inverse = inverse_matrix(U11);

    //L_21
    Mfloat L21 = generate_matrix(hdim);
    vector<Mfloat> to_multiply = {A, U11_inverse};
    multiply_matrices(to_multiply, {{r+hdim, c}, {0, 0}},
                      L21, {r + hdim, c}, hdim);

    //L11_-1
    Mfloat L11_inverse = inverse_matrix(L11);

    //U12
    Mfloat U12 = generate_matrix(hdim);
    to_multiply = {L11_inverse, A};
    multiply_matrices(to_multiply, {{0, 0}, {r, c+hdim}},
                      U12, {r, c+hdim}, hdim);

    //S
    Mfloat S_temp = generate_matrix(hdim);
    Mfloat S = generate_matrix(hdim);
    to_multiply = {A, U11_inverse, L11_inverse, A};
    multiply_matrices(to_multiply, {{r+hdim, c}, {0, 0}, {0, 0}, {r, c+hdim}},
                      S_temp, {0, 0}, hdim);
    subtruct_matrices(A, r+hdim, c+hdim,
                      S_temp, 0, 0,
                      S, 0, 0, hdim);

    //L22, U22
    pair<Mfloat, Mfloat> LU_S = lu_decomp_matrix(S, 0, 0, hdim);
    Mfloat U22 = LU_S.second;
    Mfloat L22 = LU_S.first;
}

pair<Mfloat, Mfloat> lu_decomp(Mfloat &M) {
    return lu_decomp_matrix(M, 0, 0, M.size());
}


void test_function(const function<void(Mfloat&)>& function_being_tested, int exp, unsigned submeasure_no = 1) {
    clock_t t, total_t = 0;
    int size = pow(2, exp);
    Mfloat M = generate_matrix(size, true);
    Mfloat res = generate_matrix(size);

    cout << exp;
    for (int att=0; att < submeasure_no; att++) {
        t = clock();
        function_being_tested(M);
        t = clock() - t;
        cout << ';' << ((float)t) / CLOCKS_PER_SEC;
        total_t += t;
    }
}



void measure(const function<void(Mfloat&)>& function_being_tested, int submeasure_no, vector<int> &exps) {
    cout << "k";

    for (int e : exps) {
        test_function(function_being_tested, e, submeasure_no);
        cout << '\n';
    }

}

int main() {
    vector<int> small_exps = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    measure(inverse_matrix, 1, small_exps);
}
