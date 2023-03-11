#include <iostream>
#include <time.h>
#include <math.h>
#include <vector>

using namespace std;

vector<vector<int>> generate_matrix(int n, bool random = true){
    srand(time(0));
    vector<vector<int>> matrix(n);
    for (int i=0; i<n; i++) {
        matrix[i] = vector<int>(n);
        for (int j=0; j<n; j++){
            matrix[i][j] = 0;
            if (random) {
                matrix[i][j] = rand() % 10;
            }
        }
    }
    return matrix;
}

void multiply_matrices(vector<vector<int>> first, vector<vector<int>> second, vector<vector<int>> &res_matrix, int n) {
    for (int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            int sum = 0;
            for(int k=0; k<n; k++){
                sum += first[i][k] * second[k][j];
            }
            res_matrix[i][j] = sum;
        }
    }
}






int main() {
    int max_exp = 10;
    for (int i=1; i<max_exp; i++){
        clock_t t;
        int size = pow(2, i);
        vector<vector<int>> first = generate_matrix(size);
        vector<vector<int>> second = generate_matrix(size);
        vector<vector<int>> res = generate_matrix(size, false);
        t = clock();
        multiply_matrices(first, second, res, size);
        t = clock() -t;
        printf ("%d %f\n", i, ((float)t)/CLOCKS_PER_SEC);
    }
}
