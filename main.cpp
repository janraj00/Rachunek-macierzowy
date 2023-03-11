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

vector<long long> multiply_matrices(vector<vector<int>> first, vector<vector<int>> second, vector<vector<int>> &res_matrix, int n) {
    long long adds =0;
    long long muls = 0;
    for (int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            int sum = 0;
            for(int k=0; k<n; k++){
                sum += first[i][k] * second[k][j];
                adds++;
                muls++;
            }
            res_matrix[i][j] = sum;
        }
    }
    return vector<long long>{adds, muls};
}

void run_traditional_pipeline(int max_exp, vector<vector<long long>> &operations, vector<float> &times){
    operations.push_back(vector<long long>{0, 0});
    times.push_back(0.0);
    for (int i=1; i<=max_exp; i++){
        clock_t t;
        int size = pow(2, i);
        vector<vector<int>> first = generate_matrix(size);
        vector<vector<int>> second = generate_matrix(size);
        vector<vector<int>> res = generate_matrix(size, false);
        t = clock();
        vector<long long> ops = multiply_matrices(first, second, res, size);
        operations.push_back(ops);
        t = clock() -t;
        float time = ((float) t)/CLOCKS_PER_SEC;
        times.push_back(time);
    }
}

void print_results(int max_exp, vector<vector<long long>> &operations, vector<float> &times) {
    for (int i = 1; i <= max_exp; i++) {
        cout << i << " " << times[i] << " " << operations[i][0] << " " << operations[i][1] << endl;
    }
}


int main() {
    int max_exp = 10;
    vector<vector<long long>> operations;
    vector<float> times;
    run_traditional_pipeline(max_exp, operations, times);
    print_results(max_exp, operations, times);
}
