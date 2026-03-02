#include <iostream>
#include <vector>
#include <cmath>

int main() {

    // Get number of equations
    int n;
    std::cout << "Enter the number of equations: ";
    std::cin >> n;


    std::vector<std::vector<double>> matrix(n, std::vector<double>(n)); // matrix contrusction to solve

    // Matrix input
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            std::cout << "Enter element [" << i << "][" << j << "]: ";
            std::cin >> matrix[i][j];
        }
    }

    std::vector<double> b(n); // right side of augemented matrix

    for(int i = 0; i < n; i++) {
        std::cout << "Enter b[" << i << "]: ";
        std::cin >> b[i];
    }

    // Swap rows to give diagonal dominance
    for(int i = 0; i < n; i++) {

        int maxRow = i;
        for(int k = i+1; k < n; k++) {
            if(abs(matrix[k][i]) > abs(matrix[maxRow][i]))
                maxRow = k;
        }
        
        std::swap(matrix[i], matrix[maxRow]);
        std::swap(b[i], b[maxRow]);

    }

    // Check if the matrix can even be solved with either method
    for(int i = 0; i < n; i++) {
        if(abs(matrix[i][i]) < 1e-10) {
            std::cout << "Matrix has infinite/no solutions.";
            return 1;
        }
    }

    std::vector<double> x(n, 0.0); // x starts at zero

    // Jacobi method
    double tolerance = 1e-10; // tolerance for convergence
    int jacobi_iterations = 0;
    for(int iter = 0; iter < 100000; iter++) {
        std::vector<double> old_x = x; // save old x
        std::vector<double> new_x(n, 0.0); // new x starts at zero

        // update x values
        for(int i = 0; i < n; i++) {
            double sum = b[i];
            for(int j = 0; j < n; j++) {
                if(i != j) {
                    sum -= matrix[i][j] * x[j];
                }
            }
            new_x[i] = sum / matrix[i][i];
        }

        x = new_x;

        // check if it is within tolerance to produce a solution
        double error = 0.0;
        for(int i = 0; i < n; i++) {
            error = std::max(error, std::abs(x[i] - old_x[i]));
        }

        if(error < tolerance) {
                std::cout << "Converged in " << iter + 1 << " iterations with Jacobi method\n";
                jacobi_iterations = iter + 1;
                break;
            }

        if(error > 1e10) {
            std::cout << "Diverging! Try reordering your equations.\n";
            return 1;
        }

        if(iter == 99999) {
            std::cout << "Did not converge within 100000 iterations.\n";
            return 1;
        }

    }

    // Gauss Seidel method
    x.assign(n, 0.0); // make x zero again
    int gauss_iterations = 0;
    for(int iter = 0; iter < 100000; iter++) {
        std::vector<double> old_x = x;

        for(int i = 0; i < n; i++) {
            double sum = b[i];
            for(int j = 0; j < n; j++) {
                if(i != j) {
                    sum -= matrix[i][j] * x[j];
                }
            }
            x[i] = sum / matrix[i][i];
        }

        double error = 0.0;
        for(int i = 0; i < n; i++) {
            error = std::max(error, std::abs(x[i] - old_x[i]));
        }

        if(error < tolerance) {
                std::cout << "Converged in " << iter + 1 << " iterations with Gauss Seidel method\n";
                gauss_iterations = iter + 1;
                break;
            }

        if(error > 1e10) {
            std::cout << "Diverging! Try reordering your equations.\n";
            return 1;
        }

        if(iter == 99999) {
            std::cout << "Did not converge within 100000 iterations.\n";
            return 1;
        }


    }

    // Output solution
    for(int i = 0; i < n; i++) {
        double val = std::abs(x[i]) < tolerance ? 0.0 : x[i];
        std::cout << "x[" << i << "] = " << val << "\n";
    }

    // find how much faster Gauss Seidel is than Jacobi
    int diff = jacobi_iterations - gauss_iterations;
    double pct = (double)diff / jacobi_iterations * 100;
    std::cout << "Gauss Seidel converged in " << pct << "% fewer iterations\n";

    return 0;
}