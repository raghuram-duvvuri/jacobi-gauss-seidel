# jacobi-gauss-seidel

A C++ program that solves systems of linear equations using two iterative methods: **Jacobi** and **Gauss-Seidel**.

## How to Compile and Run

```bash
g++ solver.cpp -o solver
./solver
```

## Example

For the system:
```
4x + y + z = 6
x + 4y + z = 6
x + y + 4z = 6
```

Input:
```
Enter the number of equations: 3
Enter element [0][0]: 4
Enter element [0][1]: 1
Enter element [0][2]: 1
Enter element [1][0]: 1
Enter element [1][1]: 4
Enter element [1][2]: 1
Enter element [2][0]: 1
Enter element [2][1]: 1
Enter element [2][2]: 4
Enter b[0]: 6
Enter b[1]: 6
Enter b[2]: 6
```

Output:
```
Converged in 21 iterations with Jacobi method
Converged in 8 iterations with Gauss Seidel method
x[0] = 1
x[1] = 1
x[2] = 1
Gauss Seidel converged in 61.9% fewer iterations
```

## How It Works

Both methods solve **Ax = b** by rearranging each equation to isolate one variable and iterating until the solution stops changing.

- **Jacobi** — computes all new values using only the previous iteration's values, then updates all at once
- **Gauss-Seidel** — updates each variable immediately, so newer values are used within the same iteration, leading to faster convergence

## Notes

- The program automatically reorders rows to improve diagonal dominance before solving
- Convergence is not guaranteed for all matrices — both methods work best when the matrix is diagonally dominant
- Tolerance is set to `1e-10`
