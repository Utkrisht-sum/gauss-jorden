#include <stdio.h>
#include <math.h>

#define EPSILON 1e-9

// Function to perform Gaussian elimination
void gaussElimination(int rows, int cols, double matrix[][cols], double solution[]) {
    int i, j, k, maxRow;
    double ratio;
    int rank = 0;  // To track the rank of the matrix

    // Gaussian Elimination process
    for (i = 0; i < rows; i++) {
        // Find the pivot row and swap with the current row
        maxRow = i;
        for (k = i + 1; k < rows; k++) {
            if (fabs(matrix[k][i]) > fabs(matrix[maxRow][i])) {
                maxRow = k;
            }
        }

        // Swap the rows if needed
        for (j = 0; j < cols; j++) {
            double temp = matrix[i][j];
            matrix[i][j] = matrix[maxRow][j];
            matrix[maxRow][j] = temp;
        }

        // Check if the pivot element is too small (i.e., almost zero)
        if (fabs(matrix[i][i]) < EPSILON) {
            continue;  // Skip this row for elimination
        }

        rank++;  // This row has contributed to the rank of the matrix

        // Make all rows below this one 0 in the current column
        for (k = i + 1; k < rows; k++) {
            ratio = -matrix[k][i] / matrix[i][i];
            for (j = i; j < cols; j++) {
                matrix[k][j] += ratio * matrix[i][j];
            }
        }
    }

    // Check if there are infinite solutions or no solution
    for (i = rank; i < rows; i++) {
        if (fabs(matrix[i][cols - 1]) > EPSILON) {
            printf("No solution exists.\n");
            return;
        }
    }

    if (rank < cols - 1) {
        printf("Warning: Infinite solutions exist.\n");
        return;
    }

    // Perform back-substitution to find the solutions
    for (i = cols - 2; i >= 0; i--) {
        solution[i] = matrix[i][cols - 1] / matrix[i][i];
        for (k = i - 1; k >= 0; k--) {
            matrix[k][cols - 1] -= matrix[k][i] * solution[i];
        }
    }

    // Display the solutions
    printf("Solutions:\n");
    for (i = 0; i < cols - 1; i++) {
        printf("x%d = %.6f\n", i + 1, solution[i]);
    }
}

int main() {
    int rows, vars, i, j;

    // Input number of variables and equations
    printf("Enter the number of equations: ");
    scanf("%d", &rows);

    printf("Enter the number of variables: ");
    scanf("%d", &vars);

    // Augmented matrix has rows and (vars + 1) columns
    double matrix[rows][vars + 1], solution[vars];

    // Input the augmented matrix
    printf("Enter the augmented matrix (coefficients and constants):\n");
    for (i = 0; i < rows; i++) {
        for (j = 0; j <= vars; j++) {
            scanf("%lf", &matrix[i][j]);
        }
    }

    // Call the Gaussian Elimination function
    gaussElimination(rows, vars + 1, matrix, solution);

    return 0;
}
