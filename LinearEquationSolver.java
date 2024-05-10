import java.util.Arrays;
//Muhammad Najib Pradana
//21120122120006
public class LinearEquationSolver {
    public static double[] solveByInverse(double[][] A, double[] b) {

        int n = A.length;
        double[][] inverse = inverse(A, n);
        double[] x = new double[n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                x[i] += inverse[i][j] * b[j];
            }
        }

        return x;
    }

    private static double[][] inverse(double[][] A, int n) {
        double[][] inverse = new double[n][n];

        // Membuat matriks identitas
        for (int i = 0; i < n; i++) {
            Arrays.fill(inverse[i], 0);
            inverse[i][i] = 1;
        }

        // Menghitung matriks balikan
        for (int i = 0; i < n; i++) {
            double temp = 1.0 / A[i][i];

            for (int j = 0; j < n; j++) {
                A[i][j] *= temp;
                inverse[i][j] *= temp;
            }

            for (int x = 0; x < n; x++) {
                if (x != i) {
                    double lam = A[x][i];

                    for (int j = 0; j < n; j++) {
                        A[x][j] -= lam * A[i][j];
                        inverse[x][j] -= lam * inverse[i][j];
                    }
                }
            }
        }

        return inverse;
    }

    public static double[] solveBySolveLU(double[][] A, double[] b) {

        int n = A.length;
        double[][] L = new double[n][n];
        double[][] U = new double[n][n];

        // Dekomposisi LU
        for (int i = 0; i < n; i++) {
            // Matriks L
            for (int j = 0; j <= i; j++) {
                double sum = 0;
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * U[k][j];
                }
                L[i][j] = A[i][j] - sum;
            }

            // Matriks U
            for (int j = i; j < n; j++) {
                double sum = 0;
                for (int k = 0; k < i; k++) {
                    sum += L[i][k] * U[k][j];
                }
                U[i][j] = (A[i][j] - sum) / L[i][i];
            }
        }

        // Menyelesaikan Ly = b
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += L[i][j] * y[j];
            }
            y[i] = b[i] - sum;
        }

        // Menyelesaikan Ux = y
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < n; j++) {
                sum += U[i][j] * x[j];
            }
            x[i] = (y[i] - sum) / U[i][i];
        }

        return x;
    }

    public static double[] solveByCrout(double[][] A, double[] b) {

        int n = A.length;
        double[][] L = new double[n][n];
        double[][] U = new double[n][n];

        // Dekomposisi Crout
        for (int j = 0; j < n; j++) {
            L[j][j] = 1;

            for (int i = j; i < n; i++) {
                double sum = 0;
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * U[k][j];
                }
                U[j][i] = A[j][i] - sum;
            }

            for (int i = j + 1; i < n; i++) {
                double sum = 0;
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * U[k][j];
                }
                L[i][j] = (A[i][j] - sum) / U[j][j];
            }
        }

        // Menyelesaikan Ly = b
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += L[i][j] * y[j];
            }
            y[i] = b[i] - sum;
        }

        // Menyelesaikan Ux = y
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < n; j++) {
                sum += U[i][j] * x[j];
            }
            x[i] = (y[i] - sum) / U[i][i];
        }

        return x;
    }
    public static void main(String[] args) {
        // Kasus uji 1
        double[][] A1 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
        double[] b1 = {6, 15, 24};

        // Kasus uji 2
        double[][] A2 = {{1, 1, 1}, {2, 3, 1}, {1, 2, 4}};
        double[] b2 = {3, 5, 7};

        // Kasus uji 3
        double[][] A3 = {{3, -1, 2}, {1, 2, 1}, {2, 1, 3}};
        double[] b3 = {5, 7, 11}; // Inisialisasi vektor b3

        double[] x1_inverse = solveByInverse(A1, b1);
        double[] x1_LU = solveBySolveLU(A1, b1);
        double[] x1_Crout = solveByCrout(A1, b1);

        // Penyelesaian untuk kasus uji 2
        double[] x2_inverse = solveByInverse(A2, b2);
        double[] x2_LU = solveBySolveLU(A2, b2);
        double[] x2_Crout = solveByCrout(A2, b2);

        // Penyelesaian untuk kasus uji 3
        double[] x3_inverse = solveByInverse(A3, b3);
        double[] x3_LU = solveBySolveLU(A3, b3);
        double[] x3_Crout = solveByCrout(A3, b3);

        // Menampilkan hasil penyelesaian
        System.out.println("Hasil Penyelesaian:");
        System.out.println("Kasus Uji 1:");
        System.out.println("Metode Matriks Balikan: " + Arrays.toString(x1_inverse));
        System.out.println("Metode Dekomposisi LU Gauss: " + Arrays.toString(x1_LU));
        System.out.println("Metode Dekomposisi Crout: " + Arrays.toString(x1_Crout));

        System.out.println("Kasus Uji 2:");
        System.out.println("Metode Matriks Balikan: " + Arrays.toString(x2_inverse));
        System.out.println("Metode Dekomposisi LU Gauss: " + Arrays.toString(x2_LU));
        System.out.println("Metode Dekomposisi Crout: " + Arrays.toString(x2_Crout));

        System.out.println("Kasus Uji 3:");
        System.out.println("Metode Matriks Balikan: " + Arrays.toString(x3_inverse));
        System.out.println("Metode Dekomposisi LU Gauss: " + Arrays.toString(x3_LU));
        System.out.println("Metode Dekomposisi Crout: " + Arrays.toString(x3_Crout));
    }

    }
