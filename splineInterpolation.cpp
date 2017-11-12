#include <iostream>
#include <cmath>

double f(double x, int k, int m) {
	return pow(tan(M_PI * pow(x, m) / 4), k);
}

// метод прогонки
double *sweep(double **matrix, int n) {
	double* a = new double[n];
    double* b = new double[n];
    double* x = new double[n + 1];

    a[0] = b[0] = x[0] = x[n] = 0.0;

    for (int i = 0; i < n - 1; i++) {
        a[i + 1] = -matrix[2][i] / (matrix[1][i] + matrix[0][i] * a[i]);
        b[i + 1] = (matrix[3][i] - matrix[0][i] * b[i]) / (matrix[1][i] + matrix[0][i] * b[i]);
    }
    
    for (int i = n; i > 1; i--)
        x[i - 1] = a[i - 1] * x[i] + b[i - 1];

    delete[] a;
    delete[] b;

    return x;
}

double interpolate(double x0, int n, double *x, double *a, double *b, double *c, double *d) {
	int i = 0;
    while (x[i + 1] < x0 && i < n - 1)
    	i++;

    double dx = (x0 - x[i]);

	return a[i] + b[i] * dx + c[i] * pow(dx, 2) + d[i] * pow(dx, 3);
}


void freeMemory(double **a, double **b, double **c, double **d, double **x, double **y, double ***matrix) {
	delete[] *a;
	delete[] *b;
	delete[] *c;
	delete[] *d;
	delete[] *x;
	delete[] *y;

	for (int i = 0; i < 4; i++)
		delete[] (*matrix)[i];

	delete[] *matrix;
}

int main() {
	int n, k, m;
	double left, right;

	printf("Enter a: ");
	scanf("%lf", &left);
	printf("Enter b: ");
	scanf("%lf", &right);

	printf("Enter n: ");
	scanf("%d", &n);
	printf("Enter k: ");
	scanf("%d", &k);
	printf("Enter m: ");
	scanf("%d", &m);

	double *x = new double[n + 1];
	double *y = new double[n + 1];

	double h = (right - left) / n;

	for (int i = 0; i <= n; i++) {
		x[i] = left + i * h;
		y[i] = f(x[i], k, m);
	}
  
    double** matrix = new double*[4];
    for (int i = 0; i < 4; i++) {
    	matrix[i] = new double[n - 1];

    	for (int j = 0; j < n - 1; j++)
    		matrix[i][j] = 0;
    }

    for (int i = 0; i < n - 1; i++) {
        matrix[0][i] = matrix[2][i] = h;
        matrix[1][i] = 4 * h;
        matrix[3][i] = 3 * (y[i + 2] - 2 * y[i + 1] + y[i]) / h;
    }
    
    double *c = sweep(matrix, n);
    double *a = new double[n + 1];
    double *b = new double[n + 1];
	double *d = new double[n + 1];

	for (int i = 0; i < n; i++) {
		a[i] = y[i];
		b[i] = (y[i + 1] - y[i]) / h - h * (c[i + 1] + c[i]) / 3.0;
        d[i] = ((c[i + 1] - c[i]) / 3 / h);
	}

    double x0;
    printf("Enter x0: ");
    scanf("%lf", &x0);

	printf("y(%lf) = %.10lf\n", x0, interpolate(x0, n, x, a, b, c, d));

	printf("\nEnter h for table: ");
	scanf("%lf", &h);

	printf("+=======+================+===================+=====================+\n");
	printf("|   x   |    Real f(x)   | Interpolated f(x) |    |f(x) - s(x)|    |\n");	
	printf("+=======+================+===================+=====================+\n");
	for (double x0 = left; x0 <= right; x0 += h) {
		double f0 = f(x0, k, m);
		double s0 = interpolate(x0, n, x, a, b, c, d);

		printf("| %5.2lf | %14.12lf | %17.12lf | %19.17lf |\n", x0, f0, s0, fabs(f0 - s0));
	}
	printf("+=======+================+===================+=====================+\n");

	freeMemory(&a, &b, &c, &d, &x, &y, &matrix);
}