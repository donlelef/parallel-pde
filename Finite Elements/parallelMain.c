#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define TOL 10e-3

void printMatrix(double *matrix, int xDim, int yDim, FILE* fout) {
	int i = 0, j = 0;
	for (i = 0; i < xDim; i++) {
		for (j = 0; j < yDim; j++) {
			fprintf(fout, "%lf ", *(matrix + i * yDim + j));
		}
		fprintf(fout, "\n");
	}
	fclose(fout);
}

void printSparseMatrix(double *matrix, int xDim, int yDim, FILE* fout) {
	int i = 0, j = 0;
	for (i = 0; i < xDim; i++) {
		for (j = 0; j < yDim; j++) {
			double element = *(matrix + i * yDim + j);
			if (element != 0.0) {
				fprintf(fout, "(%d %d) %lf\n", i, j, *(matrix + i * yDim + j));
			}
		}
	}
	fclose(fout);
}

double triangleArea(double vertices[][2]) {
	double area = 0;
	area = 0.5
			* ((vertices[1][0] - vertices[0][0])
					* (vertices[2][1] - vertices[0][1])
					- (vertices[2][0] - vertices[0][0])
							* (vertices[1][1] - vertices[0][1]));
	area = fabs(area);
	return area;
}

void localVector(double vertices[][2], double f, double localB[]) {
	int i;
	double area = triangleArea(vertices);
	for (i = 0; i < 3; i++) {
		localB[i] = (1.0 / 3) * area * f;
	}
}

double scalarProduct(double x[], double y[], int length) {
	int i = 0;
	double res = 0;
	for (i = 0; i < length; i++) {
		res += x[i] * y[i];
	}
	return res;
}

void localStiffnessMatrix(double vertices[][2], double localW[][3]) {
	int i = 0, j = 0;
	double edges[3][2];

	edges[0][0] = vertices[2][0] - vertices[1][0];
	edges[0][1] = vertices[2][1] - vertices[1][1];
	edges[1][0] = vertices[0][0] - vertices[2][0];
	edges[1][1] = vertices[0][1] - vertices[2][1];
	edges[2][0] = vertices[1][0] - vertices[0][0];
	edges[2][1] = vertices[1][1] - vertices[0][1];

	double area = triangleArea(vertices);

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			localW[i][j] = 1 / (4 * area)
					* scalarProduct(edges[i], edges[j], 2);
		}
	}
}

int lineCount(FILE* file) {
	int lines = 0;
	char string[100];
	rewind(file);

	while (fgets(string, sizeof(string), file) != NULL) {
		lines++;
	}
	return lines;
}

double* readMatrixFile(FILE* file, int rows, int columns) {
	int j = 0, i = 0;
	char string[100];
	rewind(file);

	double* matrix = (double*) malloc(rows * columns * sizeof(double));
	while (fgets(string, sizeof(string), file) != NULL) {
		for (j = 0; j < columns; j++) {
			*(matrix + i * columns + j) = atof(
					strtok(j == 0 ? string : NULL, ","));
		}
		i++;
	}
	return matrix;
}

double* zeros(int dim) {
	int i = 0;
	double* vector = (double*) malloc(sizeof(double) * dim);
#pragma omp parallel for private(i) schedule(static)
	for (i = 0; i < dim; i++) {
		vector[i] = 0;
	}
	return vector;
}

void gaussianElimination(int n, double* matrix, double b[]){
	int i, j, k;
	double aux;

	/* Gaussian elimination */
	for (i = 0; i < (n - 1); i++) {
#pragma omp parallel for private(aux, j, k) schedule(dynamic)
		for (j = (i + 1); j < n; j++) {
			aux = *(matrix + n * j + i) / *(matrix + n * i + i);
			for (k = i; k < n; k++) {
				*(matrix + n * j + k) -= (aux * (*(matrix + n * i + k)));
			}
			b[j] -= (aux * b[i]);
		}
	}
}

double* backSubstitution(int n, double* matrix, double b[]) {

	int i, j;
	double temp;
	double* x = zeros(n);

	/* Back substitution */
	x[n - 1] = b[n - 1] / *(matrix + n * (n - 1) + n - 1);
	for (i = (n - 2); i >= 0; i--) {
		temp = b[i];
		for (j = (i + 1); j < n; j++) {
			temp -= (*(matrix + n * i + j) * x[j]);
		}
		x[i] = temp / *(matrix + n * i + i);
	}

	return x;
}

void assignDirichletCondition(int meshSize, double g, double* meshPoints,
		double* w, double* b) {
	int i, j;

#pragma omp parallel for private(i, j) schedule(dynamic)
	for (i = 0; i < meshSize; i++) {
		double x = *(meshPoints + 2 * i);
		double y = *(meshPoints + 2 * i + 1);
		if (fabs(x * x + y * y - 1) < TOL) {
			for (j = 0; j < meshSize; j++) {
				*(w + meshSize * i + j) = (i == j) ? 1 : 0;
			}
			b[i] = g;
		}
	}
}

double* solvePDE(char fileP[], char fileT[], FILE* result) {

	int vertexSize, globalVertex2, globalVertex, meshSize;
	int tri = 0, i = 0, j = 0;
	double f = -4, g = 1;
	double vertices[3][2] = { { 0 } };
	double localW[3][3] = { { 0 } };
	double localB[3] = { 1 };

	FILE* meshFile = fopen(fileP, "r");
	FILE* vertexFile = fopen(fileT, "r");

	vertexSize = lineCount(vertexFile);
	meshSize = lineCount(meshFile);

	/* Read files: START */
	double startIO = omp_get_wtime();
	double* meshPoints = readMatrixFile(meshFile, lineCount(meshFile), 2);
	double* vertexNumbers = readMatrixFile(vertexFile, lineCount(vertexFile),
			3);
	double ioTime = omp_get_wtime() - startIO;
	/* Read files: END */

	/* Prepare matrices: START */
	double startZeros = omp_get_wtime();
	double* w = zeros(meshSize * meshSize);
	double* b = zeros(meshSize);
	double zerosTime = omp_get_wtime() - startZeros;
	/* Prepare matrices: START */

	/* Assembling: START */
	double startComputation = omp_get_wtime();
	for (tri = 0; tri < vertexSize; tri++) {
		for (i = 0; i < 3; i++) {
			globalVertex = (int) *(vertexNumbers + tri * 3 + i) - 1;
			vertices[i][0] = *(meshPoints + globalVertex * 2 + 0);
			vertices[i][1] = *(meshPoints + globalVertex * 2 + 1);
		}

		localStiffnessMatrix(vertices, localW);
		localVector(vertices, f, localB);

		for (i = 0; i < 3; i++) {
			globalVertex = (int) *(vertexNumbers + tri * 3 + i) - 1;
			b[globalVertex] += localB[i];
			for (j = 0; j < 3; j++) {
				globalVertex2 = (int) *(vertexNumbers + tri * 3 + j) - 1;
				*(w + globalVertex * meshSize + globalVertex2) += localW[i][j];
			}
		}
	}
	double assemblingTime = omp_get_wtime() - startComputation;
	/* Assembling: END */

	/* Dirichlet condition: START */
	double startDirichlet = omp_get_wtime();
	assignDirichletCondition(meshSize, g, meshPoints, w, b);
	double dirichletTime = omp_get_wtime() - startDirichlet;
	/* Dirichlet condition: END */

	/* Gaussian Elimination: START */
	double startGauss = omp_get_wtime();
	gaussianElimination(meshSize, w, b);
	double gaussTime = omp_get_wtime() - startGauss;
	/* Gaussian Elimination: END */

	/* Backward Substitution: START */
	double startBack = omp_get_wtime();
	double* u = backSubstitution(meshSize, w, b);
	double backTime = omp_get_wtime() - startBack;
	/* Backward Substitution: END */

	double computationTime = omp_get_wtime() - startComputation;
	/* Computation: END */

	printf("Elapsed time for computation: %lf seconds.\n", computationTime);
	printf("Elapsed time for IO: %lf seconds.\n", ioTime);
	printf("Elapsed time for zeros: %lf seconds.\n", zerosTime);
	printf("Elapsed time for assembling: %lf seconds.\n", assemblingTime);
	printf("Elapsed time for dirichlet: %lf seconds.\n", dirichletTime);
	printf("Elapsed time for gaussian elimination: %lf seconds.\n", gaussTime);
	printf("Elapsed time for back substitution: %lf seconds.\n", backTime);

	fprintf(result, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n", computationTime, ioTime, zerosTime, assemblingTime, dirichletTime, gaussTime, backTime);

	return u;
}

int main(int argc, char **argv) {

	int i;
	char dirP[50] = "data/p";
	char dirT[50] = "data/t";
	char extension[50] = ".csv";
	char fileP[100];
	char fileT[100];
	int files = atoi(argv[1]);

	omp_set_num_threads(4);

	FILE* result = fopen("resParallel.csv", "w");
	fprintf(result, "computationTime, ioTime, zerosTime, assemblingTime, dirichletTime, gaussTime, backTime\n");

	for (i = 0; i < files; i++) {
		sprintf(fileP, "%s%d%s", dirP, i, extension);
		sprintf(fileT, "%s%d%s", dirT, i, extension);
		printf("Iteration %d\n", i);
		solvePDE(fileP, fileT, result);
	}

	fclose(result);

	return 0;
}
