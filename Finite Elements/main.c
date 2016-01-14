#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
		localB[i] = 1 / 3 * area * f;
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
			localW[i][j] = 1 / (4 * area) * scalarProduct(edges[i], edges[j], 2);
		}
	}
}

int lineCount(FILE* file) {
	int lines = 0;
	char string[20];
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
	for (i = 0; i < dim; i++) {
		vector[i] = 0;
	}
	return vector;
}

int main(int argc, char **argv) {
	int vertexSize, globalVertex2, globalVertex, meshSize;
	int tri = 0, i = 0, j = 0;
	double f = -4;
	double vertices[3][2] = { { 0 } };
	double localW[3][3] = { { 0 } };
	double localB[3] = { 0 };
	FILE* meshFile = fopen("data/p.csv", "r");
	FILE* vertexFile = fopen("data/t.csv", "r");

	vertexSize = lineCount(vertexFile);
	meshSize = lineCount(meshFile);
	double* meshPoints = readMatrixFile(meshFile, lineCount(meshFile), 2);
	double* vertexNumbers = readMatrixFile(vertexFile, lineCount(vertexFile), 3);

	double* w = zeros(meshSize * meshSize);
	double* b = zeros(meshSize);

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
				globalVertex2 = (int) *(vertexNumbers + tri * 3 + j)
						- 1;
				*(w + globalVertex * meshSize + globalVertex2) +=
						localW[i][j];
			}
		}
	}
	printMatrix(w, meshSize, meshSize, fopen("prova.txt", "w"));
	for(i = 0; i < 32; i++){
		printf("%lf\n", b[i]);
	}
	return 0;
}

