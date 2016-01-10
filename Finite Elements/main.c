#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double determinant(double matrix[][2]) {
	double det = ((matrix[1][0] - matrix[0][0]) * (matrix[2][1] - matrix[0][1])
			- (matrix[2][0] - matrix[0][0]) * (matrix[1][1] - matrix[0][1]));
	return det;
}

double triangleArea(double vertices[][2]) {
	double area = 0;
	area = 0.5 * determinant(vertices);
	return area;
}

void localVector(double vertices[][2], double f, double localB[]) {
	int i;
	double area = triangleArea(vertices);
	for (i = 0; i < 3; i++) {
		localB[i] = 1 / 3 * area * f;
	}
}

double scalarProduct(double x[], double y[]) {
	int i = 0;
	double res = 0;
	for (i = 0; i < sizeof(x) / sizeof(double); i++) {
		res += x[1] * y[i];
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
			localW[i][j] = 1 / (4 * area) * scalarProduct(edges[i], edges[j]);
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

void printMatrix(double *matrix, int xDim, int yDim, char* fileName) {
	int i = 0, j = 0;
	FILE* fout = fopen(fileName, "w");
	for (i = 0; i < xDim; i++) {
		for (j = 0; j < yDim; j++) {
			fprintf(fout, "%lf ", *(matrix + i * yDim + j));
		}
		fprintf(fout, "\n");
	}
	fclose(fout);
}

int main(int argc, char **argv) {

	int f = -4, vertexSize, globalVertex2, globalVertex;
	int tri = 0, localVert = 0, localVertex2 = 0;
	double vertices[3][2] = { { 0 } };
	double localW[3][3] = { { 0 } };
	double localB[3] = { 0 };
	FILE* meshFile = fopen("data/p.csv", "r");
	FILE* vertexFile = fopen("data/t.csv", "r");

	double* meshPoints = readMatrixFile(meshFile, lineCount(meshFile), 2);
	vertexSize = lineCount(vertexFile);
	double* vertexNumbers = readMatrixFile(vertexFile, vertexSize, 3);

	double* w = zeros(vertexSize * vertexSize);
	double* b = zeros(vertexSize);

	for (tri = 0; tri < vertexSize; tri++) {
		for (localVert = 0; localVert < 3; localVert++) {
			globalVertex = (int) *(vertexNumbers + tri * 3 + localVert);
			vertices[localVert][0] = *(meshPoints + globalVertex * 2 + 0);
			vertices[localVert][1] = *(meshPoints + globalVertex * 2 + 1);
		}

		localStiffnessMatrix(vertices, localW);
		localVector(vertices, f, localB);

		for (localVert = 0; localVert < 3; localVert++) {
			globalVertex = *(vertexNumbers + tri * 3 + localVert);
			b[globalVertex] += localB[localVert];
			for (localVertex2 = 0; localVertex2 < 3; localVertex2++) {
				globalVertex2 = (int) *(vertexNumbers + tri * 3 + localVertex2);
				*(w + globalVertex * vertexSize + globalVertex2) +=
						localW[localVert][localVertex2];
			}
		}
	}
	printMatrix(w, vertexSize, vertexSize, "prova.txt");

	return 0;

}

