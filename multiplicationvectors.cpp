#include <stdio.h>
#include <iostream>
using namespace std;

void freeSpaceMatrix(double** res, int row) {
	if (res) {
		for (int i = 0; i < row; i++) { if (res[i]) { delete[] res[i]; } }
		delete[] res;
	}
}

double** multiMatrixes(double a[], double b[], const int row, int const column) {
	double** res = new double* [row];
	for (int i = 0; i < row; i++)
		res[i] = new double[row];
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < row; j++)
		{
			res[i][j] = 0;
			for (int k = 0; k < column; k++)
				res[i][j] += (a[i * column + k] * b[k * row + j]);
		}
	}
	return res;
}



double* multiVectors(double* a, double* b, int size) {
	double* res = new double[size];
	for (int i = 0; i < size; i++)
		res[i] = a[i] * b[i];
	return res;
}

int main() {
	const int row = 4;
	const int column = 2;
	double vector1[] = { 0.,1.,0.};
	double vector2[] = { 1.,1.,1.};
	double matrix1[row][column] = { {1., 0.},
									{0., 1.},
									{0., 0.},
									{1., 1} };
	double matrix2[column][row] = { {2., 2., 2., 1},
									{2., 2., 2., 1} };
	double *resVector;
	double **resMatrix;

	/**
	for (int i = 0; i < row; i++) {
		cout << endl;
		for (int j = 0; j < column; j++)
			cout << matrix2[i][j] << "  ";
	}
	*/


	cout << "multiplication of vectors"<<endl;
	int size = sizeof(vector1) / sizeof(*(vector1));
	resVector = multiVectors(vector1, vector2, size);
	for (int i = 0; i < size; i++) {
		cout << resVector[i] << "  " ;
	}
	cout << endl;
	delete[] resVector;


	cout << "multiplication of matrixes " << endl;
	resMatrix = multiMatrixes(&matrix1[0][0], &matrix2[0][0], row, column);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < row; j++)
			cout << resMatrix[i][j] << "  ";
		cout << endl;
	}
	cout << endl;
	freeSpaceMatrix(resMatrix, row);

	return 0;
}
