#include <stdio.h>
#include <iostream>
#include <time.h>
#include <omp.h>
#include <iomanip>
#include <cstdlib>
using namespace std;
//#pragma comment(linker, "/STACK:8000000")
//#pragma comment(linker, "/HEAP:8000000")

float **res;
int N, nstep;
float eps, h;

typedef float(*func)(float x, float y);
float** SolveEquationDirikhle(func G, func F, int N, float eps);
float F(float x, float y);
float G(float x, float y);

void freeSpaceMatrix(float** res, int row) {
	if (res) {
		for (int i = 0; i < row; i++) { if (res[i]) { delete[] res[i]; } }
		delete[] res;
	}
}


//float** multiMatrixesForce(float a[], float b[], const int row, int const column) {
//	float** res = new float* [row];
//	float sum;
//	for (int i = 0; i < row; i++)
//		res[i] = new float[column];
//
//	}
//	for (int i = 0; i < row; i++) {
//		for (int j = 0; j < column; j++)
//		{
//			sum = 0.;
//			for (int k = 0; k < column; k++) {
//				cout << "i: " << i << "j: " << j << "k: "<< k<<endl;
//				sum += (aa[i][k] * bb[k][j]);
//			}
//			res[i][j] = sum;
//		}
//	}
//	return res;
//}
//




float* multiVectors(float* a, float* b, int size) {
	float* res = new float[size];
	for (int i = 0; i < size; i++)
		res[i] = a[i] * b[i];
	return res;
}

int main() {
	const int row = 2;
	const int column = 2;
	const int Brow = 100;
	const int Bcolumn = 100;
	//float *Bmatrix1[Brow];
	//float *Bmatrix2[Bcolumn];
	float vector1[] = { 0.,1.,0.};
	float vector2[] = { 1.,1.,1.};
	float matrix1[row][column] = { {1., 0.},
									{0., 1.} };
	float matrix2[column][row] = { {2., 2.},
									{2., 2.}};
	float *resVector;
	float sum;

	/**
	for (int i = 0; i < row; i++) {
		cout << endl;
		for (int j = 0; j < column; j++)
			cout << matrix2[i][j] << "  ";
	}
	*/
	// create dinamic arrays matrix1, matrix2
	float** Bmatrix1 = new float*[Brow];
	for (int i = 0; i < Brow; i++)
		Bmatrix1[i] = new float[Bcolumn];
	float** Bmatrix2 = new float* [Brow];
	for (int i = 0; i < Brow; i++)
		Bmatrix2[i] = new float[Bcolumn];
	// create dinamic array result
	float** resMatrix = new float* [Brow];
	for (int i = 0; i < Brow; i++)
		resMatrix[i] = new float[Bcolumn];
	// filling dinamic arrays
	srand(time(0));
	for (int i = 0; i < Brow; i++) {
		for (int j = 0; j < Bcolumn; j++) {
			Bmatrix1[i][j] = rand() % 10 - rand() % 10;
			Bmatrix2[i][j] = rand() % 10 - rand() % 10;
		}
	}

	//multiplication of matrixes fast
	cout << "multiplication of matrixes " << endl;
	double start_time_force = omp_get_wtime();
	for (int i = 0; i < Brow; i++) {
		for (int j = 0; j < Bcolumn; j++)
		{
			sum = 0.;
			for (int k = 0; k < column; k++) {
				//cout << "i: " << i << "j: " << j << "k: " << k << endl;
				sum += (Bmatrix1[i][k] * Bmatrix2[k][j]);
			}
			resMatrix[i][j] = sum;
		}
	}
	double end_time_force = omp_get_wtime();


	//multiplication of matrixes slow
	double start_time = omp_get_wtime();
	for (int j = 0; j < Brow; j++) {
		for (int i = 0; i < Bcolumn; i++)
		{
			sum = 0.;
			for (int k = 0; k < column; k++) {
				//cout << "i: " << i << "j: " << j << "k: " << k << endl;
				sum += (Bmatrix1[i][k] * Bmatrix2[k][j]);
			}
			resMatrix[i][j] = sum;
		}
	}
	double end_time = omp_get_wtime();
	freeSpaceMatrix(resMatrix, Bcolumn);
	freeSpaceMatrix(Bmatrix1, Bcolumn);
	freeSpaceMatrix(Bmatrix2, Bcolumn);

	//multiplication of vectors
	cout << "size of vector: ";
	int size = sizeof(vector1) / sizeof(*(vector1));
	cout << size << endl;

	cout << "multiplication of vectors"<<endl;
	resVector = multiVectors(vector1, vector2, size);
	for (int i = 0; i < size; i++) {
		cout << resVector[i] << "  " ;
	}
	cout << endl;
	delete[] resVector;


	/*double start_time_force = clock();
	multiMatrixes((float*)Bmatrix1, (float*)Bmatrix2, Brow, Bcolumn);
	double end_time_force = clock();*/

	cout << "Time general:" << end_time - start_time << endl;
	cout << "Time force:" << end_time_force - start_time_force << endl;

	/*//output result matrixes
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < row; j++)
			cout << resMatrix[i][j] << "  ";
		cout << endl;
	}
	cout << endl;
	freeSpaceMatrix(resMatrix, row);*/

	//The task Dirikhle method Gauss-Zeidal
	N = 99;
	eps = 0.0001;
	nstep = 10;

	double start_t = omp_get_wtime();
	res = SolveEquationDirikhle(&G, &F, N, eps);
	double finish_t = omp_get_wtime();
	cout << "Time = " << finish_t - start_t <<endl;

	for (int i = 0; i < N + 2; i++) {
		for (int j = 0; j < N + 2; j++) {
			//cout << "i: " << i << "j: " << j << endl;
			if (j % nstep == 0) {
				if (i % nstep == 0)
					cout << setw(10) << setprecision(3) << res[i][j];
			}
		}
		if (i % nstep == 0)
			cout << endl;
	}
	

	return 0;
}



float F(float x, float y)
{
	return sin(x) * sin(y);
}
float G(float x, float y)
{
	if (x == 0.) return 1. - 2. * y;
	if (x == 1.) return -1. + 2. * y;
	if (y == 0.) return 1. - 2. * x;
	if (y == 1.) return -1. + 2. * x;
}


float** SolveEquationDirikhle(func G, func F, int N, float eps)
{

	h = 1.0 / (N + 1);
	float const h = 1.0 / (N + 1);
	float** f = new float* [N];
	for (int i = 0; i < N; i++)
	{
		f[i] = new float[N];
		for (int j = 0; j < N; j++)
			f[i][j] = F((i + 1) * h, (j + 1) * h);
	}
	float** u = new float* [N + 2];
	for (int i = 1; i < N + 1; i++)
	{
		u[i] = new float[N + 2];
		for (int j = 1; j < N + 1; j++)
			u[i][j] = 0.;
		u[i][0] = G(i * h, 0.);
		u[i][N + 1] = G(i * h, (N + 1) * h);
	}
	u[0] = new float[N + 2];
	u[N + 1] = new float[N + 2];
	for (int j = 0; j < N + 2; j++)
	{
		u[0][j] = G(0, j * h);
		u[N + 1][j] = G((N + 1) * h, j * h);
	}

	for (int i = 0; i < N+2; i++) {
		for (int j = 0; j < N + 2 + 1; j++) {
			//cout << "i: " << i << "j: " << j << endl;
			if (j % nstep == 0) {
				if (i % nstep == 0)
					cout << setw(7) << setprecision(3) << u[i][j];
			}
			
		}
		if (i % nstep == 0)
			cout << endl;
	}

	float max;
	int Iter = 0;
	do
	{
		Iter++;
		max = 0;
		for (int i = 1; i < N + 1; i++)
			for (int j = 1; j < N + 1; j++)
			{
				float u0 = u[i][j];
				u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + 
					u[i][j + 1] - h * h * f[i - 1][j - 1]);
				float d = abs(u[i][j] - u0);
				if (d > max)
					max = d;
			}
	} while (max > eps);

	cout << "IterCnt = " << Iter << endl;

	return u;
}
	//uOutput();
	

//void uOutput()
//{
//	for (int i = 0; i < nstep; i++)
//	{
//		for (int j = 0; j < nstep; j++) {
//			cout << "i: " << i << "j: " << j << endl;
//			cout << u[ + int((N + 1) / (nstep - 1))][j + int((N + 1) / (nstep - 1))] << " ";
//		}
//		cout << endl;
//	}
//}

