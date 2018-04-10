#include "HLS/hls.h"
#include "HLS/math.h"

#include <stdio.h> //printf
#include <stdlib.h> //srand, rand

#define TEST_SEED 3
#define M 9
#define N 8

component void r8mat_svd_linpack_9_8(short m, short n, ihc::stream_in<float> &a, ihc::stream_out<float> &U, ihc::stream_out<float> &S,
	ihc::stream_out<float> &V);
// print Matrix in col major order
void printMat(short m, short n, float* A){
	for(short i = 0; i < m; i++){
		for(short j = 0; j < n; j++){
			printf("%f\t", A[i + j * m]);
		}
		printf("\n");
	}
	printf("\n");
}

int main(void) {
	// All Matrix are in col major order
	float A[M * N];
	srand(TEST_SEED);

	// Generate test inputs
	for (short i = 0; i < M * N; ++i) {
		// generate inputs between 0.0 and 1.0
		float x = (float) rand() / (((float) RAND_MAX) / 1.0f);
		A[i] = x;
	}

	printMat(9,8,A);



	// fill stream with input data and prepare output streams.
  	ihc::stream_out<float> U;
	ihc::stream_out<float> S;
	ihc::stream_out<float> V;
	ihc::stream_in<float>  A_stream;

	for (short row = 0; row < M; row++)
	{
		for (short col = 0; col < N; col++)
		{
			A_stream.write(A[row * N + col]);
		}
	}

	 // invoke the component and collect the results
	 r8mat_svd_linpack_9_8(M, N, A_stream, U, S, V);
	 printf("svd finish\n");
	 float computed_U_matrix[M*M];
	 float computed_S_matrix[M*N];
	 float computed_V_matrix[N*N];

	 for(short i = 0; i < M; i++){
		for(short j = 0; j < M; j++){
			// row major order
			computed_U_matrix[i * M + j] = U.read();
		}
	 }
	 printMat(9,9,computed_U_matrix);

	 for(short i = 0; i < M; i++){
		for(short j = 0; j < N; j++){
			// row major order
			computed_S_matrix[i * N + j] = S.read();
		}
	 }
	 printMat(9,8,computed_S_matrix);

	 for(short i = 0; i < N; i++){
		for(short j = 0; j < N; j++){
			// row major order
			computed_V_matrix[i * N + j] = V.read();
		}
	 }
	 printMat(8,8,computed_V_matrix);
	 return 0;
}