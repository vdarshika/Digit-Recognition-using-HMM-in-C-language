#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "float.h"
#include "LBG.h"

const int N = 5;
const int M = 32; // M represents number of observations per state.
int T = 0;		  // T represents number of observations.
const int C = 12;
const int S = 320;

long double A_matrix[N][N] = {0};
long double A_matrix_new[N][N] = {0};
long double A_matrix_avg[N][N] = {0};

long double B_matrix[N][M] = {0};
long double B_matrix_new[N][M] = {0};
long double B_matrix_avg[N][M] = {0};

long double Pi_matrix[N] = {0};
long double Pi_matrix_new[N] = {0};
long double Pi_matrix_avg[N] = {0};

long double alpha[120][N] = {0};
long double beta[120][N] = {0};
long double delta_i[120][N] = {0};
long double Xi[120][N][N] = {0};
long double gamma[120][N] = {0};

int q_star[120] = {0};
int psi_i[120][N] = {0};

long double P_star = 0;
long double p_new = 0;
long double p_old = 0;
long double prob_obs_seq_given_lambda = 0.0;
int Frames[10][20] = {0};
long double universe[10][20][120][12] = {0};
long double frame_array[150][S] = {0};
long double frame_array_temp[150][S] = {0};
long double reference_Ci_array[120][C + 1] = {0};
int observation_sequence[120] = {0};
int observation_sequence_test[120] = {0};

///////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// MODELS FOR EVERY DIGIT //////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
long double A_matrix_0[N][N] = {0};
long double B_matrix_0[N][M] = {0};
long double Pi_matrix_0[N] = {0};

long double A_matrix_1[N][N] = {0};
long double B_matrix_1[N][M] = {0};
long double Pi_matrix_1[N] = {0};

long double A_matrix_2[N][N] = {0};
long double B_matrix_2[N][M] = {0};
long double Pi_matrix_2[N] = {0};

long double A_matrix_3[N][N] = {0};
long double B_matrix_3[N][M] = {0};
long double Pi_matrix_3[N] = {0};

long double A_matrix_4[N][N] = {0};
long double B_matrix_4[N][M] = {0};
long double Pi_matrix_4[N] = {0};

long double A_matrix_5[N][N] = {0};
long double B_matrix_5[N][M] = {0};
long double Pi_matrix_5[N] = {0};

long double A_matrix_6[N][N] = {0};
long double B_matrix_6[N][M] = {0};
long double Pi_matrix_6[N] = {0};

long double A_matrix_7[N][N] = {0};
long double B_matrix_7[N][M] = {0};
long double Pi_matrix_7[N] = {0};

long double A_matrix_8[N][N] = {0};
long double B_matrix_8[N][M] = {0};
long double Pi_matrix_8[N] = {0};

long double A_matrix_9[N][N] = {0};
long double B_matrix_9[N][M] = {0};
long double Pi_matrix_9[N] = {0};
void initialisation_of_hmm()
{
	int index = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			A_matrix[i][j] = 0;
		}
	}
	A_matrix[0][0] = A_matrix[1][1] = A_matrix[2][2] = A_matrix[3][3] = 0.8;
	A_matrix[4][4] = 1;
	A_matrix[0][1] = A_matrix[1][2] = A_matrix[2][3] = A_matrix[3][4] = 0.2;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			B_matrix[i][j] = 0.03125;
		}
	}

	for (int i = 0; i < N; i++)
	{
		if (i == 0)
		{
			Pi_matrix[i] = 1;
		}
		else
		{
			Pi_matrix[i] = 0;
		}
	}
}

void forward_procedure(int T)
{
	int i, j;

	for (i = 0; i < T; i++)
	{
		for (j = 0; j < N; j++)
		{
			alpha[i][j] = 0;
		}
	}
	prob_obs_seq_given_lambda = 0.0;
	//INITIALIZATION
	for (i = 0; i < N; i++)
	{
		alpha[0][i] = Pi_matrix[i] * B_matrix[i][observation_sequence[0]];
	}
	//INDUCTION
	for (i = 1; i < T; i++)
	{

		for (j = 0; j < N; j++)
		{
			long double sum = 0;
			for (int k = 0; k < N; k++)
			{
				sum = sum + (alpha[i - 1][k] * A_matrix[k][j]);
			}
			alpha[i][j] = sum * B_matrix[j][observation_sequence[i]];
		}
	}

	//TERMINATION
	for (i = 0; i < N; i++)
	{
		prob_obs_seq_given_lambda = prob_obs_seq_given_lambda + alpha[T - 1][i];
	}

	//printf("P(observation_sequence/Lambda)--->%g\n", prob_obs_seq_given_lambda);

	// printf("---------THE ALPHA MATRIX----------------\n\n");
	// printf("---------THE ALPHA MATRIX----------------\n\n");

	// for (i = 0; i < T; i++)
	// {
	// 	printf("row--->%d\n",i);
	// 	for (j = 0; j < N; j++)
	// 	{
	// 		printf("%g\t\t", alpha[i][j]);
	// 	}
	// 	printf("\n");
	// }
}

long double forward_procedure_test(int T, long double A_matrix[N][N], long double B_matrix[N][M], long double Pi_matrix[N])
{
	int i, j;

	for (i = 0; i < T; i++)
	{
		for (j = 0; j < N; j++)
		{
			alpha[i][j] = 0;
		}
	}
	prob_obs_seq_given_lambda = 0.0;
	//INITIALIZATION
	for (i = 0; i < N; i++)
	{
		alpha[0][i] = Pi_matrix[i] * B_matrix[i][observation_sequence_test[0]];
	}
	//INDUCTION
	for (i = 1; i < T; i++)
	{

		for (j = 0; j < N; j++)
		{
			long double sum = 0;
			for (int k = 0; k < N; k++)
			{
				sum = sum + (alpha[i - 1][k] * A_matrix[k][j]);
			}
			alpha[i][j] = sum * B_matrix[j][observation_sequence_test[i]];
		}
	}

	//TERMINATION
	for (i = 0; i < N; i++)
	{
		prob_obs_seq_given_lambda = prob_obs_seq_given_lambda + alpha[T - 1][i];
	}
	return prob_obs_seq_given_lambda;
}

void backward_procedure(int T)
{
	int i, j;
	for (i = 0; i < T; i++)
	{
		for (j = 0; j < N; j++)
		{
			beta[i][j] = 0;
		}
	}

	//INITIALIZATION
	for (i = 0; i < N; i++)
	{
		beta[T - 1][i] = 1;
	}

	//INDUCTION
	for (i = T - 2; i >= 0; i--)
	{
		for (j = 0; j < N; j++)
		{
			long double sum = 0.0;
			for (int k = 0; k < N; k++)
			{
				sum = sum + (A_matrix[j][k] * B_matrix[k][observation_sequence[i + 1]] * beta[i + 1][k]);
			}
			beta[i][j] = sum;
		}
	}
	// printf("\n\n-------------THE BETA MATRIX----------------\n\n");
	// for (i = 0; i < T; i++)
	// {
	// 	printf("T--->%d\n", i);
	// 	for (j = 0; j < N; j++)
	// 	{
	// 		printf("%g\t", beta[i][j]);
	// 	}
	// 	printf("\n");
	// }
}
void compute_gamma(int T)
{
	//printf("-----------------------ganmmammamamamamamam----------------\n");
	memset(gamma, 0, sizeof(gamma));
	long double temp = 0.0;
	for (int k = 0; k < T; k++)
	{
		for (int i = 0; i < N; i++)
		{
			temp += alpha[k][i] * beta[k][i];
		}
		for (int i = 0; i < N; i++)
		{
			gamma[k][i] = (alpha[k][i] * beta[k][i]) / temp;
			//		printf("%g   ", gamma[k][i]);
		}
		temp = 0.0;
		//	printf("\n");
	}
}
//PROBLEM-2
void viterbi_algorithm(int T)
{
	int i = 0;
	int j = 0;
	p_new = 0;
	memset(q_star, 0, sizeof(q_star));
	for (i = 0; i < T; i++)
	{
		for (j = 0; j < N; j++)
		{
			psi_i[i][j] = 0;
			delta_i[i][j] = 0;
		}
	}

	//Initialisation
	for (j = 0; j < N; j++)
	{
		delta_i[0][j] = Pi_matrix[j] * B_matrix[j][observation_sequence[0]];
		psi_i[0][j] = 0;
	}

	// //Recursion
	// printf("\n--------------------------delta-----------------------------\n\n");
	int index_delta_max = 0;
	for (i = 1; i < T; i++)
	{
		for (j = 0; j < N; j++)
		{
			long double max_delta = 0;
			index_delta_max = 0;

			for (int k = 0; k < N; k++)
			{
				long double x;

				x = delta_i[i - 1][k] * A_matrix[k][j];

				if (x > max_delta)
				{
					max_delta = x;
					index_delta_max = k;
				}
			}
			delta_i[i][j] = max_delta * B_matrix[j][observation_sequence[i]];
			psi_i[i][j] = index_delta_max;
		}
	}

	//Termination
	int max_index_P_star = 0;

	for (i = 0; i < N; i++)
	{
		if (delta_i[T - 1][i] > p_new)
		{
			p_new = delta_i[T - 1][i];
			q_star[T - 1] = i;
		}
	}
	// printf("P*---->%g\n", p_new);
	// printf("Q*---->%d\n", q_star[T - 1]);
	p_old = p_new;
	//Backtracking
	for (i = T - 2; i >= 0; i--)
	{
		q_star[i] = psi_i[i + 1][q_star[i + 1]];
	}

	// printf("\n\nOPTIMAL STATE SEQUENCE: \n");
	// for (i = 0; i < T; i++)
	// {
	// 	printf("%d  ", q_star[i]);
	// }
	// printf("\n\n");
	// printf("\n-------------------------------------------------------------------------------------------------------------------\n\n");
}
// PROBLEM-3
void bawn_welch_algorithm(int T)
{
	int t = 0;
	int i = 0;
	int j = 0;

	//printf("\n\n---------PROBLEM 3----------------\n\n");
	memset(Xi, 0, sizeof(Xi));
	memset(A_matrix_new, 0, sizeof(A_matrix_new));
	memset(B_matrix_new, 0, sizeof(B_matrix_new));
	memset(Pi_matrix_new, 0, sizeof(Pi_matrix_new));

	long double temp = 0.0;
	long double temp2 = 0.0;
	for (t = 0; t <= T - 2; t++)
	{
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				temp += alpha[t][i] * A_matrix[i][j] * B_matrix[j][observation_sequence[t + 1]] * beta[t + 1][j];
			}
		}
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				Xi[t][i][j] = (alpha[t][i] * A_matrix[i][j] * B_matrix[j][observation_sequence[t + 1]] * beta[t + 1][j]) / temp;
			}
		}
		temp = 0.0;
	}

	compute_gamma(T);
	// printf("---------------------------------------------------------------------------------------------------------\n");
	// printf("---------------------------------------------------------------------------------------------------------\n");
	// printf("---------------------------------------------------------------------------------------------------------\n");

	// printf("----------------------------NEW MODEL VALUES------------------------------------------\n");

	// printf("############--------NEW PI VALUES----------##########\n");
	for (i = 0; i < N; i++)
	{
		Pi_matrix_new[i] = gamma[0][i];
		// printf("%g\t", Pi_matrix_new[i]);
	}

	//  printf("\n\n############--------NEW A matrix VALUES----------##########\n\n");
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			long double temp_xi = 0;
			long double temp_gamma = 0;
			for (t = 0; t < T - 1; t++)
			{
				temp_xi = temp_xi + Xi[t][i][j];
				temp_gamma = temp_gamma + gamma[t][i];
			}
			A_matrix_new[i][j] = temp_xi / temp_gamma;
			//printf("%g\t", A_matrix_new[i][j]);7
		}
		//printf("\n");
	}
	//  printf("\n\n############--------NEW B matrix VALUES----------##########\n\n");
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < M; j++)
		{
			long double temp_xi = 0;
			long double temp_xi_total = 0;
			for (t = 0; t < T; t++)
			{
				if ((observation_sequence[t]) == j)
				{
					temp_xi = temp_xi + gamma[t][i];
				}
				temp_xi_total = temp_xi_total + gamma[t][i];
			}

			B_matrix_new[i][j] = temp_xi / temp_xi_total;
			if (B_matrix_new[i][j] == 0)
			{
				B_matrix_new[i][j] = 1e-030;
			}
		}
	}
	for (i = 0; i < N; i++)
	{
		long double sum = 0;
		for (j = 0; j < M; j++)
		{
			//printf("%g   ", B_matrix_new[i][j]);
		}
	}
	// // printf("----------------------------------------------------------------------------------------------\n");
	// printf("--------------------------UPDATE THE VALUES OF LAMBDA-----------------------------------------\n");
	// printf("----------------------------------------------------------------------------------------------\n");

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			A_matrix[i][j] = A_matrix_new[i][j];
		}
	}
	for (i = 0; i < N; i++)
	{
		long double sum = 0;
		for (j = 0; j < M; j++)
		{
			B_matrix[i][j] = B_matrix_new[i][j];
		}
	}
	for (i = 0; i < N; i++)
	{
		Pi_matrix[i] = Pi_matrix_new[i];
	}

	// printf("----------------------------------------------------------------------------------------------\n");
	// printf("--------------------------STOCHASTIC-----------------------------------------\n");
	// printf("----------------------------------------------------------------------------------------------\n");

	for (int i = 0; i < N; i++)
	{
		//	printf("row-->%d------->", i);
		long double sum = 0;
		long double max = 0.0;
		int index = -1;
		for (int j = 0; j < N; j++)
		{
			if (A_matrix[i][j] > max)
			{
				max = A_matrix[i][j];
				index = j;
			}
			sum = sum + A_matrix[i][j];
		}
		if (sum < 1.0)
		{
			long double diff = 1.0 - sum;
			max = max + diff;
			A_matrix[i][index] = max;
		}
		else if (sum > 1.0)
		{
			long double diff = sum - 1.0;
			max = max - diff;
			A_matrix[i][index] = max;
		}
		else
		{
			continue;
		}
		//	printf("%lf\n", sum);
	}

	// printf("----------------------------------------------------------------------------------------------\n");
	// printf("--------------------------STOCHASTIC-----------------------------\n");
	// printf("----------------------------------------------------------------------------------------------\n");

	for (int i = 0; i < N; i++)
	{
		//	printf("row-->%d------->", i);
		long double sum = 0;
		long double max = 0.0;
		int index = -1;
		for (int j = 0; j < M; j++)
		{
			if (B_matrix[i][j] > max)
			{
				max = B_matrix[i][j];
				index = j;
			}
			sum = sum + B_matrix[i][j];
		}
		if (sum < 1.0)
		{
			long double diff = 1.0 - sum;
			max = max + diff;
			B_matrix[i][index] = max;
		}
		else if (sum > 1.0)
		{
			long double diff = sum - 1.0;
			max = max - diff;
			B_matrix[i][index] = max;
		}
		else
		{
			continue;
		}
		//printf("%lf\n", sum);
	}

	for (i = 0; i < N; i++)
	{
		long double sum = 0;
		//printf("row-->%d------->", i);
		for (j = 0; j < N; j++)
		{
			//sum=sum+B_matrix[i][j];
			// printf("%g   ", A_matrix[i][j]);
		}
		//printf("%g\n",sum);
		//printf("\n");
	}
	for (i = 0; i < N; i++)
	{
		long double sum = 0;
		//printf("row-->%d------->", i);
		for (j = 0; j < M; j++)
		{
			sum = sum + B_matrix[i][j];
			//printf("%g   ", B_matrix[i][j]);
		}
		//printf("%g\n",sum);
		//printf("\n\n");
	}
	/*
	for (i = 0; i < N; i++)
	{
		long double sum = 0;
			printf("row-->%d------->", i);
		for (j = 0; j < M; j++)
		{
			sum=sum+B_matrix[i][j];
			 //printf("%g   ", B_matrix[i][j]);
		}
		printf("%g\n",sum);
		printf("\n");
	}
	*/
}
double *Compute_Ri(double *sample_data, double *R_ref) //function to calculate the value of Ri.
{
	int i = 0, j = 0;
	double sum = 0.0;
	for (i = 0; i <= 12; i++)
	{
		sum = 0.0;
		for (j = 0; j < 320 - i; j++)
		{
			sum = sum + (sample_data[j] * sample_data[j + i]);
		}
		R_ref[i] = sum;
	}

	return R_ref;
}
double *Compute_Ai(double *A_ref, double *R) //function to calculate the value of Ai.
{
	double alpha[13][13], k[13], E[13];
	double sum = 0.0;
	int i = 0, j = 0;

	E[0] = R[0];
	for (i = 1; i <= 12; i++)
	{

		if (i == 1)
		{
			k[1] = R[1] / R[0];
		}
		else
		{
			sum = 0.0;
			for (j = 1; j <= i - 1; j++)
			{
				sum = sum + (alpha[i - 1][j] * R[i - j]);
			}
			k[i] = ((R[i] - sum) / E[i - 1]);
		}
		alpha[i][i] = k[i];
		for (int j = 1; j <= i - 1; j++)
		{
			alpha[i][j] = alpha[i - 1][j] - (k[i] * alpha[i - 1][i - j]);
		}
		E[i] = (1 - (k[i] * k[i])) * E[i - 1];
	}

	for (int j = 1; j <= 12; j++)
	{
		A_ref[j] = alpha[12][j];
	}
	return A_ref;
}

double *Compute_Ci(double *C_ref, double *A, double *R) //function to calculate the value of Ci.
{
	int i = 0, j = 0;
	double sum = 0.0;

	C_ref[0] = log(R[0] * R[0]);
	for (i = 1; i <= 12; i++)
	{
		sum = 0.0;
		for (j = 1; j <= i - 1; j++)
		{
			sum = sum + (((double)j / (double)i) * (C_ref[j] * A[i - j]));
		}
		C_ref[i] = A[i] + sum;
	}
	return C_ref;
}

double *Apply_window(double *W_ref, double *C) //function to apply the raised sine window in the Ci values calculated.
{
	int i = 0;
	double x;
	for (i = 1; i <= 12; i++)
	{
		x = (3.14 * i) / 12;
		W_ref[i] = C[i] * (1 + sin(x));
	}
	return W_ref;
}
double Compute_dc_shift() //function to calculate the DC SHIFT
{
	FILE *fp;
	char ch;
	int count = 0;
	double value = 0, sum = 0.0;
	double shift = 0.0;
	fp = fopen("dc_shift_file.txt", "r"); //opens the file that contain the value of silence

	while (!feof(fp))
	{
		fscanf(fp, "%lf", &value); //reading values from the file
		sum = sum + value;		   //adding up all the valeus of the file
	}

	fseek(fp, 0, SEEK_SET); //setting pointer to the initial position 0

	for (ch = getc(fp); ch != EOF; ch = getc(fp))
	{
		if (ch == '\n')
			count++; //counting the number of lines present the file
	}
	fclose(fp);

	shift = (sum / count); //averaging the sum out

	return shift;
}
void create_universe()
{
	FILE *fp;
	fp = fopen("universe.txt", "w");
	FILE *fp1 = fopen("number_of_frames.txt", "w");
	for (int digit = 0; digit < 10; digit++)
	{
		printf("--------------------------------------------------\n");
		printf("------------DIGIT %d --------------\n", digit);
		printf("--------------------------------------------------\n");

		for (int record = 0; record < 20; record++)
		{
			//printf("------------recording %d --------------\n",record+1);

			int i = 0, j = 0, k = 0;
			T = 0;
			char recording[20];
			long double value, max_value = 0;
			memset(frame_array, 0, sizeof(frame_array));
			memset(reference_Ci_array, 0, sizeof(reference_Ci_array));
			sprintf(recording, "digit_%d_%d.txt", digit, record + 1);
			//printf("%s\n",recording);
			//system("Recording_module.exe 3 sample.wav sample.text");
			long double dc_shift = Compute_dc_shift();

			FILE *file_pointer;
			FILE *file_pointer_temp;
			file_pointer = fopen(recording, "r");

			if (!file_pointer)
			{
				perror("Error in opening the file !!!\n");
				return;
			}
			while (!feof(file_pointer))
			{
				fscanf(file_pointer, "%Lf", &value);
				if (value > max_value)
				{
					max_value = value;
				}
			}
			fseek(file_pointer, 0, SEEK_SET);

			file_pointer_temp = fopen(recording, "r");
			while (!feof(file_pointer))
			{
				for (j = 0; j < S; j++)
				{
					fscanf(file_pointer, "%Lf", &frame_array[T][j]);
				}
				fseek(file_pointer_temp, 0, SEEK_SET);
				for (k = 0; k < (T + 1) * 100; k++)
				{
					fscanf(file_pointer_temp, "%Lf", &value);
					continue;
				}
				//printf("----------------->%lf\n",value);

				file_pointer = file_pointer_temp;
				T++;
			}

			if (T > 120)
			{
				T = 120;
			}
			//printf("NUMber of frames %d\n",T);
			Frames[digit][record] = T;
			fprintf(fp1, "%d\t", T);
			printf("%d\t", T);
			//printf(" \n\n\nThe max value = %g\n", max_value);

			for (i = 0; i < T; i++)
			{
				for (j = 0; j < S; j++)
				{
					//printf("%g   ",frame_array[i][j]);
				}
				//printf("\n\n\n");
			}
			//Normalizing the values of the file.
			for (i = 0; i < T; i++)
			{
				for (j = 0; j < S; j++)
				{
					frame_array[i][j] = ((frame_array[i][j] - dc_shift) * 5000.0) / max_value;
				}
			}
			for (i = 0; i < T; i++)
			{
				for (j = 0; j < S; j++)
				{
					//printf("%g   ",frame_array[i][j]);
				}
				//printf("\n");
			}
			double R_ref[13];
			double A_ref[13];
			double C_ref[13];
			double W_ref[13];
			for (int frame = 0; frame < T; frame++)
			{
				double temp[320];
				for (int i = 0; i < 320; i++)
				{
					temp[i] = frame_array[frame][i]; //store the 320 samples pf each frame in a sample array to calculate the Ri,Ci and Ai values
				}

				double *R = Compute_Ri(temp, R_ref); //compute teh value of Ri
				double *A = Compute_Ai(A_ref, R);	 //compute teh value of Ai
				double *C = Compute_Ci(C_ref, A, R); //compute teh value of Ci
				double *W = Apply_window(W_ref, C);	 //Apply the raised sine window on the Ci values.

				for (int sample = 1; sample <= 12; sample++)
				{
					reference_Ci_array[frame][sample] = W[sample]; //storing the ci values for 10 records 5 frames 12 sample ci values.
					universe[digit][record][frame][sample] = W[sample];
					fprintf(fp, "%Lf\t", universe[digit][record][frame][sample]);
				}
				fprintf(fp, "\n");
			}
		}
		fprintf(fp1, "\n");
	}
}
void create_codebook()
{
	int i = 0, j = 0;
	int iteration = 0;
	long double temp;
	long double centroid[12] = {0};

	FILE *file_pointer; //file pointer read the file.

	file_pointer = fopen("universe.txt", "r"); //open the file.

	//check whether the file is opened correctly or not.
	if (file_pointer == NULL)
	{
		printf("\n Error in opening the file!!");
		return;
	}

	// read the file and store the values in the 2d array data_vectors.
	for (i = 0; i < total_vectors; i++)
	{
		for (j = 0; j < total_ci; j++)
		{
			fscanf(file_pointer, "%Lf", &temp);
			data_vectors[i][j] = temp;
		}
	}

	//Calculating the optimal centroid of the entire universe by summing all the data vectors and dividing it by the total number of data vectors.
	for (i = 0; i < total_ci; i++)
	{
		for (j = 0; j < total_vectors; j++)
		{
			centroid[i] = centroid[i] + data_vectors[j][i];
		}
		centroid[i] = centroid[i] / (double)total_vectors;
	}

	//Initialisation of code_book for LBG algorithm
	for (i = 0; i < total_ci; i++)
	{
		code_book[0][i] = centroid[i];
	}

	printf("\n\n\n\t\t\t\t\t\t\t\t\t\t Darshika Verma\n\t\t\t\t\t\t\t\t\t\t 214101014\n\t\t\t\t\t\t\t\t\t\t Assignment-4\n\n\n\n");

	printf("---------------------------------------------------------------------- INITIAL  CODEBOOK OF LBG  ALGORITHM ------------------------------------------------------------------------------\n\n");

	for (i = 0; i < clusters; i++)
	{
		printf("   ");
		for (j = 0; j < total_ci; j++)
		{
			printf("%Lf\t", code_book[i][j]);
		}
		printf("\n\n");
	}

	LBG_algorithm();
	FILE *file_pointer_codebook = fopen("OPTIMAL_CODEBOOK.txt", "w");

	printf("---------------------------------------------------------------------- OPTIMAL  CODEBOOK OF LBG  ALGORITHM ----------------------------------------------------------------------\n\n");
	for (i = 0; i < clusters; i++)
	{
		//printf("  ");
		for (j = 0; j < total_ci; j++)
		{
			printf("%Lf\t", code_book[i][j]);
			fprintf(file_pointer_codebook, "%Lf\t", code_book[i][j]);
		}
		printf("\n");
		fprintf(file_pointer_codebook, "\n");
	}

	printf("\n-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n");

	printf("\tFINAL DISTORTION ------> %Lf\n\n", new_distortion);
}
void training()
{
	for (int digit = 0; digit < 10; digit++)
	{
		memset(A_matrix_avg, 0, sizeof(A_matrix_avg));
		memset(B_matrix_avg, 0, sizeof(B_matrix_avg));
		memset(Pi_matrix_avg, 0, sizeof(Pi_matrix_avg));
		for (int record = 0; record < 20; record++)
		{
			//printf("\n---->recording %d\n\n", record + 1);
			long double tokhura_weights[13] = {0.0, 1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
			long double distance[33] = {0};
			long double tokhura_distance = INT_MAX, final_distance = 0.0;
			int obs_index = -1;
			for (int frame = 0; frame < Frames[digit][record]; frame++)
			{
				tokhura_distance = INT_MAX;
				for (int obs = 0; obs < M; obs++)
				{
					distance[obs] = 0;
					final_distance = 0.0;
					for (int index = 1; index <= C; index++)
					{
						final_distance = final_distance + tokhura_weights[index] * ((code_book[obs][index] - universe[digit][record][frame][index]) * (code_book[obs][index] - universe[digit][record][frame][index]));
					}
					distance[obs] = final_distance;
					if (final_distance < tokhura_distance)
					{
						tokhura_distance = final_distance;
						obs_index = obs;
					}
				}
				observation_sequence[frame] = obs_index;
			}
			// printf("------------------------observations sequence--------------\n");
			// for (int i = 0; i < Frames[digit][record]; i++)
			// {
			// 	printf("%d   ", observation_sequence[i]);
			// }

			initialisation_of_hmm(); //Initialisation of HMM with the bakis model
			p_old = 0;
			p_new = 0;
			P_star = 0;
			for (int i = 0; i < 100; i++)
			{
				forward_procedure(Frames[digit][record]);
				backward_procedure(Frames[digit][record]);
				viterbi_algorithm(Frames[digit][record]);
				if (p_new < p_old)
				{
					P_star = p_old;
					break;
				}
				bawn_welch_algorithm(Frames[digit][record]);
			}
			// printf("\n\nP(observation_sequence/Lambda)--->%g\n", prob_obs_seq_given_lambda);
			// printf("P*---->%g\n", p_old);
			// printf("Q*---->%d\n", q_star[Frames[digit][record] - 1]);

			// printf("\n\nOPTIMAL STATE SEQUENCE: \n");
			// for (int i = 0;i < Frames[digit][record]; i++)
			// {
			// 	printf("%d  ", q_star[i]);
			// }
			// printf("\n\n");
			// printf("\n-------------------------------------------------------------------------------------------------------------------\n\n");

			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					A_matrix_avg[i][j] = A_matrix_avg[i][j] + A_matrix[i][j];
				}
			}
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < M; j++)
				{
					B_matrix_avg[i][j] = B_matrix_avg[i][j] + B_matrix[i][j];
				}
			}

			for (int i = 0; i < N; i++)
			{
				Pi_matrix_avg[i] = Pi_matrix_avg[i] + Pi_matrix[i];
			}
		}

		// printf("\n-----------------------------------------------------------------------\n");
		// printf("\n-----------------AVERAGED MODEL %d------------------------\n", digit);
		// printf("-----------------------------------------------------------------------\n");
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				A_matrix_avg[i][j] = A_matrix_avg[i][j] / 20.0;
				// printf("%g   ", A_matrix_avg[i][j]);
			}
			// printf("\n");
		}
		// printf("\n\n");
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				B_matrix_avg[i][j] = B_matrix_avg[i][j] / 20.0;
				// printf("%g   ", B_matrix_avg[i][j]);
			}
			//	printf("\n\n\n");
		}
		//printf("\n\n");
		for (int i = 0; i < N; i++)
		{
			Pi_matrix_avg[i] = Pi_matrix_avg[i] / 20.0;
			//	printf("%g    ", Pi_matrix_avg[i]);
		}
		//printf("\n\n");

		for (int itr = 0; itr < 4; itr++)
		{
			// printf("------------------------------------------------------\n");
			// printf("-------iteration   %d -----------------\n", itr);
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					A_matrix[i][j] = A_matrix_avg[i][j];
					//printf("%g   ", A_matrix[i][j]);
				}
				//printf("\n");
			}
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < M; j++)
				{
					B_matrix[i][j] = B_matrix_avg[i][j];
					//printf("%g   ", B_matrix[i][j]);
				}
				//printf("\n\n");
			}
			for (int i = 0; i < N; i++)
			{
				Pi_matrix[i] = Pi_matrix_avg[i];
				//printf("%g    ", Pi_matrix[i]);
			}

			memset(A_matrix_avg, 0, sizeof(A_matrix_avg));
			memset(B_matrix_avg, 0, sizeof(B_matrix_avg));
			memset(Pi_matrix_avg, 0, sizeof(Pi_matrix_avg));
			for (int record = 0; record < 20; record++)
			{
				//printf("\n---->recording %d\n\n", record + 1);
				long double tokhura_weights[13] = {0.0, 1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
				long double distance[33] = {0};
				long double tokhura_distance = INT_MAX, final_distance = 0.0;
				int obs_index = -1;
				for (int frame = 0; frame < Frames[digit][record]; frame++)
				{
					tokhura_distance = INT_MAX;
					for (int obs = 0; obs < M; obs++)
					{
						distance[obs] = 0;
						final_distance = 0.0;
						for (int index = 1; index <= C; index++)
						{
							final_distance = final_distance + tokhura_weights[index] * ((code_book[obs][index] - universe[digit][record][frame][index]) * (code_book[obs][index] - universe[digit][record][frame][index]));
						}
						distance[obs] = final_distance;
						if (final_distance < tokhura_distance)
						{
							tokhura_distance = final_distance;
							obs_index = obs;
						}
					}
					observation_sequence[frame] = obs_index;
				}
				// printf("------------------------observations sequence--------------\n");
				// for (int i = 0; i < Frames[digit][record]; i++)
				// {
				// 	printf("%d   ", observation_sequence[i]);
				// }
				p_old = 0;
				p_new = 0;
				P_star = 0;
				for (int i = 0; i < 100; i++)
				{
					forward_procedure(Frames[digit][record]);
					backward_procedure(Frames[digit][record]);
					viterbi_algorithm(Frames[digit][record]);
					if (p_new <= p_old)
					{
						P_star = p_old;
						break;
					}
					bawn_welch_algorithm(Frames[digit][record]);
				}
				// printf("\n\nP(observation_sequence/Lambda)--->%g\n", prob_obs_seq_given_lambda);
				// printf("P*---->%g\n", P_star);
				// printf("Q*---->%d\n", q_star[Frames[digit][record] - 1]);

				// printf("\n\nOPTIMAL STATE SEQUENCE: \n");
				// for (int i = 0; i < Frames[digit][record]; i++)
				// {
				// 	printf("%d  ", q_star[i]);
				// }
				// printf("\n\n");
				// printf("\n-------------------------------------------------------------------------------------------------------------------\n\n");

				for (int i = 0; i < N; i++)
				{
					for (int j = 0; j < N; j++)
					{
						A_matrix_avg[i][j] = A_matrix_avg[i][j] + A_matrix[i][j];
					}
				}
				for (int i = 0; i < N; i++)
				{
					for (int j = 0; j < M; j++)
					{
						B_matrix_avg[i][j] = B_matrix_avg[i][j] + B_matrix[i][j];
					}
				}

				for (int i = 0; i < N; i++)
				{
					Pi_matrix_avg[i] = Pi_matrix_avg[i] + Pi_matrix[i];
				}
			}

			// printf("\n-----------------------------------------------------------------------\n");
			// printf("\n-----------------AVERAGED MODEL %d------------------------\n", digit);
			// printf("-----------------------------------------------------------------------\n");
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					A_matrix_avg[i][j] = A_matrix_avg[i][j] / 20.0;
					//printf("%g   ", A_matrix_avg[i][j]);
				}
				//printf("\n");
			}
			//printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < M; j++)
				{
					B_matrix_avg[i][j] = B_matrix_avg[i][j] / 20.0;
					//printf("%g   ", B_matrix_avg[i][j]);
				}
				//printf("\n\n\n");
			}
			//printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				Pi_matrix_avg[i] = Pi_matrix_avg[i] / 20.0;
				//printf("%g    ", Pi_matrix_avg[i]);
			}
			//printf("\n\n");
		}

		printf("\n-----------------------------------------------------------------------\n");
		printf("\n-----------------FINAL MODEL %d------------------------\n", digit);
		printf("-----------------------------------------------------------------------\n");

		if (digit == 0)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					A_matrix_0[i][j] = A_matrix_avg[i][j];
					printf("%g    ", A_matrix_0[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < M; j++)
				{
					B_matrix_0[i][j] = B_matrix_avg[i][j];
					printf("%g    ", B_matrix_0[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				Pi_matrix_0[i] = Pi_matrix_avg[i];
				printf("%g    ", Pi_matrix_0[i]);
			}
		}
		else if (digit == 1)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					A_matrix_1[i][j] = A_matrix_avg[i][j];
					printf("%g    ", A_matrix_0[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < M; j++)
				{
					B_matrix_1[i][j] = B_matrix_avg[i][j];
					printf("%g    ", B_matrix_1[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				Pi_matrix_1[i] = Pi_matrix_avg[i];
				printf("%g    ", Pi_matrix_1[i]);
			}
		}
		else if (digit == 2)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					A_matrix_2[i][j] = A_matrix_avg[i][j];
					printf("%g    ", A_matrix_2[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < M; j++)
				{
					B_matrix_2[i][j] = B_matrix_avg[i][j];
					printf("%g    ", B_matrix_2[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				Pi_matrix_2[i] = Pi_matrix_avg[i];
				printf("%g    ", Pi_matrix_2[i]);
			}
		}
		else if (digit == 3)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					A_matrix_3[i][j] = A_matrix_avg[i][j];
					printf("%g    ", A_matrix_3[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < M; j++)
				{
					B_matrix_3[i][j] = B_matrix_avg[i][j];
					printf("%g    ", B_matrix_3[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				Pi_matrix_3[i] = Pi_matrix_avg[i];
				printf("%g    ", Pi_matrix_3[i]);
			}
		}
		else if (digit == 4)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					A_matrix_4[i][j] = A_matrix_avg[i][j];
					printf("%g    ", A_matrix_4[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < M; j++)
				{
					B_matrix_4[i][j] = B_matrix_avg[i][j];
					printf("%g    ", B_matrix_4[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				Pi_matrix_4[i] = Pi_matrix_avg[i];
				printf("%g    ", Pi_matrix_4[i]);
			}
		}
		else if (digit == 5)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					A_matrix_5[i][j] = A_matrix_avg[i][j];
					printf("%g    ", A_matrix_5[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < M; j++)
				{
					B_matrix_5[i][j] = B_matrix_avg[i][j];
					printf("%g    ", B_matrix_5[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				Pi_matrix_5[i] = Pi_matrix_avg[i];
				printf("%g    ", Pi_matrix_5[i]);
			}
		}
		else if (digit == 6)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					A_matrix_6[i][j] = A_matrix_avg[i][j];
					printf("%g    ", A_matrix_6[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < M; j++)
				{
					B_matrix_6[i][j] = B_matrix_avg[i][j];
					printf("%g    ", B_matrix_6[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				Pi_matrix_6[i] = Pi_matrix_avg[i];
				printf("%g    ", Pi_matrix_6[i]);
			}
		}
		else if (digit == 7)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					A_matrix_7[i][j] = A_matrix_avg[i][j];
					printf("%g    ", A_matrix_7[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < M; j++)
				{
					B_matrix_7[i][j] = B_matrix_avg[i][j];
					printf("%g    ", B_matrix_7[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				Pi_matrix_7[i] = Pi_matrix_avg[i];
				printf("%g    ", Pi_matrix_7[i]);
			}
		}
		else if (digit == 8)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					A_matrix_8[i][j] = A_matrix_avg[i][j];
					printf("%g    ", A_matrix_9[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < M; j++)
				{
					B_matrix_8[i][j] = B_matrix_avg[i][j];
					printf("%g    ", B_matrix_8[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				Pi_matrix_8[i] = Pi_matrix_avg[i];
				printf("%g    ", Pi_matrix_8[i]);
			}
		}
		else if (digit == 9)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					A_matrix_9[i][j] = A_matrix_avg[i][j];
					printf("%g    ", A_matrix_9[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < M; j++)
				{
					B_matrix_9[i][j] = B_matrix_avg[i][j];
					printf("%g    ", B_matrix_9[i][j]);
				}
				printf("\n\n");
			}
			printf("\n\n");
			for (int i = 0; i < N; i++)
			{
				Pi_matrix_9[i] = Pi_matrix_avg[i];
				printf("%g    ", Pi_matrix_9[i]);
			}
		}
		printf("\nP* for model____%d ---->%g\n", digit, P_star);
	}
}

void testing_and_accuracy()
{
	int i = 0, j = 0, k = 0;
	int total_accuracy = 0;
	for (int digit = 0; digit < 1; digit++)
	{
		int accuracy = 0;
		for (int record = 0; record < 1; record++)
		{
			system("Recording_module.exe 3 sample.wav sample.text");
			FILE *fp = fopen("o.txt", "r");
			FILE *fp1 = fopen("o.txt", "r");
			long double array[150][320] = {0};
			long double value = 0;
			long double max_value = 0.0;
			while (!feof(fp))
			{
				fscanf(fp, "%Lf", &value);
				if (value > max_value)
				{
					max_value = value;
				}
			}
			fseek(fp, 0, SEEK_SET);
			printf("max-->%Lf\n", max_value);
			while (!feof(fp))
			{
				fscanf(fp, "%Lf", &value);
				if (value > 200)
				{
					break;
				}
			}
			fp1 = fp;
			int num_frame = 0;
			for (int j = 0; j < 120; j++)
			{
				for (int i = 0; i < 320; i++)
				{
					fscanf(fp, "%Lf", &value);
					array[j][i] = value;
					if (feof(fp))
					{
						num_frame = j;
						break;
					}
				}
				if (feof(fp))
				{
					break;
				}
				for (int k = 0; k < 80; k++)
				{
					fscanf(fp1, "%Lf", &value);
					continue;
				}
				fp = fp1;
			}
			printf("--------%d-------\n", num_frame);
			long double dc_shift = Compute_dc_shift();
			for (i = 0; i < num_frame; i++)
			{
				for (j = 0; j < S; j++)
				{
					array[i][j] = ((array[i][j] - dc_shift) * 5000.0) / max_value;
				}
			}
			double R_ref[13];
			double A_ref[13];
			double C_ref[13];
			double W_ref[13];
			long double reference[150][13] = {0};
			for (i = 0; i < num_frame; i++)
			{
				double temp[320];
				for (int j = 0; j < 320; j++)
				{
					temp[i] = array[i][j]; //store the 320 samples pf each frame in a sample array to calculate the Ri,Ci and Ai values
				}
				double *R = Compute_Ri(temp, R_ref); //compute teh value of Ri
				double *A = Compute_Ai(A_ref, R);	 //compute teh value of Ai
				double *C = Compute_Ci(C_ref, A, R); //compute teh value of Ci
				double *W = Apply_window(W_ref, C);	 //Apply the raised sine window on the Ci values.

				for (int sample = 1; sample <= 12; sample++)
				{
					reference[i][sample] = W[sample];
				}
			}
			for (int i = 0; i < num_frame; i++)
			{
				for (int j = 1; j <= 12; j++)
				{
					printf("%Lf   ", reference[i][j]);
				}
				printf("\n\n");
			}

			long double tokhura_weights[13] = {0.0, 1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
			long double distance[33] = {0};
			long double tokhura_distance = INT_MAX, final_distance = 0.0;
			int obs_index = -1;
			for (int frame = 0; frame < num_frame; frame++)
			{
				tokhura_distance = INT_MAX;
				for (int obs = 0; obs < M; obs++)
				{
					distance[obs] = 0;
					final_distance = 0.0;
					for (int index = 1; index <= C; index++)
					{
						final_distance = final_distance + tokhura_weights[index] * ((code_book[obs][index] - reference[frame][index]) * (code_book[obs][index] - reference[frame][index]));
					}
					distance[obs] = final_distance;
					if (final_distance < tokhura_distance)
					{
						tokhura_distance = final_distance;
						obs_index = obs;
					}
				}
				observation_sequence_test[frame] = obs_index;
			}
			printf("----------------------observation sequence--------------------------------\n");
			for (int i = 0; i < num_frame; i++)
			{
				printf("%d    ", observation_sequence_test[i]);
			}
			printf("\n\n");

			long double probability_models[10] = {0};
			for (int i = 0; i < 10; i++)
			{
				if (i == 0)
				{
					probability_models[i] = forward_procedure_test(Frames[0][record], A_matrix_0, B_matrix_0, Pi_matrix_0);
					// printf("prob_0----> %g\n", probability_models[i]);
				}
				else if (i == 1)
				{
					probability_models[i] = forward_procedure_test(Frames[0][record], A_matrix_1, B_matrix_1, Pi_matrix_1);
					// printf("prob_1----> %g\n", probability_models[i]);
				}
				else if (i == 2)
				{
					probability_models[i] = forward_procedure_test(Frames[0][record], A_matrix_2, B_matrix_2, Pi_matrix_2);
					// printf("prob_2----> %g\n", probability_models[i]);
				}
				else if (i == 3)
				{
					probability_models[i] = forward_procedure_test(Frames[0][record], A_matrix_3, B_matrix_1, Pi_matrix_3);
					// printf("prob_3----> %g\n", probability_models[i]);
				}
				else if (i == 4)
				{
					probability_models[i] = forward_procedure_test(Frames[0][record], A_matrix_4, B_matrix_4, Pi_matrix_4);
					// printf("prob_4----> %g\n", probability_models[i]);
				}
				else if (i == 5)
				{
					probability_models[i] = forward_procedure_test(Frames[0][record], A_matrix_5, B_matrix_5, Pi_matrix_5);
					// printf("prob_5----> %g\n", probability_models[i]);
				}
				else if (i == 6)
				{
					probability_models[i] = forward_procedure_test(Frames[0][record], A_matrix_6, B_matrix_6, Pi_matrix_6);
					// printf("prob_6----> %g\n", probability_models[i]);
				}
				else if (i == 7)
				{
					probability_models[i] = forward_procedure_test(Frames[0][record], A_matrix_7, B_matrix_7, Pi_matrix_7);
					// printf("prob_7----> %g\n", probability_models[i]);
				}
				else if (i == 8)
				{
					probability_models[i] = forward_procedure_test(Frames[0][record], A_matrix_8, B_matrix_8, Pi_matrix_8);
					// printf("prob_8----> %g\n", probability_models[i]);
				}
				else if (i == 9)
				{
					probability_models[i] = forward_procedure_test(Frames[0][record], A_matrix_9, B_matrix_9, Pi_matrix_9);
					// printf("prob_9----> %g\n", probability_models[i]);
				}
			}
			long double max_prob = 0.0;
			int detect_digit = -1;
			for (int i = 0; i < 10; i++)
			{
				if (probability_models[i] > max_prob)
				{
					max_prob = probability_models[i];
					detect_digit = i;
				}
			}
			printf("Detected digit---->  %d\n", detect_digit);
			if (detect_digit == digit)
			{
				accuracy++;
				total_accuracy++;
			}
		}
		//printf("\n\n ACCURACY OF DIGIT %d------>%lf\n", digit, accuracy * 10.0);
	}
	//printf("\n\nTOTAL  ACCURACY ----->%d\n", total_accuracy);
}
int _tmain(int argc, _TCHAR *argv[])
{
	create_universe();
	//create_codebook();
	// FILE *fp, *fp1, *fp2;
	// long double values = 0.0;
	// int temp = 0;
	// fp = fopen("universe.txt", "r");
	// for (int digit = 0; digit < 10; digit++)
	// {
	// 	for (int record = 0; record < 20; record++)
	// 	{
	// 		for (int frame = 0; frame < Frames[digit][record]; frame++)
	// 		{
	// 			for (int index = 1; index <= C; index++)
	// 			{
	// 				fscanf(fp, "%Lf", &universe[digit][record][frame][index]);
	// 			}
	// 		}
	// 	}
	// }

	// fp1 = fopen("codebook_ex.txt", "r");
	// for (int i = 0; i < 32; i++)
	// {
	// 	for (int j = 1; j <= 12; j++)
	// 	{
	// 		fscanf(fp1, "%Lf", &code_book[i][j]);
	// 	}
	// }

	// fp2 = fopen("number_of_frames", "r");
	// for (int i = 0; i < 10; i++)
	// {
	// 	for (int j = 0; j < 20; j++)
	// 	{
	// 		fscanf(fp1, "%d", &Frames[i][j]);
	// 		printf("%d   ", Frames[i][j]);
	// 	}
	// 	printf("\n");
	// }
	// training();
	// testing_and_accuracy();

	system("PAUSE");
	return 0;
}