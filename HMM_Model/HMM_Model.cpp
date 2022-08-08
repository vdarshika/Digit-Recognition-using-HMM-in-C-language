#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "float.h"
#include "LBG.h"

const int N = 5;  // N represents number of states
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

long double A_matrix_test[N][N] = {0};
long double B_matrix_test[N][M] = {0};
long double Pi_matrix_test[N] = {0};

//Load the inital model of the Lambda
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

// Function to calulate the probability of the observation sequence given the model.
void forward_procedure(int T)
{
	int i, j;

	// Reloading the model with 0 values.
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
}

// Function to calulate the probability of the observation sequence given the model for the testing data.
long double forward_procedure_test(int T, int digit)
{
	int i, j;
	FILE *file_pointer_A;
	FILE *file_pointer_B;
	FILE *file_pointer_Pi;
	char buff_A[500];
	char buff_B[500];
	char buff_Pi[500];
	sprintf(buff_A, "Final_model\\\\Digit_%d\\\\A_matrix.txt", digit);
	sprintf(buff_B, "Final_model\\\\Digit_%d\\\\B_matrix.txt", digit);
	sprintf(buff_Pi, "Final_model\\\\Digit_%d\\\\Pi_matrix.txt", digit);

	memset(A_matrix_test, 0, sizeof(A_matrix_test));
	memset(B_matrix_test, 0, sizeof(B_matrix_test));
	memset(Pi_matrix_test, 0, sizeof(Pi_matrix_test));

	long double temp = 0.0;

	// Loading the final model inorder to calculate the probabilities.
	file_pointer_A = fopen(buff_A, "r");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			fscanf(file_pointer_A, "%Lf", &A_matrix_test[i][j]);
		}
	}
	fclose(file_pointer_A);

	file_pointer_B = fopen(buff_B, "r");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			fscanf(file_pointer_B, "%Lf", &B_matrix_test[i][j]);
		}
	}
	fclose(file_pointer_B);

	file_pointer_Pi = fopen(buff_Pi, "r");
	for (int i = 0; i < N; i++)
	{
		fscanf(file_pointer_Pi, "%Lf", &Pi_matrix_test[i]);
	}
	fclose(file_pointer_A);

	// Reloading the model with 0 values.
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
		alpha[0][i] = Pi_matrix_test[i] * B_matrix_test[i][observation_sequence_test[0]];
	}
	//INDUCTION
	for (i = 1; i < T; i++)
	{

		for (j = 0; j < N; j++)
		{
			long double sum = 0;
			for (int k = 0; k < N; k++)
			{
				sum = sum + (alpha[i - 1][k] * A_matrix_test[k][j]);
			}
			alpha[i][j] = sum * B_matrix_test[j][observation_sequence_test[i]];
		}
	}

	//TERMINATION
	for (i = 0; i < N; i++)
	{
		prob_obs_seq_given_lambda = prob_obs_seq_given_lambda + alpha[T - 1][i];
	}
	return prob_obs_seq_given_lambda;
}

// Computation of the beta matrix.
void backward_procedure(int T)
{
	int i, j;

	// Reloading the model with 0 values.
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
}

// Computation of the gamma matrix using in the solution of problem 3.
void compute_gamma(int T)
{
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
		}
		temp = 0.0;
	}
}

//Solution to PROBLEM-2
void viterbi_algorithm(int T)
{
	int i = 0;
	int j = 0;

	P_star = 0;

	//Reloading the matrices with 0 values.
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
		if (delta_i[T - 1][i] > P_star)
		{
			P_star = delta_i[T - 1][i];
			q_star[T - 1] = i;
		}
	}

	//Backtracking
	for (i = T - 2; i >= 0; i--)
	{
		q_star[i] = psi_i[i + 1][q_star[i + 1]];
	}
}

//Solution to  PROBLEM-3
void bawn_welch_algorithm(int T)
{
	int t = 0;
	int i = 0;
	int j = 0;

	//Reloading the matrices with values 0.
	memset(Xi, 0, sizeof(Xi));
	memset(A_matrix_new, 0, sizeof(A_matrix_new));
	memset(B_matrix_new, 0, sizeof(B_matrix_new));
	memset(Pi_matrix_new, 0, sizeof(Pi_matrix_new));

	long double temp = 0.0;
	long double temp2 = 0.0;

	// Computation of the Xi matrix used in problem 3.
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

	compute_gamma(T); // Compute tha gamma matrix.

	// Computation of the new model lambda(A,B,Pi)
	// Compute Pi_Bar.
	for (i = 0; i < N; i++)
	{
		Pi_matrix_new[i] = gamma[0][i];
	}

	// Compute A_Bar.
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
		}
	}

	// Compute B_Bar.
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

	//--------------------------UPDATE THE VALUES OF LAMBDA ---------------------------------------

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

	// Check whether the matrices are STochastic or not.
	for (int i = 0; i < N; i++)
	{
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
	}

	for (int i = 0; i < N; i++)
	{
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
	}
}

//function to calculate the value of Ri.
double *Compute_Ri(double *sample_data, double *R_ref)
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

//function to calculate the value of Ai.
double *Compute_Ai(double *A_ref, double *R)
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

//function to calculate the value of Ci.
double *Compute_Ci(double *C_ref, double *A, double *R)
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

//function to apply the raised sine window in the Ci values calculated.
double *Apply_window(double *W_ref, double *C)
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

//function to calculate the DC SHIFT
double Compute_dc_shift()
{
	FILE *fp;
	char ch;
	int count = 0;
	double value = 0, sum = 0.0;
	double shift = 0.0;
	fp = fopen("dc_shift_file.txt", "r"); //opens the file that contain the value of silence
	if (!fp)
	{
		printf("error at 686\n");
	}

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

// Function to create the universe for generation of codebook.
void create_universe()
{
	FILE *fp;
	fp = fopen("Final_Universe.txt", "w+");
	if (!fp)
	{
		printf("error at 714\n");
	}
	FILE *fp1 = fopen("Final_Frames.txt", "w+");
	if (!fp1)
	{
		printf("error at 686\n");
	}
	for (int digit = 0; digit < 10; digit++)
	{
		
		printf("--------------------------------------------------\n");
		printf("------------DIGIT %d --------------\n", digit);
		printf("--------------------------------------------------\n");

		for (int record = 0; record < 20; record++)
		{
			int i = 0, j = 0, k = 0;
			T = 0;
			char recording[500];
			long double value, max_value = 0;
			memset(frame_array, 0, sizeof(frame_array));
			memset(reference_Ci_array, 0, sizeof(reference_Ci_array));
			sprintf(recording, "Digit_%d\\\\digit_%d_%d.txt", digit, digit, record + 1);

			long double dc_shift = Compute_dc_shift();

			FILE *file_pointer;
			FILE *file_pointer_temp;
			file_pointer = fopen(recording, "r");
			if (!file_pointer)
			{
				printf("error at 747\n");
			}

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
			if (!file_pointer_temp)
			{
				printf("error at 768\n");
			}
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

				file_pointer = file_pointer_temp;
				T++;
			}

			if (T > 120)
			{
				T = 120;
			}
			Frames[digit][record] = T;
			fprintf(fp1, "%d\t", T);

			//Normalizing the values of the file.
			for (i = 0; i < T; i++)
			{
				for (j = 0; j < S; j++)
				{
					frame_array[i][j] = ((frame_array[i][j] - dc_shift) * 5000.0) / max_value;
				}
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
					//printf("%Lf\t",universe[digit][record][frame][sample]);
				}
				fprintf(fp, "\n");
				//printf("\n");
			}
		}
		fprintf(fp1, "\n");
	}
	fclose(fp);
	fclose(fp1);
}

// Function for generation of codebook.
void create_codebook()
{
	int i = 0, j = 0;
	int iteration = 0;
	long double temp;
	long double centroid[12] = {0};

	FILE *file_pointer; //file pointer read the file.

	file_pointer = fopen("Final_Universe.txt", "r"); //open the file.
	if (!file_pointer)
	{
		printf("error at 864\n");
	}

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
	fclose(file_pointer);
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
	FILE *file_pointer_codebook = fopen("Final_Codebook.txt", "w+");
	if (!file_pointer_codebook)
	{
		printf("error at 918\n");
	}

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
	fclose(file_pointer_codebook);
	printf("\n-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n");

	printf("\tFINAL DISTORTION ------> %Lf\n\n", new_distortion);
}

// Training the model for digit recognition.
void training()
{
	for (int digit = 0; digit < 10; digit++)
	{
		// Files to save the final model.
		FILE *file_pointer_A;
		FILE *file_pointer_B;
		FILE *file_pointer_Pi;
		char buff_A[500];
		char buff_B[500];
		char buff_Pi[500];
		sprintf(buff_A, "Final_model\\\\Digit_%d\\\\A_matrix.txt", digit);
		sprintf(buff_B, "Final_model\\\\Digit_%d\\\\B_matrix.txt", digit);
		sprintf(buff_Pi, "Final_model\\\\Digit_%d\\\\Pi_matrix.txt", digit);
		

		memset(A_matrix_avg, 0, sizeof(A_matrix_avg));
		memset(B_matrix_avg, 0, sizeof(B_matrix_avg));
		memset(Pi_matrix_avg, 0, sizeof(Pi_matrix_avg));
		for (int record = 0; record < 20; record++)
		{
			char x[200];
		sprintf(x,"obs_%d_%d.txt",digit,record+1);
		//printf("%s",x);
		FILE *fpo=fopen(x,"w+");
			p_old = 0;

			// Code to calculate the obervation sequence for the file.
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
				fprintf(fpo,"%d  ",(observation_sequence[frame]+1));
			}
			fprintf(fpo,"\n");
			fclose(fpo);

			initialisation_of_hmm(); //Initialisation of HMM with the bakis model

			//execute the code for problem 1,2 and 3.
			for (int i = 0; i < 100; i++)
			{
				forward_procedure(Frames[digit][record]);
				backward_procedure(Frames[digit][record]);
				viterbi_algorithm(Frames[digit][record]);
				if (P_star < p_old)
				{
					P_star = p_old;
					break;
				}
				p_old = P_star;
				bawn_welch_algorithm(Frames[digit][record]);
			}

			//Summation of all the models for the same digit.
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
		//fprintf(fpo,"\n");
		
		//-----------------AVERAGED MODEL ------------------------
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				A_matrix_avg[i][j] = A_matrix_avg[i][j] / 20.0;
			}
		}
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				B_matrix_avg[i][j] = B_matrix_avg[i][j] / 20.0;
			}
		}
		for (int i = 0; i < N; i++)
		{
			Pi_matrix_avg[i] = Pi_matrix_avg[i] / 20.0;
		}

		// printf("-------------------------------------converge the model for better accuracy for recognition----------------------------------------------------------\n");
		for (int itr = 0; itr < 4; itr++)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					A_matrix[i][j] = A_matrix_avg[i][j];
				}
			}
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < M; j++)
				{
					B_matrix[i][j] = B_matrix_avg[i][j];
				}
			}
			for (int i = 0; i < N; i++)
			{
				Pi_matrix[i] = Pi_matrix_avg[i];
			}

			memset(A_matrix_avg, 0, sizeof(A_matrix_avg));
			memset(B_matrix_avg, 0, sizeof(B_matrix_avg));
			memset(Pi_matrix_avg, 0, sizeof(Pi_matrix_avg));
			for (int record = 0; record < 20; record++)
			{

				//Calculate the observation sequence.
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

				for (int i = 0; i < 100; i++)
				{

					forward_procedure(Frames[digit][record]);
					backward_procedure(Frames[digit][record]);
					viterbi_algorithm(Frames[digit][record]);
					if (P_star < p_old)
					{
						P_star = p_old;
						break;
					}
					p_old = P_star;
					bawn_welch_algorithm(Frames[digit][record]);
				}

				//Summation of the models for same digit.
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
			//------------------Compute the averaged models
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					A_matrix_avg[i][j] = A_matrix_avg[i][j] / 20.0;
				}
			}
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < M; j++)
				{
					B_matrix_avg[i][j] = B_matrix_avg[i][j] / 20.0;
				}
			}
			for (int i = 0; i < N; i++)
			{
				Pi_matrix_avg[i] = Pi_matrix_avg[i] / 20.0;
			}
		}

		//-----------------FINAL MODEL %d------------------------

		file_pointer_A = fopen(buff_A, "w");
		file_pointer_B = fopen(buff_B, "w");
		file_pointer_Pi = fopen(buff_Pi, "w");

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				fprintf(file_pointer_A, "%g    ", A_matrix_avg[i][j]);
			}
			fprintf(file_pointer_A, "\n");
		}
		fclose(file_pointer_A);

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				fprintf(file_pointer_B, "%g    ", B_matrix_avg[i][j]);
			}
			fprintf(file_pointer_B, "\n");
		}
		fclose(file_pointer_B);

		for (int i = 0; i < N; i++)
		{
			fprintf(file_pointer_Pi, "%g    ", Pi_matrix_avg[i]);
		}
		fclose(file_pointer_Pi);

		printf("\n\tP* for Model for Digit %d -------> %g\n", digit, P_star);
	}
}

//Testing of the recorded files.
void testing_and_accuracy_recorded_files()
{
	printf("\n\n-------RESULT OF THE TESTING OF RECORDED FILES--------------\n\n\n");
	int total_accuracy = 0;
	for (int digit = 0; digit < 10; digit++)
	{
		int accuracy = 0;
		char buff[500];
		sprintf(buff, "Test_Probability\\\\digit_%d.txt", digit);
		FILE *fp = fopen(buff, "w+");
		for (int record = 11; record < 20; record++)
		{
			int x = record;
			if (digit == 8)
			{
				record = 0;
			}
			long double tokhura_weights[13] = {0.0, 1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
			long double tokhura_distance = INT_MAX, final_distance = 0.0;
			int obs_index = -1;
			memset(observation_sequence_test, 0, sizeof(observation_sequence_test));
			for (int frame = 0; frame < Frames[digit][record]; frame++)
			{
				tokhura_distance = INT_MAX;
				for (int obs = 0; obs < M; obs++)
				{
					//compute the observation sequence for the recorded file.

					final_distance = 0.0;
					for (int index = 1; index <= C; index++)
					{
						final_distance = final_distance + tokhura_weights[index] * ((code_book[obs][index] - universe[digit][record][frame][index]) * (code_book[obs][index] - universe[digit][record][frame][index]));
					}
					if (final_distance < tokhura_distance)
					{
						tokhura_distance = final_distance;
						obs_index = obs;
					}
				}
				observation_sequence_test[frame] = obs_index;
			}

			// Compute the probability of the observation sequence given the model.
			long double probability_models[10] = {0};
			for (int i = 0; i < 10; i++)
			{

				probability_models[i] = forward_procedure_test(Frames[digit][record], i);
			}
			long double max_prob = 0.0;
			int detect_digit = -1;
			fprintf(fp, "\t----------------------------------------------------------------------------\n");
			fprintf(fp, "\t---------------------- TEST FILE  %d ---------------------------------------\n", x);
			fprintf(fp, "\t----------------------------------------------------------------------------\n\n");
			for (int i = 0; i < 10; i++)
			{
				if (probability_models[i] > max_prob)
				{
					max_prob = probability_models[i];
					detect_digit = i;
				}
				fprintf(fp, "\n\t P(Observation Seq/ Lambda) for Model of Digit %d -------->  %g", i, probability_models[i]);
			}
			if (detect_digit == digit)
			{
				accuracy++;
				total_accuracy++;
			}
			fprintf(fp, "\n\n\n\t~~~~~~~~~~~~~~~~~~~>> Actual Digit = %d\n", digit);
			fprintf(fp, "\t~~~~~~~~~~~~~~~~~~~>> Detected Digit = %d\n\n", detect_digit);
			record = x;
		}
		printf("\n\n\tACCURACY OF DIGIT %d------>%d%%\n", digit, accuracy * 10);
		fprintf(fp, "\n\n ACCURACY OF DIGIT %d------>%d%%", digit, accuracy * 10);
		fclose(fp);
	}
	printf("\n\n\t-----------------------------------------------------------------------------------------------------------\n");
	printf("\tTOTAL  ACCURACY ----->%d%%\n", total_accuracy);
	printf("\t-----------------------------------------------------------------------------------------------------------\n");
	printf("\n\tPlease check the probabilities corresponding to each record file in the TEST PROBABILITY folder!!!!!!\n\n");
}

// Live testing on the files.
void testing_and_accuracy_live_recording()
{
	int i = 0, j = 0, k = 0;
	int total_accuracy = 0;

	int accuracy = 0;

	system("Recording_module.exe 3 test.wav word.txt");
	FILE *fp = fopen("word.txt", "r");
	if (!fp)
	{
		printf("error at 1612\n");
	}
	FILE *fp1 = fopen("word.txt", "r");
	if (!fp)
	{
		printf("error at 1617\n");
	}
	long double temp_data[50000] = {0};
	long double f_array[150][320] = {0};
	long double value = 0;
	long double max_value = 0.0;
	int index = 0;
	while (!feof(fp))
	{
		fscanf(fp, "%Lf", &value);
		temp_data[index] = value;
		index++;
		if (value > max_value)
		{
			max_value = value;
		}
	}
	fclose(fp);

	int start = -1;
	int end = -1;
	int num_frame = 0;
	for (int i = 0; i < 48000; i++)
	{
		if (temp_data[i] > 500)
		{
			start = i;
			break;
		}
	}
	int x = start;
	for (int i = 0; i < 100; i++)
	{
		for (int j = 0; j < 320; j++)
		{

			f_array[i][j] = temp_data[x];
			x++;
		}
		x = start + i * 80;
	}

	long double dc_shift = Compute_dc_shift();
	for (i = 0; i < 100; i++)
	{
		for (j = 0; j < S; j++)
		{
			f_array[i][j] = ((f_array[i][j] - dc_shift) * 5000.0) / max_value;
		}
	}

	double R_ref[13];
	double A_ref[13];
	double C_ref[13];
	double W_ref[13];
	long double reference[150][13] = {0};
	for (i = 0; i < 100; i++)
	{
		double temp[320] = {0};
		for (int j = 0; j < 320; j++)
		{
			temp[j] = f_array[i][j]; //store the 320 samples pf each frame in a sample array to calculate the Ri,Ci and Ai values
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

	// Compute the observation sequence for the test file.
	long double tokhura_weights[13] = {0.0, 1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
	long double distance[33] = {0};
	long double tokhura_distance = INT_MAX, final_distance = 0.0;
	int obs_index = -1;
	for (int frame = 0; frame < 100; frame++)
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

	//Comparision of the probabilities obtained from all models.
	long double probability_models[10] = {0};
	for (int i = 0; i < 10; i++)
	{
		probability_models[i] = forward_procedure_test(100, i);
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
		printf("\n\tP(Observation Seq/ Lambda)for Model %d------->  %g", i, probability_models[i]);
	}
	printf("\n\n\n\tDetected digit---->  %d\n\n\n", detect_digit);
}

int _tmain(int argc, _TCHAR *argv[])
{
	//create_universe();
	//create_codebook();

	printf("\n\n\n\t\t\t Darshika Verma\n\t\t\t 214101014\n\t\t\t Assignment-5\n\t\t\t DIGIT RECOGNITION \n\n\n\n");

	FILE *fp, *fp1, *fp2;
	long double values = 0.0;
	int temp = 0;

	// Read the universe
	fp = fopen("Final_Universe.txt", "r");
	if (!fp)
	{
		printf("File is not found at line  no ,1788");
	}
	for (int digit = 0; digit < 10; digit++)
	{
		for (int record = 0; record < 20; record++)
		{
			for (int frame = 0; frame < 85; frame++)
			{
				for (int index = 1; index <= 12; index++)
				{
					fscanf(fp, "%Lf", &universe[digit][record][frame][index]);
				}
			}
		}
	}
	fclose(fp);

	//Read the number of frames for each recorded file.
	fp2 = fopen("Final_Frames.txt", "r");
	if (!fp2)
	{
		printf("File is not found at line  no  1813");
	}
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 20; j++)
		{
			fscanf(fp2, "%d", &Frames[i][j]);
		}
	}
	fclose(fp2);

	//Read the codebook.
	fp1 = fopen("Final_Codebook.txt", "r");
	if (!fp1)
	{
		printf("File is not found at line  no ,1829");
	}
	for (int i = 0; i < 32; i++)
	{
		for (int j = 1; j <= 12; j++)
		{
			fscanf(fp1, "%Lf", &code_book[i][j]);
		}
	}

	fclose(fp1);

	printf("\tThe Models are trained.\n\n");

	training();

	printf("\n\n\tEnter 1 to perform testing on recorded files.\n\n");
	printf("\tEnter 2 to perform live testing on files.\n\n");
	printf("\tEnter the number: ");

	int x;
	scanf("%d", &x);
	if (x == 1)
	{
		testing_and_accuracy_recorded_files();
	}
	else if (x == 2)
	{
		int test;
		printf("\tEnter the number of times you want to perform live testing?  ");
		scanf("%d", &test);
		while (test--)
		{
			testing_and_accuracy_live_recording();
		}
	}
	else
	{
		printf("\tPlease enter the valid number\n\n");
	}

	system("PAUSE");
	return 0;
}