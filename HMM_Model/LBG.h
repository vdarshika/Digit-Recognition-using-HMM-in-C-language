#ifndef lbg_h
#define lbg_h
// 214101014_LBG.cpp : Defines the entry point for the console application.
//
/*
									SPEECH PROCESSING ASSIGNMENT-5
										 Digit Recognition Algorithm

									NAME: DARSHIKA VERMA
									ROLL_NO:214101014 

									*/

//include header files that are required in the program.
#include "stdafx.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"

const int total_vectors=18859;					 //represents the total data vectors present in the universe.
const int total_ci=12;							//represents the total number of cepstral coefficients.
const int clusters=32;						    //number of clusters into which the universe is to be divided.

long double data_vectors[18859][12];					// 2D matrix to store the values of data vectors present in the universe. 
long double density_of_clusters[32][18859]={0};		// 2D matrix to store the indexes of the data vectors present in each cluster.
long double code_book[32][12];						// 2D matrix to store the codebook.
long double temp_code_book[32][12]={0.0};			// 2D matrix to store the temporary codebook values inorder to split the codebook.

long double initial_distortion=9999.0;
long double new_distortion=0.0;
int final_code_book_size=1;
long double epsilon=0.03;

// This functions takes a data vector as a input and find the distance of that data vector from each cluster of codebook
// and return the index of that cluster which has the minimum distance to that data vector.
int compute_distance(long double data[])
{
	int i=0,j=0,index;
	long double dist=999;
	long double sum=0.0;

	//given tokhura weights 
	long double tokhura_weights[12]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

	// calculating distance for each cluster.
	for(i=0;i<clusters;i++)
	{
		sum=0.0;
		for(j=0;j<total_ci;j++)
		{
			sum=sum+(tokhura_weights[j]*((data[j]-code_book[i][j])*(data[j]-code_book[i][j])));
		}
		
		sum=sum/12.0;
		
		//finding the index with the minimum distance with the data vector.
		if(sum<dist)
		{
			dist=sum;
			index=i;
		}
	}

	return index;
}

// This function takes a row of the density_of_cluster array and then find the total number of data vectors stored in that cluster.
// and return the total number of data vectors stored in the chosen cluster.
long int compute_existing_number_of_vectors(long double vectors[])
{
	long int count=0;
	for(int i=0;i<18859;i++)
	{
		// 0 represents that there is no further index of any data vector present in the chosen cluster.
		if(vectors[i]==0)
			break;

		count++;
	}
	return count;
}

//This function takes all the data vectors present in the universe and calculates the minimum distance using tokhura distance with the cluster
// and then calculate the distortion for all the data vectors and then calculate the average of the distortion.
long double compute_distortion()
{
	int i=0,j=0,k=0;

	long double distortion=0.0;
	long double min_distance=9999.0;
	long double distance=0.0;
	long int count=0;

	long double tokhura_weights[12]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
	
	for(i=0;i<total_vectors;i++)
	{
		min_distance=999.0;
		
		for(j=0;j<8;j++)
		{
			distance=0.0;
			count=0;
			// calculate the distance for each cluster.
			for(k=0;k<12;k++)
			{
				distance+=((data_vectors[i][k]-code_book[j][k])*(data_vectors[i][k]-code_book[j][k])*tokhura_weights[k]);
			}
			
			//calculating the minimum distance.
			if(distance<min_distance)
			{
				min_distance=distance;
			}
		}
		distortion=distortion+min_distance;
	}
	//calculate the average of distortion.
	distortion=distortion/18859.0;

	return distortion;
}

// This function re-initialise the values of the density_of_clusters.
// So that whenever new data vectors are assigned to the cluster the previous clusters are removed.
void empty_density_of_clusters()
{
	int i=0,j=0;
	for(i=0;i<clusters;i++)
	{
		for(j=0;j<total_vectors;j++)
		{
			density_of_clusters[i][j]=0; //removing the assigned data vectors to the clusters.
		}
	}
}

// This function updates the centroids value of each cluster present in the codebook.
void updation_of_centroids()
{
	int i=0,j=0,k=0,vector_index=0,count=0;

	long double centroid[12];		// array that stores the calculated centroid value.
	for(i=0;i<clusters;i++)
	{
		memset(centroid, 0, sizeof(centroid));		//initialises the centroid array with 0 values.
		count=0;	//represents the total number of data vectors present in the cluster.

		for(j=0;j<total_vectors;j++)
		{
			vector_index=density_of_clusters[i][j];
			if(vector_index!=0)
			{
				for(k=0;k<12;k++)
				{
					centroid[k]=centroid[k]+data_vectors[vector_index][k];
				}
				count++;
			}
			else
			{
				break;
			}
		}

		//updating the value of centroid in the codebook.
		for(k=0;k<12;k++)
		{
			// average out the value of centroid calculated by total number of data vectors present in the cluster.
			code_book[i][k]=(centroid[k]/(double)count);
		}
	}
}

// This function classify all the data vectors present in the universe into their respective clusters based on the
//  minimum distance of the cluster to the data vector.
void classify_data_vectors()
{
	int i=0,j=0;

	for(i=0;i<total_vectors;i++)
		{
			long int cluster=0;
			long double data[12];
			
			// convert the 2d array into 1d array so as to calculate the distance of that vector from each cluster.
			for(j=0;j<total_ci;j++)
				{
					data[j]=data_vectors[i][j];
				}
			
			// finding the cluster where the data vector is to be assigned.
			cluster=compute_distance(data);

			// finding the index where the index of the data vector is to be stored in the density of cluster array.
			int entry_index= compute_existing_number_of_vectors(density_of_clusters[cluster]);

			// assigning the index of the data vector in the respective cluster.
			density_of_clusters[cluster][entry_index]=i;
		}
}

// This function performs the kmeans algorithm and call all the required fucntion for finding the optimal codebook.
void k_means_algorithm()
{
	int iteration=0;
	
	// assign all data vectors to their respective clusters.
	classify_data_vectors();

	// calculate the initial distortion of the codebook.
	new_distortion=compute_distortion();
	
	do
	{
		//update the values of the codebook with the new values to achieve the optimality in the code book.
		updation_of_centroids();

		// clear the assigned values from the clusters.
		empty_density_of_clusters();

		// assign all data vectors to their respective clusters.
		classify_data_vectors();

		initial_distortion=new_distortion;
		iteration++;
		new_distortion=compute_distortion();
		//printf("\n %Lf\n",new_distortion);
	}
	while((initial_distortion-new_distortion)>0.0001);

}

// This function splits the codebook and doubles the size of the codebook.
void split_code_book()
{
	int i=0,j=0,k=0;
	
	//initialise the temp_code_book with 0.0 values.
	memset(temp_code_book, 0.0, sizeof(temp_code_book));

	// these temp arrays are to store the values of the centroid that is split.
	long double temp_split_array_1[12]={0.0};
	long double temp_split_array_2[12]={0.0};
		
	//store the value of codebook into the temporary codebook.
	for(i=0;i<final_code_book_size;i++)
	{
		for(j=0;j<total_ci;j++)
		{
			temp_code_book[i][j]=code_book[i][j];
		}
	}

	//this shows where to store the newly generated centroid in the codebook.
	int index=0;

	for(i=0;i<final_code_book_size;i++)
	{
		//split the codebook into two parts.
		for(j=0;j<total_ci;j++)
		{
			temp_split_array_1[j]=code_book[i][j]*(1+epsilon);
			temp_split_array_2[j]=code_book[i][j]*(1-epsilon);
		}

		//storing the values into the codebook.
		for(k=0;k<total_ci;k++)
		{
			code_book[index][k]=temp_split_array_1[k];
		}

		index++;

		//storing the values into the codebook.
		for(k=0;k<total_ci;k++)
		{
			code_book[index][k]=temp_split_array_2[k];
		}
		index++;
	}
}


//This function performs the LBG algorithm and nake the function call to all the required functions.
void LBG_algorithm()
{
	printf("\n-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n");
	
	//loop continues until we reach the codebook of required size.
	while(final_code_book_size<=32)
	{
		printf("   Distortion for Codebook of Size %d -->> %Lf\n\n",final_code_book_size,new_distortion);

		if(final_code_book_size>=32)
		{
			break;
		}

		//function call to the split fucntion until the required size codebook is obtained.
		split_code_book();

		//function call to k means algorithm to find the optimal codebook after performing the split operation.
		k_means_algorithm();

		final_code_book_size=2*final_code_book_size;
	}
	printf("\n-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n");
	
	return;
}



#endif