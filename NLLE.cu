#include <windows.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cuda.h>
#include <time.h>
#include <string.h>

int NPoints = 90000 ;
int LDimension = 3;
int Mcolumn = 90;
int HMcolumn = Mcolumn;

int TCycle = (int)(1e3);
int PTCycle = (int)(1e2);
float Rc1 = 1.0;
float LLL = (float)(TCycle * 1.0);

int Region = 2 ;
char fname_cx[100] =  "C:\\Users\\Safanajim\\Dropbox\\mFSPEx";//3.txt"  ;
//--------------------------------------------------------------------
__host__ void storeX(float *x, int N, int LDimension)
{/*  	
	FILE *fp1; 
 fp1 = fopen("H:\\Reports\\Step_by_Step_MDS\\MDS4.txt", "w");	 
	for (int i=0; i<N*LDimension; i++)
       fprintf(fp1,"%f,",x[i]);        
  fclose(fp1);
 */
 FILE *fp1;   	  
	char fname[200]="";
	strcat(fname, fname_cx);	
	char Region_str[3] = "";
	itoa(Region, Region_str, 10);
	strcat(fname, Region_str);
	strcat(fname, ".txt");
	//printf("xxxx) %s\n", fname); getchar();
 fp1 = fopen(fname, "w");	 
 int j=0;
	for (int i=0; i<N; i++)
			for (int k=0; k<LDimension; k++)
				{fprintf(fp1,"%f,",x[j]);  				
				 //printf("%d) %f\n", j, x[j]);
				 j++;				 
				}
 fclose(fp1); 
 
}
//---------------------------------------------------------------------
//----------------------------------------------------------------------------------
__host__ void SPE_CPU(int a, float *OData, float lambda,float *x, float rcut, int Dimension, int LDimension, int N)
{   //int b = blockIdx.x * blockDim.x + threadIdx.x;
				//if (b < N)   
				for(int b=0; b<N; b++)
				{			float sum = 0.0;
					   for(int k=0; k<LDimension; k++)
					      sum += (float)((x[a+ N*k] - x[b+ N*k]) * (x[a+N*k] - x[b+N*k]));												  
					   float dab = sqrt(sum);
					   float rab = dab;
			    															              															
					   sum = 0.0;															
					   for(int k=0; k<Dimension; k++)
								   sum += (float)((OData[a + N*k] - OData[b+ N*k]) * (OData[a+ N*k] - OData[b+ N*k]));	
					   rab = sqrt(sum);
							
			 	   if (((rab <= rcut) || ((rab>rcut) && (dab < rab))))
			     //if ((dab<=rcut) || ((dab > rcut)&&(dab<rab)))																	
								   {   float T;
												   T = (float)(lambda * (rab - dab) / (dab + 1e-8));																																		          								
												   for (int k=0; k<LDimension; k++)
                   x[b+N*k] += (float)(T * (x[b+N*k] - x[a+N*k]));    																																																																													
									   }													  
       } 
       //__syncthreads();
  } 
//----------------------------------------------------------------------------------
__host__ void initial_X(float *x, int N, float rcut, int LDimension)
{   int q = RAND_MAX;
    for (int i = 0; i < N*LDimension; i++) 
       x[i] = (float)(1.0 * rand() / q); 
}
//--------------------------------------------------------------------------------
__host__ void stress(float *x, float *OData, int N, float &S, float rcut, int Dimension, int LDimension)
{		int i, a, b;		
		float b1=0, b2=0;
		for (i=0; i<N; i++)
		//for (j=i+1; j<N; j++)
		{ 
			a = (int)(rand() * (float)(N-1) / (RAND_MAX+1.0));
			while(1) {
					b = (int)(rand() * (float)(N-1) / (RAND_MAX+1.0));
					if (b == a) continue;
					else break;
		}
		//a = i;
		//b = j;      
		float sum = 0.0;
        for(int k=0; k<LDimension; k++)
			sum += (float)((x[a+N*k] - x[b+N*k]) * (x[a+N*k] - x[b+N*k]));
        float dab = sqrt(sum);
								
		int k;
		sum = 0.0;
		for(k=0; k<Dimension; k++)
								sum += (float)((OData[a + N*k] - OData[b+ N*k]) * (OData[a+ N*k] - OData[b+ N*k]));	
		float rab = sqrt(sum); 
		//(abs(rab) > 0)&&
		if (((rab <= rcut) || ((rab>rcut)&&(dab < rab))))
		{  b1 += (float)((dab - rab) * (dab - rab)/(1e-8 + rab));
		   b2 += rab ;}
		}
		S = (float) (b1/(1e-8+b2));
}
//---------------------------------------------------------------------------------
__host__ void get_MD(float *OData, int N, int Dimension)
{
  int i, j;	
	float x1[90000][100]; 
 int ik;
	for (ik=0; ik<Dimension; ik++)
		{  char fname[100] =  "H:\\head\\Moffet_test"; //3\\Gulf";		
		   char Region_str[3] = "";
	    itoa(Region, Region_str, 10);
		   strcat(fname, Region_str);
		   strcat(fname, "\\Gulf");
		   
		   //char fname[100] =  "H:\\head\\Moffet_test0\\Gulf";
		
		
     int fno2= ik+1;   			
		   char buffer[3]="";
		   itoa(fno2, buffer, 10);
		   strcat(buffer,".txt");
		   strcat(fname,buffer);
		   if (ik==0)
			     printf("%3d) %s\n", ik, fname);
		   FILE *fp1;
		   float Value;
		   fp1 = fopen(fname,"r");
		   for(i=0; i<N; i++)
			   {  fscanf(fp1,"%f,", &Value);
			  	   x1[i][ik]= Value;
			   }
		   fclose(fp1);
			}
   j=0;
   for (int ik=0; ik<Dimension; ik++)
		 for (int i=0; i<N; i++)
			 { OData[j] = x1[i][ik];
					 j++;
			 }	
		 //printf("Reading originl space is done\n");			
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                      
int main() 
{   
				int Dimension = Mcolumn; 
    int i;
    int N = NPoints;
    int k = 1;
    int ncycle = k * TCycle;
    int PrintN = PTCycle;      
    float rcut= Rc1 , lambda0 = 1.0;	   
    	
    srand ( (unsigned int)time(NULL) );
    
				float *OData, *dOData;
    OData = (float*) malloc(Dimension*N*sizeof(float));
    cudaMalloc(&dOData, Dimension*N* sizeof(float));
    
    get_MD(OData, N, Dimension);
    cudaMemcpy(dOData, OData, Dimension*N*sizeof(float), cudaMemcpyHostToDevice);   	
    printf("---- read Md Matrix  was done ----\n");		 
				//########################################################################################
				
    float *x, *xd;
				x=(float*)malloc(N*LDimension*sizeof(float)); 
				cudaMalloc(&xd,N*LDimension*sizeof(float));	   
									
				initial_X(x, N, rcut, LDimension);    
				cudaMemcpy(xd,x, LDimension*N*sizeof(float),cudaMemcpyHostToDevice);				

				int ThreadsPerBlock = 512;
				int BlocksNeeded = (N+ThreadsPerBlock -1)/ThreadsPerBlock ;
				dim3 dimGrid( BlocksNeeded );
				dim3 dimBlock( ThreadsPerBlock );
				/*
    printf("---------------- FSPE Starts----------------\n");
    clock_t Start0 = clock();
    //diff=(clock()-start)/(double) CLOCKS_PER_SEC
    printf("Starting time is %f\n", (float)Start0/CLOCKS_PER_SEC);
    */
				
				float S1 = LLL; 			
				float Rc_LL;				
				for(i=0; i<ncycle; i++)
					{ //int a = i;					
							int a = rand()% N;		
							float lambda = lambda0 - lambda0 * i / ncycle; 			
							//float lambda = lambda0/ (1.0+i); 	
							float LL = LLL/(1+i);			
							//float LL = LLL - LLL * i / ncycle; 				            						
							Rc_LL = LL ;
						 FSPE_GPU<<<dimGrid, dimBlock>>>(a, dOData, lambda, xd, Rc_LL, Dimension, LDimension, N);  	
						 //SPE_CPU(a, OData, lambda, x, Rc_LL, Dimension, LDimension, N);  							 						 
							if ((i% PrintN ==0) && (i !=0))
								{  cudaMemcpy(x, xd, LDimension*N*sizeof(float), cudaMemcpyDeviceToHost);	
								   S1 = 0.0;           
											stress(x, OData, N, S1, Rc_LL, Dimension, LDimension);
											printf("%d stress=%f\n", i, S1);
								}
			}	 
		/*	
		clock_t Start2 = clock();		
  printf("Ending  time is %f\n", (float)Start2/CLOCKS_PER_SEC);
  float diff = ((float)Start2-(float)Start0);///1000000;
  float seconds = diff /  CLOCKS_PER_SEC;
  printf("------------- Diff. time is %f\n", seconds);
  */
  
		cudaMemcpy(x, xd, LDimension*N*sizeof(float), cudaMemcpyDeviceToHost);	   
		storeX(x, N, LDimension);
		printf("----------------  Program End  -----------------\n");
		free(OData);
		free(x);		
		cudaFree(dOData);
		cudaFree(xd);
		
  //getchar();	  
		return 0;
}