#include <math.h>
#include <stdlib.h>
#include <fftw3.h>
#define N 16
int D=1;
int* dimensions;
int L=1000;
const double pi=M_PI;
const double x_min=0.0;
const double x_max=1.0;
const double y_min=0.0;
const double y_max=1.0;

double f(double x,double y)
{
	
	return 8.0*pi*pi*cos(4.0*pi*y)*(cos(4.0*pi*x) - sin(4.0*pi*x)) -16.0*pi*pi*(sin(4.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*y) + 
                       sin(2.0*pi*x)*sin(2.0*pi*x)* cos(4.0*pi*y));
}


void freq(int n,float *arr)
{
	// float *arr;
	// arr=malloc(sizeof(float)*n);
	
		/* code */
		int i=0;
		while(i<n/2)
		{
			arr[i]=(float)i/n;
			// printf("%f\n",arr[i] );
			i+=1;
		}
		while(i<n)
		{
			arr[i]=(float)(i-n)/n;
			i+=1;
		}
		// return arr;
	
}
void Creategrid(double **  grid,double ** f_grid)
{
	
	//boundary
	//x=0, f=0 x=1 f=0, y=0 == y=1
	for (int i = 0; i < N+2; ++i)
	{
		grid[i][N+1]=0;
		grid[0][i]=0;
		grid[N+1][i]=0;
		grid[i][0]=0;
	}
	
	double dx=(x_max-x_min)/N;
	double dy=(y_max-y_min)/N;
	for (int i = 0; i < N; ++i)
	{
		/* code */
		for (int j = 0; j< N; ++j)
		{
			/* code */

			f_grid[i][j]=f((i+0.5)*dx,(j+0.5)*dy);
			// printf("%lf %lf\n",(i+0.5)*dx,(j+0.5)*dy);
		}
	}



}

// void freq(int *kx)
// {
// 	if (N%2==0)
// 	{
// 		/* code */
// 		for (int i = 0; i < N; ++i)
// 		{
// 			/* code */

// 		}
// 	}
// }

void fft(double ** grid, double** f_grid)
{

	fftwl_complex *in,*out;
	fftwl_plan my_plan;
	in=(fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex)*N*N);
	out=(fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex)*N*N);
	int k=0;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			// printf("%lf %s\n",f_grid[i][j] ,"f_grid");
				in[k][0]=f_grid[i][j];
				in[k++][1]=0;

		}
	}	
	// int kx[N]=freq(N);
	double dx=(x_max-x_min)/N;
	double dy=(y_max-y_min)/N;
	// int ky[N]=freq(N);
	float *kx,*ky;
	kx=malloc(sizeof(float)*N);
	ky=malloc(sizeof(float)*N);
	freq(N,kx);
	freq(N,ky);
	
	my_plan=fftwl_plan_dft_2d(N,N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	fftwl_execute(my_plan);
	k=0;
	long double a,b;

	a=out[0][0];
	b=out[0][1];
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			float ki=kx[i]/dx;
			float kj=ky[j]/dy;	
			// printf("%f %f\n",ki,kj );
			// printf("%f %f %d %d####\n",ki,kj,i,j );
			in[k][0]=out[k][0]/(2*((cos((2*pi*ki)/N))-1)/(dx*dx) + (cos((2*pi*kj)/N)-1)/(dy*dy) );
			in[k][1]=out[k][1]/(2*((cos((2*pi*ki)/N))-1)/(dx*dx) + (cos((2*pi*kj)/N)-1)/(dy*dy) );
			// printf("%Lf %Lf %lf dd\n",in[k][0],in[k][1], ((cos((2*pi*ki)/N))-1)/(dx*dx) + (cos((2*pi*kj)/N)-1)/(dy*dy)   );
			k++;

		}
	
	}	
	in[0][0]=a;
	in[0][1]=b;
	my_plan=fftwl_plan_dft_2d(N,N,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
	fftwl_execute(my_plan);
	// for (int i = 0; i < N*N; ++i)
	// {
	// 		printf("%Lf %Lf %s\n",out[i][0],out[i][1],"output");
	// // 		// printf("%Lf %Lf %s\n",in[i][0],in[i][1],"input");
	// // 			// in[k][0]=f_grid[i][j];
	// // 			// in[k++][1]=0;

		
	// }	
	k=0;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N-1; ++j)
		{
			/* code */
			printf("%Lf, ",out[k][0]);
			k+=1;

		}
		printf("%Lf ",out[k][0]);
		k+=1;
		printf("\n");
	}
	fftwl_destroy_plan(my_plan);
	fftwl_free(in);
	fftwl_free(out);

}
int main(int argc, char const *argv[])
{
	
	
	
	// int dy=L/dimensions[1];
		// Create a grid
	double n;
	// n=f(0.3,0.3);
	// printf("%lf\n", n);
	double ** grid;
	double **f_grid;
	grid=calloc(N+2,sizeof(double*));
	for (int i = 0; i < N+2; ++i)
	{
		grid[i]=calloc(N+2,sizeof(double));
	}
	f_grid=calloc(N,sizeof(double));
	for (int i = 0; i < N; ++i)
	{
		f_grid[i]=calloc(N,sizeof(double));
	}
	// grid= malloc(size);
	Creategrid(grid,f_grid);
	// for (int i = 0; i < N; ++i)
	// {
	// 	for (int j = 0; j < N; ++j)
	// 	{
	// 		printf("%lf\n",f_grid[i][j] );
	// 			// in[k][0]=f_grid[i][j];
	// 			// in[k++][1]=0;

	// 	}
	// }
	fft(grid,f_grid);

	return 0;
}