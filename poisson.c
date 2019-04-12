#include <math.h>
#include <stdlib>

int D=2;
int* dimensions;
int L=1000;
double f(int x,int y)
{
	return x*x+y*y;
}
void Creategrid(double **  grid,int dx,int dy)
{
	for (int i = 0; i < dimensions[0]; ++i)
	{
		for (int j = 0; j < dimensions[1]; ++j)
		{
			grid[i][j]=feval(dx*i,dy*j);
		}
	}
}

int main(int argc, char const *argv[])
{
	
	dimensions = malloc(sizeof(int)*D);
	size_t size=sizeof(double);
	for (int i = 0; i < D; ++i)
	{
		/* code */
		dimensions[i]=100;
		size*=dimensions[i];
	}
	int dx=L/dimensions[0];
	int dy=L/dimensions[1];
		// Create a grid
	double ** grid;
	grid= malloc(size);
	Creategrid(grid,dx,dy);
	return 0;
}