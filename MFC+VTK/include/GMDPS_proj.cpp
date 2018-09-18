//#include "stdafx.h"

//共用的一些函数的实现文件
#include "GMDPS_proj.h"

void Assign_Array2(double** data,int rows,int cols,int indata)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			data[i][j]=indata;
		}
	}
}

int Sign(double number)
{
	if (number>0)
	{
		return 1;
	}else if (number<0)
	{
		return -1;
	}else
	{
		return 0;
	}
}

double** CreateArray2(int rows,int cols)
{
	double** array0=new double* [rows];
	for (int i = 0; i < rows; i++)
	{
		array0[i]=new double[cols];
		array0[i][1]=2;
	}
	return array0;
}
void DeleteArray2(double** array2,int rows,int cols)
{
	for (int i = 0; i < rows; i++)
	{
		if (array2)
		{
			delete array2[i];
		}
	}
	if(array2) delete array2;
	//array2=NULL;
}

void GetGrdMinMax(double** data,CGM_GRDDATAINFO & grdinfo)
{
	double zmin=data[0][0],zmax=data[0][0];
	for (int i = 0; i < grdinfo.GetNumber_y(); i++)
	{
		for (int j = 0; j < grdinfo.GetNumber_x(); j++)
		{
			if (zmax<data[i][j])
			{
				zmax=data[i][j];
			}
			if (zmin>data[i][j])
			{
				zmin=data[i][j];
			}
		}
	}
	grdinfo.m_Ranges[0]=zmin;
	grdinfo.m_Ranges[1]=zmax;
}

//矩阵相乘
void Mat_Multiply(double* A,double* B,double* C,int rows,int cols)
{
	for (int i = 0; i < rows; i++)
	{
		C[i]=0;
		for (int j = 0; j < cols; j++)
		{
			C[i]+=A[j+i*rows]*B[j];
		}
	}
}

double Distance(double* point1,double *point2,int dims)
{
	int distance=0;
	for (int i = 0; i < dims; i++)
	{
		distance+=pow(point1[i]-point2[i],2.0);
	}
	return sqrt(distance);
}

double Length_Vector(double* vector0, const int dims)
{
	double length = 0;
	for (int i = 0; i < dims; i++)
	{
		length += vector0[i] * vector0[i];
	}
	return sqrt(length);
}

double GetMax(double* data,const int dims)
{
	double max=data[0];
	for (int i = 0; i < dims; i++)
	{
		if (max<data[i])
		{
			max=data[i];
		}
	}
	return max;
}
double GetMin(double* data,const int dims)
{
	double min=data[0];
	for (int i = 0; i < dims; i++)
	{
		if (min>data[i])
		{
			min=data[i];
		}
	}
	return min;
}

double* Cross(double* a, double* b)
{
	double result[3];
	result[0] = a[1] * b[2] - a[2] * b[1];
	result[1] = a[2] * b[0] - a[0] * b[2];
	result[2] = a[0] * b[1] - a[1] * b[0];

	return result;
}

double VectorDot(double *a, double *b, int dims)
{
	double result = 0;
	for (int i = 0; i < dims; i++)
	{
		result += a[i] * b[i];
	}
	return result;
}