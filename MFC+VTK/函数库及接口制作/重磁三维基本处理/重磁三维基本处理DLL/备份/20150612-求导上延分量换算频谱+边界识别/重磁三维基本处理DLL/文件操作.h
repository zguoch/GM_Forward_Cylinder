
//#include "stdafx.h"
#include "vector"
using namespace std;

#ifndef FILE2
#define FILE2
typedef struct tagGridDataInfo
{
	char*	Description;							//数据的描述
	double	gridXmin, gridXmax, gridYmin, gridYmax;	//测区范围
	double	gridZmin, gridZmax;						//异常数值最小和最大
	double	deltX, deltY, ExGridX, ExGridY;			//测区点、线距, 扩充后的大小
	int		rows, columns, ExRows, ExColumns;		//数据点线数
	double	DirectProfile, DirectBaseline;			//测区测线、基线方位角
	double	dTI0, dTD0;								//地磁场方位角（倾角、偏角）
}GridDataInfo;


vector<double> ReadGrd(char* infile,GridDataInfo& datainfo);
void WriteGrd(char *filename,double* z,const int dataNumber,int m,int n,double xmin,double xmax,double ymin,double ymax);
void ReadDat(char* filename,vector<double>& x,vector<double>& y,double& dx);
void WriteDat(char* filename,vector<double> x,double* y);
void WriteDat(char* filename,const int datanum,double* x,double* y);
void WriteDat(char* filename,const int datanum,double xmin,double dx,double* y);
void WriteDat(char* filename,vector<double>x,vector<double> y);
void ShinkGrd(char* filename,const int ShinkXnum1,const int ShinkXnum2,const int ShinkYnum1,const int ShinkYnum2);

#endif