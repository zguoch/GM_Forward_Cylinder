
/////Proj.h



#include <cmath>
#include "Ext.h"
#include "Filter.h"
#include "文件操作.h"
#include "FourierFilter.h"
#include "FourierFilter_1D.h"
#include "高斯消元解方程组.h"
//#include "ContinuationInSpace.h"
#include "GMQR.h"
#include "GMIV.h"

#ifndef proj
#define proj

//#define PI 3.1415926
#define GRAV_CONSTANT 6.667		//万有引力常数


//网格信息
//struct MESHINFO
//{
//	int M,N,L;		//x,y,z方向上的网格个数
//	double dx,dy,dz;	//x,y,z方向上的网格宽度
//	double initx,inity,initz;	//起始坐标
//};

//struct MODELINFO					//模型剖分信息
//{
//	double initx,inity,initz;		//起始坐标
//	double dx,dy,dz;				//模型x,y,z方向的间距
//	int M,N,L;						//X,Y,Z方向上的模型个数
//	double z0;						//观测面高度，深度加权函数中的z0
//};
double arctg(double v);
//{
//	double ang;
//	if(v>=0)
//	{
//		ang=atan(v);
//	}
//	else if(v<0)
//	{
//		ang=atan(v)+3.1415926535897932384626433832795;
//	}
//	return ang;
//}
//////////////////////////////申明GPUKernal.cu中的函数//////////////////////////////////////////////////
// declaration, forward
extern "C" 
	void GetGPUProperty();
extern "C"
	void GPU_GetAdk(int numElements,double* h_A,double* h_B,double* h_C,GridDataInfo datainfo);

#endif
////////////////////////////////////////////////////////////////////////////////