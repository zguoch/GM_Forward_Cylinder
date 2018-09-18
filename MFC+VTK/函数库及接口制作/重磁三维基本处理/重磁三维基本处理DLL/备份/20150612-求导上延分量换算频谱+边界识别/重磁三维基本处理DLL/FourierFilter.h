#include "GMDPS_proj.h"
#include <cmath>
#include "Ext.h"
#include "Filter.h"
#include "文件操作.h"
#include "FourierFilter_1D.h"
#include "高斯消元解方程组.h"
//#include "ContinuationInSpace.h"
#include "GMQR.h"
#include "GMIV.h"
//网格信息
struct MESHINFO
{
	int M, N, L;		//x,y,z方向上的网格个数
	double dx, dy, dz;	//x,y,z方向上的网格宽度
	double initx, inity, initz;	//起始坐标
};

struct MODELINFO					//模型剖分信息
{
	double initx, inity, initz;		//起始坐标
	double dx, dy, dz;				//模型x,y,z方向的间距
	int M, N, L;						//X,Y,Z方向上的模型个数
	double z0;						//观测面高度，深度加权函数中的z0
};

bool UpWardContinuation(double** OriginalData,double** TransData,const int rows,const int columns,const double dx,const double dy,const double rph);
bool UpWardContinuation(double* OriginalData,double* TransData,const int rows,const int columns,const double dx,const double dy,const double rph);
bool ComponentTransProcess(double** OriginalData,double** TransData,const int rows,const int columns,const double dx,
						   const double dy,const double angle_I,const double angle_D,const int OriginalStyle,const int TransStyle);
bool DownWardContinuation_Tik(double** const OriginalData,double** TransData,const int rows,const int columns,const double dx,const double dy,const double rph,const double TRP,double XiGeMe);
double AutoDownWardContinuation_Tik(double** OriginalData,double** TransData,const int rows,const int columns,
									const double dx,const double dy,const double rph,const double TRP,
									const double dTRP,const int kmax,const double XiGeMe);
bool DownWardContinuation_LandWeber(double** const OriginalData,double** TransData,const int rows,const int columns,const double dx,const double dy,const double rph,const double a_Land,const double n_Land);
bool DownWardContinuation_Compensation(double** const OriginalData,double** TransData,const int rows,const int columns,const double dx,const double dy,const double rph,const int kmax,const double TRP,const double XiGeMe);
bool PowerSpectrum(double** OriginalData,double* PowerSpectrumAmp,int mr,const int rows,const int columns,const double dx,const double dy);
bool InitFFT(GridDataInfo datainfo,int& m1,int& me,int& n1,int& ne);
bool SolveEquation(double** A,double* X,double* B,int rows,int columns);
double SelectFrequenceSegment(double* PSA,double DX,int mr,double df,double fstart,double fend,char* path);
double EquivalentLayerModel(double* PSA,int mr,double* S,double* depth,double layernum,double df,char* path);
bool HorizontalGradient(double** const OriginalData,double** TransData,const int rows,const int columns,const double dx,const double dy,const double Derivangle);
bool VerticalGradient(double** const OriginalData,double** TransData,const int rows,const int columns,const double dx,const double dy);
bool TitlAngle(double** const OriginalData,double** TransData,const int rows,const int columns,const double dx,const double dy);