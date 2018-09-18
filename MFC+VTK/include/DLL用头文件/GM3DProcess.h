//其实现在：函数库及接口制作\重磁三维正演\重磁三维正演DLL 目录下
//这是个导出的dll用的头文件

//重磁三维数据处理

#include "GMDPS_proj.h"
#ifndef GM3DPROCESS
#define GM3DPROCESS

#define GMDLL extern "C" _declspec(dllexport)
namespace GM3DProcess
{
	/*-----------------------//
	求导数
	//------------------------*/
	/*输入CGM_GRDDATA包含了网格信息，输出处理后的数据*/
	/*Orientation表示方向余弦分别为x方向，y方向，z方向*/
	GMDLL int DerivativeX(CGM_GRDDATA indata, CGM_GRDDATA* processdata);
	GMDLL int DerivativeY(CGM_GRDDATA indata, CGM_GRDDATA* processdata);
	GMDLL int DerivativeHorizontal(CGM_GRDDATA indata, double* Orientation, CGM_GRDDATA* processdata);		//Orientation元素分别为cosx，cosy
	GMDLL int DerivativeHorizontal2(CGM_GRDDATA indata, double aerfa_x, CGM_GRDDATA* processdata);			//aerfa_x表示求导方向与x正方向的夹角(顺时针):度
	GMDLL int DerivativeVertical(CGM_GRDDATA indata, CGM_GRDDATA* processdata);
	GMDLL int Derivative(CGM_GRDDATA indata, double* Orientation, CGM_GRDDATA* processdata);				//Orientation元素分别为cosx，cosy，cosz
	GMDLL int Derivative2(CGM_GRDDATA indata, double* Angle, CGM_GRDDATA* processdata);						//Angle元素分布为偏角和倾角（度）
	
	/*-----------------------//
	频谱
	//------------------------*/
	GMDLL int AmplitudeSpectrum(CGM_GRDDATA indata, CGM_GRDDATA* processdata);
	GMDLL double* MeanRadialPowerSpectrum(CGM_GRDDATA indata, int& frequenceNum);							//平均径向功率谱
	GMDLL int DeleteArray1D(double* array1d);

	/*-----------------------//
	延拓
	//------------------------*/
	GMDLL int UpwardContinuation(CGM_GRDDATA indata, double height, CGM_GRDDATA* processdata);				//height:延拓高度（m）

	/*-----------------------//
	磁异常转换
	//------------------------*/
	GMDLL int MagComponentTransT2X(CGM_GRDDATA indata, MagneticComponentTransStruct*, CGM_GRDDATA* processdata);//这里有个xy坐标系的问题，与程序里的xy对应起来
	GMDLL int MagComponentTransT2X_ext(CGM_GRDDATA indata, MagneticComponentTransStruct* magdatainfo, CGM_GRDDATA* processdata);//先扩边在计算再缩边，貌似有问题，以后再修改吧
	GMDLL int MagComponentTransT2Y(CGM_GRDDATA indata, MagneticComponentTransStruct*, CGM_GRDDATA* processdata);//
	GMDLL int MagComponentTransT2Z(CGM_GRDDATA indata, MagneticComponentTransStruct*, CGM_GRDDATA* processdata);//
	//模量
	GMDLL int MagComponentTransT2Ta(CGM_GRDDATA indata, MagneticComponentTransStruct* magdatainfo, CGM_GRDDATA* processdata);
	GMDLL int MagComponentTransT2R(CGM_GRDDATA indata, MagneticComponentTransStruct* magdatainfo, CGM_GRDDATA* processdata);
	GMDLL int MagComponentTransT2E(CGM_GRDDATA indata, MagneticComponentTransStruct* magdatainfo, CGM_GRDDATA* processdata);
	GMDLL int MagComponentTransT2Q(CGM_GRDDATA indata, MagneticComponentTransStruct* magdatainfo, CGM_GRDDATA* processdata);
	GMDLL int MagComponentTransT2L(CGM_GRDDATA indata, MagneticComponentTransStruct* magdatainfo, CGM_GRDDATA* processdata);

	/*-----------------------//
	边界识别
	//------------------------*/
	GMDLL int EdgeDetective_ImprovedTiltDerivative(CGM_GRDDATA indata, CGM_GRDDATA* processdata);//改进的Tilt导数

	/*-----------------------//
	扩边
	//------------------------*/
	GMDLL int ExtenBoundary(CGM_GRDDATA indata, CGM_GRDDATA* processdata, int method = EXTENBOUNDARY_COS);//processdata在函数外面进行初始化，也就是说已经知道要扩边多少了，此函数只是填值而已
}


#endif