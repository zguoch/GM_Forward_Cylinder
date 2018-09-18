//其实现在：函数库及接口制作\重磁三维正演\重磁三维正演DLL 目录下
//这是个导出的dll用的头文件
//规则几何形体的重磁三维正演

#include "GMDPS_proj.h"
#include "fftw3.h"
#ifndef REGULAR3DFORWARD
#define REGULAR3DFORWARD

#define GMDLL extern "C" _declspec(dllexport)


//定义不同类型模型的标志ID
#define MODEL_global 1
#define MODEL_cylinder 2
#define MODEL_cube 3

//定义不同的正演类型
#define FORWARD_V 4
#define FORWARD_Vz 5
#define FORWARD_Vx 6
#define FORWARD_Vy 7
#define FORWARD_Vxx 8
#define FORWARD_Vxy 9
#define FORWARD_Vxz 10
#define FORWARD_Vyy 11
#define FORWARD_Vyz 12
#define FORWARD_Vzz 13
#define FORWARD_ERROR_Vxx	19
#define FORWARD_ERROR_Vxy	20
#define FORWARD_ERROR_Vxz	21
#define FORWARD_ERROR_Vyy	22
#define FORWARD_ERROR_Vyz	23
#define FORWARD_ERROR_Vzz	24
#define FORWARD_ERROR_Hax	25

#define FORWARD_Ta 14
#define FORWARD_Hax 15
#define FORWARD_Hay 16
#define FORWARD_Za 17
#define FORWARD_Module 18

//圆柱体模型结构体
struct GM_CylinderInfo
{
	double Pos1[3];								//一个端点的坐标（x，y，z）
	double Pos2[3];								//另一个端点的坐标(x,y,z)，pos1和pos2完全可以描述圆柱体的位置
	double Radius;								//半径
	double Density;								//密度：重力正演时使用,单位（k/cm3）
	double CiHuaLv,Angle_I,Angle_D;				//磁化率,磁化倾角，磁化偏角：磁法正演时使用（单位：度）
	double Resolution;							//用于vtk显示时的一个变量
	COLORREF Color;								//用于VTK显示是的变量
};
//立方体模型结构体
struct GM_CubeInfo
{
	double bounds[6];							//顶点坐标范围：xmin，xmax，ymin，ymax，zmin，zmax（z轴向上为正）
	double Density;								//密度：重力正演时使用
	double CiHuaLv,Angle_I,Angle_D;				//磁化率,磁化倾角，磁化偏角：磁法正演时使用（单位：度）
	COLORREF Color;								//用于VTK显示是的变量
};

//球体模型结构体
struct GM_GlobalInfo
{
	double Radius;								//半径
	double Center[3];							//中心位置：x0，y0，z0，z0的绝对值也就是中心埋深
	double Density;								//密度：重力正演时使用
	double CiHuaLv,Angle_I,Angle_D;				//磁化率,磁化倾角，磁化偏角：磁法正演时使用（单位：度）
	double Resolution;							//用于vtk显示时的一个变量
	COLORREF Color;								//用于VTK显示是的变量
};

//规则几何体模型
struct GM_RegularGeometryInfo
{
	CGM_GRDDATAINFO grddatainfo;				//网格数据的结构体，包括了数据信息和坐标系信息
	vector<GM_CylinderInfo> cylinder_vec;		//记录圆柱体模型参数的向量
	vector<GM_CubeInfo> cube_vec;				//记录立方体模型参数的向量
	vector<GM_GlobalInfo> global_vec;			//记录球体模型参数的向量
};

//规则几何体模型正演结构体
struct RegularGeometry3DForward
{
	GM_RegularGeometryInfo model;						//模型体结构变量
	double EarthMag,EarthAngle_I,EarthAngle_D;			//记录地球磁场和磁偏角及磁倾角：正演磁的时候使用（单位：nT,度）
};
/*-----------------------------------------------将各个规则模型体的正演归在一起---------------------------------------------------
	计算的Vz单位为mGal，垂向梯度单位为E
	规则模型体重力三维正演：调用各自类型模型体的正演函数
	输入参数：grav-存放正演数据的二维数组（行对应x，列对应y，00元素对应x，y均最小的正演值，顺推），单位mGal`m,mGal,mGal/m,mGal/m/m
          forwardinfo-存放数据之外的关于观测系统和模型参数的结构体变量，这里是引用的形式，为了在正演函数里求出场值的最大值和最小值
          type-表示正演分量类型，比如FORWARD_Vz表示Vz的正演
----------------------------------------------------------------------------------------------------------------------------------*/
GMDLL int _3DRegularModel_grav(double** grav,RegularGeometry3DForward& forwardinfo,int type=FORWARD_Vz);
GMDLL int _3DRegularModel_mag(double** mag,RegularGeometry3DForward& forwardinfo,int type=FORWARD_Za);
GMDLL int _3DRegularModel_magT(double** mag,RegularGeometry3DForward& forwardinfo,int type=FORWARD_Ta);

/*-----------------------------------------------下面是各种规则形体的正演函数------------------------------------------------------*/
/*------------------
	有限长水平圆柱体
	----------------
	说明：水平圆柱体必须的参数物性参数不用说了，那是必须的；
	位置参数只需要给出两个端点的坐标即可；中心埋深是需要的；半径那也是必须的；长度可以不给；中心点坐标也不需要
------------------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------
	有限长水平圆柱体重力正演
	-------------------------
输入参数：grav-存放正演数据的二维数组（行对应x，列对应y，00元素对应x，y均最小的正演值，顺推），单位mGal`m,mGal,mGal/m,mGal/m/m
          forwardinfo-存放数据之外的关于观测系统和模型参数的结构体变量，这里是引用的形式，为了在正演函数里求出场值的最大值和最小值
          type-表示正演分量类型，比如FORWARD_Vz表示Vz的正演
-----------------------------------------------------------------------------------------------------------------------------------*/
GMDLL int _3DHorizontalCylinder_grav(double** grav,RegularGeometry3DForward& forwardinfo,int type=FORWARD_Vz);

/*
有限长水平圆柱体磁三分量正演
输入参数：mag-存放正演数据的二维数组（行对应x，列对应y，00元素对应x，y均最小的正演值，顺推）单位：nT
          forwardinfo-存放数据之外的关于观测系统和模型参数的结构体变量，这里是引用的形式，为了在正演函数里求出场值的最大值和最小值
          type-表示正演分量类型，比如FORWARD_Ta表示daiertaT的正演
-----------------------------------------------------------------------------------------------------------------------------------*/
GMDLL int _3DHorizontalCylinder_mag(double** mag,RegularGeometry3DForward& forwardinfo,int type=FORWARD_Za);

/*
有限长水平圆柱体磁Ta和模量等正演（需要调用三分量正演函数）单位：nT
输入参数：mag-存放正演数据的二维数组（行对应x，列对应y，00元素对应x，y均最小的正演值，顺推）
          forwardinfo-存放数据之外的关于观测系统和模型参数的结构体变量，这里是引用的形式，为了在正演函数里求出场值的最大值和最小值
          type-表示正演分量类型，比如FORWARD_Ta表示daiertaT的正演
-----------------------------------------------------------------------------------------------------------------------------------*/
GMDLL int _3DHorizontalCylinder_magT(double** mag,RegularGeometry3DForward& forwardinfo,int type=FORWARD_Ta);

/*----------------------------
	有限长直立线模型的重磁正演
	--------------------------
输入参数：grav-存放正演数据的二维数组（行对应x，列对应y，00元素对应x，y均最小的正演值，顺推），单位mGal`m,mGal,mGal/m,mGal/m/m
forwardinfo-存放数据之外的关于观测系统和模型参数的结构体变量，这里是引用的形式，为了在正演函数里求出场值的最大值和最小值
type-表示正演分量类型，比如FORWARD_Vz表示Vz的正演
-----------------------------------------------------------------------------------------------------------------------------------*/
GMDLL int _3DFiniteVercicalLine_mag(double** mag, RegularGeometry3DForward& forwardinfo, int type = FORWARD_Za);
GMDLL int _3DFiniteVercicalLine_grav(double** grav, RegularGeometry3DForward& forwardinfo, int type = FORWARD_Vz);
GMDLL int _3DFiniteVercicalLine_magT(double** mag, RegularGeometry3DForward& forwardinfo, int type = FORWARD_Ta);

/*----------------------------
有限长直立圆柱体模型的重磁正演
------------------------------
theta_num:0-2Pi的积分剖分个数
*/
//z轴向下为正，观测面所在高度如果在地面之上则为负值
GMDLL int _3DFiniteVercicalCylinder_grav(double** grav, RegularGeometry3DForward& forwardinfo, int type = FORWARD_Vz, int theta_num = 100);
GMDLL int _3DFiniteVercicalCylinder_mag(double** mag, RegularGeometry3DForward& forwardinfo, int type = FORWARD_Hax, int theta_num = 100);
GMDLL int _3DFiniteVercicalCylinder_magT(double** mag, RegularGeometry3DForward& forwardinfo, int type = FORWARD_Ta, int theta_num = 100);
GMDLL int _3DFiniteCylinder_Grav(double** grav, RegularGeometry3DForward& forwardinfo0, int type = FORWARD_Vz, int theta_num = 100);
GMDLL int _3DFiniteCylinder_mag(double** mag, RegularGeometry3DForward& forwardinfo0, int type = FORWARD_Hax, int theta_num = 100);
GMDLL int _3DFiniteCylinder_magT(double** mag, RegularGeometry3DForward& forwardinfo, int type = FORWARD_Ta, int theta_num = 100);
/*----------------------------
	有限长倾斜线模型的重磁正演
	--------------------------
输入参数：grav-存放正演数据的二维数组（行对应x，列对应y，00元素对应x，y均最小的正演值，顺推），单位mGal`m,mGal,mGal/m,mGal/m/m
forwardinfo-存放数据之外的关于观测系统和模型参数的结构体变量，这里是引用的形式，为了在正演函数里求出场值的最大值和最小值
type-表示正演分量类型，比如FORWARD_Vz表示Vz的正演
-----------------------------------------------------------------------------------------------------------------------------------*/
GMDLL int _3DFiniteLine_Grav(double** grav, RegularGeometry3DForward& forwardinfo, int type = FORWARD_Vz);
GMDLL int _3DFiniteLine_mag(double** mag, RegularGeometry3DForward& forwardinfo, int type=FORWARD_Za);
GMDLL int _3DFiniteLine_magT(double** mag, RegularGeometry3DForward& forwardinfo, int type = FORWARD_Ta);
/*--------
	立方体
	------
	-------------------------------------------------------------------------------------------------------------------------------*/
/*
立方体重力正演
输入参数：grav-存放正演数据的二维数组（行对应x，列对应y，00元素对应x，y均最小的正演值，顺推），单位mGal`m,mGal,mGal/m,mGal/m/m
          forwardinfo-存放数据之外的关于观测系统和模型参数的结构体变量，这里是引用的形式，为了在正演函数里求出场值的最大值和最小值
          type-表示正演分量类型，比如FORWARD_Vz表示Vz的正演
-----------------------------------------------------------------------------------------------------------------------------------*/
GMDLL int _3DCube_grav(double** grav,RegularGeometry3DForward& forwardinfo,int type=FORWARD_Vz);

/*
立方体磁三分量正演
输入参数：mag-存放正演数据的二维数组（行对应x，列对应y，00元素对应x，y均最小的正演值，顺推）单位：nT
          forwardinfo-存放数据之外的关于观测系统和模型参数的结构体变量，这里是引用的形式，为了在正演函数里求出场值的最大值和最小值
          type-表示正演分量类型，比如FORWARD_Ta表示daiertaT的正演
-----------------------------------------------------------------------------------------------------------------------------------*/
GMDLL int _3DCube_mag(double** mag,RegularGeometry3DForward& forwardinfo,int type=FORWARD_Za);

/*
立方体磁Ta和模量等正演（需要调用三分量正演函数）单位：nT
输入参数：mag-存放正演数据的二维数组（行对应x，列对应y，00元素对应x，y均最小的正演值，顺推）
          forwardinfo-存放数据之外的关于观测系统和模型参数的结构体变量，这里是引用的形式，为了在正演函数里求出场值的最大值和最小值
          type-表示正演分量类型，比如FORWARD_Ta表示daiertaT的正演
-----------------------------------------------------------------------------------------------------------------------------------*/
GMDLL int _3DCube_magT(double** mag,RegularGeometry3DForward& forwardinfo,int type=FORWARD_Ta);

/*------
	球体
	----
-------------------------------------------------------------------------------------------------------------------------------*/
/*
球体重力正演
输入参数：grav-存放正演数据的二维数组（行对应x，列对应y，00元素对应x，y均最小的正演值，顺推），单位mGal`m,mGal,mGal/m,mGal/m/m
          forwardinfo-存放数据之外的关于观测系统和模型参数的结构体变量，这里是引用的形式，为了在正演函数里求出场值的最大值和最小值
          type-表示正演分量类型，比如FORWARD_Vz表示Vz的正演
-------------------------------------------------------------------------------------------------------------------------------*/
GMDLL int _3DGlobal_grav(double** grav,RegularGeometry3DForward& forwardinfo,int type=FORWARD_Vz);
GMDLL int _3DGlobal_mag(double** mag,RegularGeometry3DForward& forwardinfo,int type=FORWARD_Za);
GMDLL int _3DGlobal_magT(double** mag,RegularGeometry3DForward& forwardinfo,int type=FORWARD_Ta);


#endif