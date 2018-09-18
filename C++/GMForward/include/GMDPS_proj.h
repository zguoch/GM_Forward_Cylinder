
//重磁数据处理系统用到的公共的数据结构及头文件等

#ifndef GMDPS_PROJ
#define GMDPS_PROJ

//#include <comdef.h>
#include "math.h"
#include <afx.h>
#include <atlstr.h>
#include <windows.h>
#include <vector>
using namespace std;

//定义常数
#define G 6.67E-11
#define PI 3.1415926535897932384626433832795
#define U0 (4.0*PI*1E-7)
#define EPS 1E-10
//------------------------------------------数据结构定义------------------------------------------------------------
//扩边方法
#define EXTENBOUNDARY_COS	0		//余弦衰减至0扩边
#define EXTENBOUNDARY_MinCurvature	1	//最小曲率扩边
#define EXTENBOUNDARY_Poly	2		//局部多项式扩边


//定义原始数据类型和转换数据类型
#define MAG_ORIGINAL_DT	0
#define MAG_ORIGINAL_DZ	1

#define MAG_TRANS_DX	0
#define MAG_TRANS_DY	1
#define MAG_TRANS_DZ	2
#define MAG_TRANS_GRAV	3
#define MAG_TRANS_MODULE_Ta	4
#define MAG_TRANS_MODULE_R	5
#define MAG_TRANS_MODULE_E	6
#define MAG_TRANS_MODULE_Q	7
#define MAG_TRANS_MODULE_L	8

//数据处理：磁场分量换算信息数据结构体
struct MagneticComponentTransStruct
{
	int OriginDataType, TransDataType;			//原始数据类型，和转换数据类型
	double MagAngle_D_Earth, MagAngle_I_Earth;//地磁偏角和倾角（度）
};
//数据处理：方向导数信息结构体
#define DERIVATIVE_NORTH	3
#define DERIVATIVE_EAST	2
#define DERIVATIVE_VERTICAL1	0
#define DERIVATIVE_VERTICAL2	1
#define DERIVATIVE_HORIZONTAL	4
struct DerivativeInfoStruct
{
	int DerivativeType;			//求导类型
	int DerivativeOrder;		//求导阶数
	double aerfa;				//求导方向与X轴（北向）正方向的夹角
};
//三维重磁数据类
/*使用这个类之前必须进行初始化，然后利用push函数来为data赋值*/
class CGM_GRDDATA
{
public:
	CGM_GRDDATA()
	{
		m_number_x = 0;
		m_number_y = 0;
		height_data = 0;
		for (int i = 0; i < 6; i++)
		{
			m_bounds[i] = 0;
		}
		m_data = new double[m_number_y*m_number_x];
	};
	~CGM_GRDDATA()
	{
		
	};

private:
	double* m_data;
	double m_bounds[6];
	int m_number_x, m_number_y;
	double height_data;//数据所在高度（观测高度，z轴向上为正）
public:
	int Initialize(int number_x, int number_y, double* bounds)
	{
		/*首先释放老的指针*/
		if (m_data)
		{
			delete m_data;
		}
		/*为变量赋值*/
		m_number_x = number_x;
		m_number_y = number_y;
		for (int i = 0; i < 6; i++)
		{
			m_bounds[i] = bounds[i];
		}
		/*重新申请数组*/
		m_data = new double[number_y*number_x];
		return 0;
	}

	int SetData(int rows, int cols,double data)
	{
		m_data[rows*m_number_x+cols] = data;
		return 0;
	}
	int GetHeight()
	{
		return height_data;
	}
	double* GetBounds()
	{
		return m_bounds;
	}

	double* GetData()
	{
		return m_data;
	}

	int GetDimension(int& number_x, int& number_y)
	{
		number_x = m_number_x;
		number_y = m_number_y;
		return 0;
	}

	int Delete()
	{
		if (m_data)
		{
			delete m_data;
		}
		return 0;
	}
};
//三维重磁数据信息（平面grd数据）类
class CGM_GRDDATAINFO
{
public:
	CGM_GRDDATAINFO()
	{
		m_AxisBounds[0]=0;
		m_AxisBounds[1]=1000;
		m_AxisBounds[2]=0;
		m_AxisBounds[3]=1200;
		m_AxisBounds[4]=-300;
		m_AxisBounds[5]=0;
		m_Ranges[0]=0;
		m_Ranges[1]=1;
		m_Number_x=128;
		m_Number_y=128;
		m_Height_data=0;
		UpdateDxDy();
	};
	~CGM_GRDDATAINFO(){};
public:
	char m_Discription[40];						//记录关于这个数据的有关描述，比如：重力Vz-正演-网格 等等
	double m_AxisBounds[6];						//记录xmin,xmax,ymin,ymax,zmin,zmax
	double m_Ranges[2];							//记录场的最小值和最大值
	double m_Height_data;							//数据高度：一般规定正常坐标系，向上为正
	double m_Dx,m_Dy;								//x，y方向的间距
private:
	int m_Number_x,m_Number_y;
public:
	int GetNumber_x()
	{
		//UpdateNumber();
		return m_Number_x;
	};
	int AutoGetNumber_x()
	{
		UpdateNumber();
		return m_Number_x;
	};
	int GetNumber_y()
	{
		//UpdateNumber();
		return m_Number_y;
	}
	int AutoGetNumber_y()
	{
		UpdateNumber();
		return m_Number_y;
	}
	void SetNumber(int number_x,int number_y)
	{
		m_Number_x=number_x;
		m_Number_y=number_y;
		UpdateDxDy();
	};
private:
	void UpdateNumber(void)
	{
		m_Number_x=(int)((m_AxisBounds[1]-m_AxisBounds[0])/m_Dx+1);
		m_Number_y=(int)((m_AxisBounds[3]-m_AxisBounds[2])/m_Dy+1);
		m_AxisBounds[1]=m_AxisBounds[0]+m_Dx*(m_Number_x-1);
		m_AxisBounds[3]=m_AxisBounds[2]+m_Dy*(m_Number_y-1);
	};
	void UpdateMax()
	{
		m_AxisBounds[1]=m_AxisBounds[0]+(m_Number_x-1)*m_Dx;
		m_AxisBounds[3]=m_AxisBounds[2]+(m_Number_y-1)*m_Dy;
	};
	void UpdateDxDy()
	{
		m_Dx=(m_AxisBounds[1]-m_AxisBounds[0])/(m_Number_x-1);
		m_Dy=(m_AxisBounds[3]-m_AxisBounds[2])/(m_Number_y-1);
	};

};

//矩阵大小
struct MatrixSize
{
	int rows;
	int cols;
};


//---------------------------共用的一些基本的函数---------------------------------------------------------------------
//为一个二维数组赋值为常数
void Assign_Array2(double** data,int rows,int cols,int indata=0);
//符号函数
int Sign(double number);
//申请二维数组
double** CreateArray2(int rows,int cols);
//释放二维数组内存
void DeleteArray2(double** array2,int rows,int cols);
//获取网格数据的最大值最小值
void GetGrdMinMax(double** data,CGM_GRDDATAINFO & grdinfo);

//M×N矩阵乘以N维向量:C=AB
void Mat_Multiply(double* A,double* B,double* C,int rows,int cols);
//两点之间的距离
double Distance(double* point1,double *point2,int dims);
double Length_Vector(double* vector0,const int dims);
//求向量最大值,最小值
double GetMax(double* data,const int dims);
double GetMin(double* data,const int dims);
//计算叉乘（三维坐标）
double* Cross(double* a,double* b);
//两向量的点乘(三维坐标)
double VectorDot(double *a,double *b,int dims=3);

#endif