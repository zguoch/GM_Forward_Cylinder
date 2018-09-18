// 重磁三维正演.cpp : 定义控制台应用程序的入口点。
//

//#include "stdafx.h"

#include "GM3DForward.h"
#include "GMFile.h"
#include "FFTN.h"
#include <comdef.h>
void _Grav_HorizontalCylinder_Test();			//测试水平有限长水平圆柱体
void _Mag_HorizontalCylinder_Test();
void _Grav_Cube_Test();							//立方体重力三维正演测试
void _Mag_Cube_Test();
void _Grav_Global_Test();						//球体重力三维正演
void _Mag_Global_Test();
void _Grav_RegularModelTest();					//测试规则模型的重力正演：包括水平圆柱体，立方体
void _Mag_RegularModelTest();
void TestFFT1();								//测试一维FFT
void TestFinitLine();
void TestFinitCylinder();
//测试隐式调用DLL
//extern __declspec(dllimport) int ADD();
int _tmain()
{
	//_Grav_HorizontalCylinder_Test();
	//_Mag_HorizontalCylinder_Test();

	//_Grav_Cube_Test();
	//_Mag_Cube_Test();

	//_Grav_Global_Test();
	//_Mag_Global_Test();

	//_Grav_RegularModelTest();
	//_Mag_RegularModelTest();
	//TestFFT1();
	//TestFinitLine();
	//_Mag_HorizontalCylinder_Test();
	TestFinitCylinder();
	system("PAUSE");
	return 0;
}
void TestFinitCylinder()
{
	RegularGeometry3DForward forwardinfo;
	//观测系统参数
	forwardinfo.model.grddatainfo.m_Dx = 1;
	forwardinfo.model.grddatainfo.m_Dy = 1;
	forwardinfo.model.grddatainfo.m_AxisBounds[0] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[1] = 300;
	forwardinfo.model.grddatainfo.m_AxisBounds[2] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[3] = 400;
	//地磁场参数
	forwardinfo.EarthMag = 50000;//nT
	forwardinfo.EarthAngle_D = 25;
	forwardinfo.EarthAngle_I = 90;
	//模型参数
	GM_CylinderInfo cylinderinfo;
	cylinderinfo.Radius = 10;
	cylinderinfo.Density = 1;
	cylinderinfo.CiHuaLv = 0.25;
	cylinderinfo.Angle_D = forwardinfo.EarthAngle_D;//无剩磁
	cylinderinfo.Angle_I = forwardinfo.EarthAngle_I;//无剩磁
	cylinderinfo.Pos1[0] = 150; cylinderinfo.Pos1[1] = 150; cylinderinfo.Pos1[2] = 15;
	cylinderinfo.Pos2[0] = 150; cylinderinfo.Pos2[1] = 150; cylinderinfo.Pos2[2] = 300;
	forwardinfo.model.cylinder_vec.push_back(cylinderinfo);
	//第二个模型
	/*cylinderinfo.Pos1[0] = 210; cylinderinfo.Pos1[1] = 100; cylinderinfo.Pos1[2] = 80;
	cylinderinfo.Pos2[0] = 280; cylinderinfo.Pos2[1] = 300; cylinderinfo.Pos2[2] = 10;
	forwardinfo.model.cylinder_vec.push_back(cylinderinfo);*/
	//二维数组存放正演数据
	double** grav = CreateArray2(forwardinfo.model.grddatainfo.AutoGetNumber_y(), forwardinfo.model.grddatainfo.AutoGetNumber_x());
	_3DFiniteCylinder_Grav(grav, forwardinfo, FORWARD_Vz);
	//_3DFiniteCylinder_mag(grav, forwardinfo, FORWARD_Hay);
	//_3DFiniteCylinder_magT(grav, forwardinfo, FORWARD_Module);
	CGMFile gmfile;
	gmfile.SetFileName(_T("测试目录\\倾斜圆柱体Vz.grd"));
	gmfile.SaveGrd(grav, forwardinfo.model.grddatainfo);
}
void TestFinitLine()
{
	RegularGeometry3DForward forwardinfo;
	//观测系统参数
	forwardinfo.model.grddatainfo.m_Dx = 1;
	forwardinfo.model.grddatainfo.m_Dy = 1;
	forwardinfo.model.grddatainfo.m_AxisBounds[0] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[1] = 300;
	forwardinfo.model.grddatainfo.m_AxisBounds[2] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[3] = 400;
	//地磁场参数
	forwardinfo.EarthMag = 50000;//nT
	forwardinfo.EarthAngle_D = -45;
	forwardinfo.EarthAngle_I = 60;
	//模型参数
	GM_CylinderInfo cylinderinfo;
	cylinderinfo.Radius = 4;
	cylinderinfo.Density = 1;
	cylinderinfo.CiHuaLv = 0.25;
	cylinderinfo.Angle_D = forwardinfo.EarthAngle_D;//无剩磁
	cylinderinfo.Angle_I = forwardinfo.EarthAngle_I;//无剩磁
	cylinderinfo.Pos1[0] = 100; cylinderinfo.Pos1[1] = 200; cylinderinfo.Pos1[2] = 4;
	cylinderinfo.Pos2[0] = 200; cylinderinfo.Pos2[1] = 200; cylinderinfo.Pos2[2] = 4;
	forwardinfo.model.cylinder_vec.push_back(cylinderinfo);
	//第二个模型
	/*cylinderinfo.Pos1[0] = 210; cylinderinfo.Pos1[1] = 100; cylinderinfo.Pos1[2] = 80;
	cylinderinfo.Pos2[0] = 280; cylinderinfo.Pos2[1] = 300; cylinderinfo.Pos2[2] = 10;
	forwardinfo.model.cylinder_vec.push_back(cylinderinfo);*/
	//二维数组存放正演数据
	double** grav = CreateArray2(forwardinfo.model.grddatainfo.AutoGetNumber_y(), forwardinfo.model.grddatainfo.AutoGetNumber_x());
	//_3DFiniteLine_Grav(grav, forwardinfo, FORWARD_Vyz);
	_3DFiniteLine_magT(grav, forwardinfo, FORWARD_Ta);
	CGMFile gmfile;
	gmfile.SetFileName(_T("测试目录\\倾斜圆柱体Ta.grd"));
	gmfile.SaveGrd(grav, forwardinfo.model.grddatainfo);
}
void TestFFT1()
{
	//1. 读取数据
	FILE* fp = NULL;
	if ((fp = fopen("65和25HZ的正弦函数叠加2048个(采样频率140HZ) (2).txt", "r")) == NULL)
	{
		printf("打开文件失败\n");
		return;
	}
	double x, y;
	vector<double>xvector, yvector;
	while (!feof(fp))
	{
		fscanf(fp, "%lf %lf", &x, &y);
		xvector.push_back(x);
		yvector.push_back(y);
	}
	if (xvector[xvector.size() - 1] == xvector[xvector.size() - 2] && yvector[yvector.size() - 1] == yvector[yvector.size() - 2])
	{
		xvector.pop_back(); yvector.pop_back();
	}
	fclose(fp);
	//2. 数据个数
	int N = xvector.size();
	printf("变换数: %d\n", N);
	//------------
	fftw_complex *out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
	//3. 调用fft1d
	FFT1d(yvector, out);
	//4. 计算频谱
	vector<double>f, spectrum;
	GetSpectrum1d(N, fabs(xvector[1] - xvector[0]), out, f, spectrum);
	//5. 计算逆变换
	vector<double>ifftoutvector;
	IFFT1d(N, out, ifftoutvector);

	//4. 输出
	FILE* fpout = NULL;
	if ((fpout = fopen("逆变换.txt", "w")) == NULL)
	{
		printf("打开输出文件失败\n");
		return;
	}
	for (int i = 0; i < ifftoutvector.size(); i++)
	{
		fprintf(fpout, "%.7f\t%.7f\n", xvector[i], ifftoutvector[i]);
	}
	fclose(fpout);



	//释放fftw指针
	fftw_free(out);

	//7.2 这是傅里叶逆变换

	//7.3 打印逆变换的值
}
void _Mag_Global_Test()
{
	RegularGeometry3DForward forwardinfo;
	//观测系统参数
	forwardinfo.model.grddatainfo.m_Dx = 10;
	forwardinfo.model.grddatainfo.m_Dy = 10;
	forwardinfo.model.grddatainfo.m_AxisBounds[0] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[1] = 1000;
	forwardinfo.model.grddatainfo.m_AxisBounds[2] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[3] = 1000;
	//地磁场参数
	forwardinfo.EarthMag = 50000;//nT
	forwardinfo.EarthAngle_D = 45;
	forwardinfo.EarthAngle_I = 45;
	//模型参数
	//--------------------第一个模型-------------------------------
	GM_GlobalInfo globalinfo;
	globalinfo.Radius = 100;
	globalinfo.Center[0] = 300;
	globalinfo.Center[1] = 600;
	globalinfo.Center[2] = -130;
	globalinfo.Angle_D = 45;
	globalinfo.Angle_I = 45;
	globalinfo.CiHuaLv = 0.2;

	forwardinfo.model.global_vec.push_back(globalinfo);

	//二维数组存放正演数据
	double** mag = CreateArray2(forwardinfo.model.grddatainfo.GetNumber_y(), forwardinfo.model.grddatainfo.GetNumber_x());
	_3DGlobal_mag(mag, forwardinfo, FORWARD_Za);//正演磁力
	//_3DGlobal_magT(mag,forwardinfo,FORWARD_Ta);
	CGMFile gmfile;
	gmfile.SetFileName(_T("测试目录\\单个球体Za.grd"));
	gmfile.SaveGrd(mag, forwardinfo.model.grddatainfo);
	//释放指针资源
	DeleteArray2(mag, forwardinfo.model.grddatainfo.GetNumber_y(), forwardinfo.model.grddatainfo.GetNumber_x());
}
void _Grav_Global_Test()
{
	RegularGeometry3DForward forwardinfo;
	//观测系统参数
	forwardinfo.model.grddatainfo.m_Dx = 10;
	forwardinfo.model.grddatainfo.m_Dy = 10;
	forwardinfo.model.grddatainfo.m_AxisBounds[0] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[1] = 1000;
	forwardinfo.model.grddatainfo.m_AxisBounds[2] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[3] = 1000;
	//地磁场参数
	forwardinfo.EarthMag = 50000;//nT
	forwardinfo.EarthAngle_D = 30;
	forwardinfo.EarthAngle_I = 45;
	//模型参数
	//--------------------第一个模型-------------------------------
	GM_GlobalInfo globalinfo;
	globalinfo.Radius = 100;
	globalinfo.Center[0] = 300;
	globalinfo.Center[1] = 600;
	globalinfo.Center[2] = -130;
	globalinfo.Density = 1;

	forwardinfo.model.global_vec.push_back(globalinfo);

	//二维数组存放正演数据
	double** grav = CreateArray2(forwardinfo.model.grddatainfo.GetNumber_y(), forwardinfo.model.grddatainfo.GetNumber_x());
	_3DGlobal_grav(grav, forwardinfo, FORWARD_Vyy);//正演重力
	CGMFile gmfile;
	gmfile.SetFileName(_T("测试目录\\单个球体Vyy.grd"));
	gmfile.SaveGrd(grav, forwardinfo.model.grddatainfo);
	//释放指针资源
	DeleteArray2(grav, forwardinfo.model.grddatainfo.GetNumber_y(), forwardinfo.model.grddatainfo.GetNumber_x());
}
void _Mag_RegularModelTest()
{
	RegularGeometry3DForward forwardinfo;
	//观测系统参数
	forwardinfo.model.grddatainfo.m_Dx = 5;
	forwardinfo.model.grddatainfo.m_Dy = 5;
	forwardinfo.model.grddatainfo.m_AxisBounds[0] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[1] = 1000;
	forwardinfo.model.grddatainfo.m_AxisBounds[2] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[3] = 1200;
	//地磁场参数
	forwardinfo.EarthMag = 50000;//nT
	forwardinfo.EarthAngle_D = 45;
	forwardinfo.EarthAngle_I = 45;
	//模型参数
	//--------水平圆柱体---------------------------
	GM_CylinderInfo cylinderinfo;
	cylinderinfo.CiHuaLv = 0.2;
	cylinderinfo.Angle_D = 5;
	cylinderinfo.Angle_I = 45;
	cylinderinfo.Radius = 50;
	cylinderinfo.Pos1[0] = 500; cylinderinfo.Pos1[1] = 800; cylinderinfo.Pos1[2] = -70;
	cylinderinfo.Pos2[0] = 500; cylinderinfo.Pos2[1] = 1100; cylinderinfo.Pos2[2] = -70;
	forwardinfo.model.cylinder_vec.push_back(cylinderinfo);

	//--------立方体----------------------------------------
	GM_CubeInfo cubeinfo;
	cubeinfo.CiHuaLv = 0.2;
	cubeinfo.Angle_D = 5;
	cubeinfo.Angle_I = 45;
	cubeinfo.bounds[0] = 300;
	cubeinfo.bounds[1] = 700;
	cubeinfo.bounds[2] = 400;
	cubeinfo.bounds[3] = 600;
	cubeinfo.bounds[4] = -350;
	cubeinfo.bounds[5] = -50;

	forwardinfo.model.cube_vec.push_back(cubeinfo);
	//---------球体模型------------------------------------------
	GM_GlobalInfo globalinfo;
	globalinfo.Radius = 100;
	globalinfo.Center[0] = 300;
	globalinfo.Center[1] = 600;
	globalinfo.Center[2] = -130;
	globalinfo.Angle_D = 45;
	globalinfo.Angle_I = 45;
	globalinfo.CiHuaLv = 0.2;

	forwardinfo.model.global_vec.push_back(globalinfo);
	//二维数组存放正演数据
	double** grav = CreateArray2(forwardinfo.model.grddatainfo.GetNumber_y(), forwardinfo.model.grddatainfo.GetNumber_x());
	//_3DRegularModel_mag(grav,forwardinfo,FORWARD_Hax);
	_3DRegularModel_magT(grav, forwardinfo, FORWARD_Ta);
	CGMFile gmfile;
	gmfile.SetFileName(_T("测试目录\\水平圆柱体+l立方体+球体Ta.grd"));
	gmfile.SaveGrd(grav, forwardinfo.model.grddatainfo);
}
void _Grav_RegularModelTest()
{
	RegularGeometry3DForward forwardinfo;
	//观测系统参数
	forwardinfo.model.grddatainfo.m_Dx = 5;
	forwardinfo.model.grddatainfo.m_Dy = 5;
	forwardinfo.model.grddatainfo.m_AxisBounds[0] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[1] = 1000;
	forwardinfo.model.grddatainfo.m_AxisBounds[2] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[3] = 1200;
	//模型参数
	//--------水平圆柱体---------------------------
	GM_CylinderInfo cylinderinfo;
	cylinderinfo.Radius = 50;
	cylinderinfo.Density = 1;
	cylinderinfo.Pos1[0] = 500; cylinderinfo.Pos1[1] = 800; cylinderinfo.Pos1[2] = -70;
	cylinderinfo.Pos2[0] = 500; cylinderinfo.Pos2[1] = 1100; cylinderinfo.Pos2[2] = -70;
	forwardinfo.model.cylinder_vec.push_back(cylinderinfo);

	//--------立方体----------------------------------------
	GM_CubeInfo cubeinfo;
	cubeinfo.Density = 1;
	cubeinfo.bounds[0] = 300;
	cubeinfo.bounds[1] = 700;
	cubeinfo.bounds[2] = 400;
	cubeinfo.bounds[3] = 600;
	cubeinfo.bounds[4] = -350;
	cubeinfo.bounds[5] = -50;
	forwardinfo.model.cube_vec.push_back(cubeinfo);
	//--------球体-----------------------------------------
	GM_GlobalInfo globalinfo;
	globalinfo.Radius = 100;
	globalinfo.Center[0] = 300;
	globalinfo.Center[1] = 600;
	globalinfo.Center[2] = -130;
	globalinfo.Density = 1;

	forwardinfo.model.global_vec.push_back(globalinfo);
	//二维数组存放正演数据
	double** grav = CreateArray2(forwardinfo.model.grddatainfo.GetNumber_y(), forwardinfo.model.grddatainfo.GetNumber_x());
	_3DRegularModel_grav(grav, forwardinfo, FORWARD_Vzz);

	CGMFile gmfile;
	gmfile.SetFileName(_T("测试目录\\水平圆柱体+l立方体+球体vzz.grd"));
	gmfile.SaveGrd(grav, forwardinfo.model.grddatainfo);
}
void _Mag_Cube_Test()
{
	RegularGeometry3DForward forwardinfo;
	//观测系统参数
	forwardinfo.model.grddatainfo.m_Dx = 10;
	forwardinfo.model.grddatainfo.m_Dy = 10;
	forwardinfo.model.grddatainfo.m_AxisBounds[0] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[1] = 1000;
	forwardinfo.model.grddatainfo.m_AxisBounds[2] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[3] = 1000;
	//地磁场参数
	forwardinfo.EarthMag = 50000;//nT
	forwardinfo.EarthAngle_D = 5;
	forwardinfo.EarthAngle_I = 45;
	//模型参数
	//--------------------第一个模型-------------------------------
	GM_CubeInfo cubeinfo;
	cubeinfo.CiHuaLv = 0.2;
	cubeinfo.Angle_D = 5;
	cubeinfo.Angle_I = 45;
	cubeinfo.bounds[0] = 300;
	cubeinfo.bounds[1] = 700;
	cubeinfo.bounds[2] = 400;
	cubeinfo.bounds[3] = 600;
	cubeinfo.bounds[4] = -350;
	cubeinfo.bounds[5] = -50;

	forwardinfo.model.cube_vec.push_back(cubeinfo);

	//二维数组存放正演数据
	double** mag = CreateArray2(forwardinfo.model.grddatainfo.GetNumber_y(), forwardinfo.model.grddatainfo.GetNumber_x());
	//_3DCube_mag(mag,forwardinfo,FORWARD_Za);//正演磁力
	_3DCube_magT(mag, forwardinfo, FORWARD_Module);
	CGMFile gmfile;
	gmfile.SetFileName(_T("测试目录\\单个立方体模量.grd"));
	gmfile.SaveGrd(mag, forwardinfo.model.grddatainfo);
	//释放指针资源
	DeleteArray2(mag, forwardinfo.model.grddatainfo.GetNumber_y(), forwardinfo.model.grddatainfo.GetNumber_x());
}
void _Grav_Cube_Test()
{
	RegularGeometry3DForward forwardinfo;
	//观测系统参数
	forwardinfo.model.grddatainfo.m_Dx = 1;
	forwardinfo.model.grddatainfo.m_Dy = 1;
	forwardinfo.model.grddatainfo.m_AxisBounds[0] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[1] = 40;
	forwardinfo.model.grddatainfo.m_AxisBounds[2] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[3] = 30;
	//地磁场参数
	forwardinfo.EarthMag = 50000;//nT
	forwardinfo.EarthAngle_D = 30;
	forwardinfo.EarthAngle_I = 45;
	//模型参数
	//--------------------第一个模型-------------------------------
	GM_CubeInfo cubeinfo;
	cubeinfo.Density = 2.67;
	cubeinfo.bounds[0] = 10;
	cubeinfo.bounds[1] = 30;
	cubeinfo.bounds[2] = 10;
	cubeinfo.bounds[3] = 20;
	cubeinfo.bounds[4] = -20;
	cubeinfo.bounds[5] = -10;

	forwardinfo.model.cube_vec.push_back(cubeinfo);

	//二维数组存放正演数据
	double** grav = CreateArray2(forwardinfo.model.grddatainfo.GetNumber_y(), forwardinfo.model.grddatainfo.GetNumber_x());
	_3DCube_grav(grav, forwardinfo, FORWARD_Vz);//正演重力
	CGMFile gmfile;
	gmfile.SetFileName(_T("测试目录\\单个立方体Vz.grd"));
	gmfile.SaveGrd(grav, forwardinfo.model.grddatainfo);
	printf("计算完毕\n");
	//释放指针资源
	DeleteArray2(grav, forwardinfo.model.grddatainfo.GetNumber_y(), forwardinfo.model.grddatainfo.GetNumber_x());
}
void _Mag_HorizontalCylinder_Test()
{
	RegularGeometry3DForward forwardinfo;
	//观测系统参数
	forwardinfo.model.grddatainfo.m_Dx = 1;
	forwardinfo.model.grddatainfo.m_Dy = 1;
	forwardinfo.model.grddatainfo.m_AxisBounds[0] = -0;//北向切面-90
	forwardinfo.model.grddatainfo.m_AxisBounds[1] = 300;//300
	forwardinfo.model.grddatainfo.m_AxisBounds[2] = -0;//东向切面-90
	forwardinfo.model.grddatainfo.m_AxisBounds[3] = 400;//400 300
	//地磁场参数
	forwardinfo.EarthMag = 50000;//nT
	forwardinfo.EarthAngle_D = 25;
	forwardinfo.EarthAngle_I = 60;
	
	//模型参数
	//--------------------第一个模型-------------------------------
	GM_CylinderInfo cylinderinfo;
	cylinderinfo.Radius = 50;
	cylinderinfo.CiHuaLv = 0.25;
	cylinderinfo.Angle_D = forwardinfo.EarthAngle_D;//无剩磁
	cylinderinfo.Angle_I = forwardinfo.EarthAngle_I;//无剩磁
	cylinderinfo.Density = 1;//g/cm3
	cylinderinfo.Pos1[0] = 150; cylinderinfo.Pos1[1] = 200; cylinderinfo.Pos1[2] = 10;
	cylinderinfo.Pos2[0] = 150; cylinderinfo.Pos2[1] = 200; cylinderinfo.Pos2[2] = 200;
	forwardinfo.model.cylinder_vec.push_back(cylinderinfo);
	//----------------------------第二个模型-------------------------------------------------
	//cylinderinfo.Pos1[0] = 170; cylinderinfo.Pos1[1] = 179; cylinderinfo.Pos1[2] = 10;
	//cylinderinfo.Pos2[0] = 253; cylinderinfo.Pos2[1] = 95; cylinderinfo.Pos2[2] = 10;
	//forwardinfo.model.cylinder_vec.push_back(cylinderinfo);
	//---------------------------------------------------------------------------------------
	//二维数组存放正演数据
	double** mag = CreateArray2(forwardinfo.model.grddatainfo.AutoGetNumber_y(), forwardinfo.model.grddatainfo.AutoGetNumber_x());
	//_3DHorizontalCylinder_mag(mag,forwardinfo,FORWARD_Za);//正演三分量
	//_3DHorizontalCylinder_magT(mag, forwardinfo, FORWARD_Ta);//正演Ta
	//_3DFiniteVercicalLine_grav(mag, forwardinfo, FORWARD_Vx);
	//_3DFiniteVercicalLine_magT(mag, forwardinfo, FORWARD_Module);
	_3DFiniteVercicalCylinder_grav(mag, forwardinfo, FORWARD_Vx);
	//_3DFiniteVercicalCylinder_mag(mag, forwardinfo, FORWARD_Za);
	//_3DFiniteVercicalCylinder_magT(mag, forwardinfo, FORWARD_Ta);
	CGMFile gmfile;
	gmfile.SetFileName(_T("测试目录\\有限长直立圆柱体Vx.grd"));
	gmfile.SaveGrd(mag, forwardinfo.model.grddatainfo);
	//释放指针资源
	DeleteArray2(mag, forwardinfo.model.grddatainfo.AutoGetNumber_y(), forwardinfo.model.grddatainfo.AutoGetNumber_x());
}
void _Grav_HorizontalCylinder_Test()
{
	RegularGeometry3DForward forwardinfo;
	//观测系统参数
	forwardinfo.model.grddatainfo.m_Dx = 5;
	forwardinfo.model.grddatainfo.m_Dy = 5;
	forwardinfo.model.grddatainfo.m_AxisBounds[0] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[1] = 300;
	forwardinfo.model.grddatainfo.m_AxisBounds[2] = 0;
	forwardinfo.model.grddatainfo.m_AxisBounds[3] = 400;
	//模型参数
	GM_CylinderInfo cylinderinfo;
	cylinderinfo.Radius = 0.4;
	cylinderinfo.Density = 1;
	cylinderinfo.Pos1[0] = 100; cylinderinfo.Pos1[1] = 195; cylinderinfo.Pos1[2] = 10;
	cylinderinfo.Pos2[0] = 218; cylinderinfo.Pos2[1] = 275; cylinderinfo.Pos2[2] = 10;
	forwardinfo.model.cylinder_vec.push_back(cylinderinfo);
	//----------------------------第二个模型-------------------------------------------------
	cylinderinfo.Pos1[0] = 170; cylinderinfo.Pos1[1] = 179; cylinderinfo.Pos1[2] = 10;
	cylinderinfo.Pos2[0] = 253; cylinderinfo.Pos2[1] = 179;
	forwardinfo.model.cylinder_vec.push_back(cylinderinfo); cylinderinfo.Pos2[2] = 10;
	//---------------------------------------------------------------------------------------
	forwardinfo.model.cylinder_vec.push_back(cylinderinfo);
	//二维数组存放正演数据
	double** grav = CreateArray2(forwardinfo.model.grddatainfo.GetNumber_y(), forwardinfo.model.grddatainfo.GetNumber_x());
	_3DHorizontalCylinder_grav(grav, forwardinfo, FORWARD_Vzz);

	CGMFile gmfile;
	gmfile.SetFileName(_T("测试目录\\水平圆柱体vzz.grd"));
	gmfile.SaveGrd(grav, forwardinfo.model.grddatainfo);
}