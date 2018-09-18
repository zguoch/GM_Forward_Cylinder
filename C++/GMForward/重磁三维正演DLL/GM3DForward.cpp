#include "GM3DForward.h"
//规则几何形体的重磁三维正演

#include "GMDPS_proj.h"

int _3DRegularModel_magT(double** mag,RegularGeometry3DForward& forwardinfo,int type)
{
	int number_x=forwardinfo.model.grddatainfo.GetNumber_x();
	int number_y=forwardinfo.model.grddatainfo.GetNumber_y();
	double** mag_cylinder=CreateArray2(number_y,number_x);
	double** mag_cube=CreateArray2(number_y,number_x);
	double** mag_global=CreateArray2(number_y,number_x);
	Assign_Array2(mag_cube,number_y,number_x,0);
	Assign_Array2(mag_cylinder,number_y,number_x,0);
	Assign_Array2(mag_global,number_y,number_x,0);
	Assign_Array2(mag,number_y,number_x,0);
	
	//调用水平圆柱体正演函数
	_3DHorizontalCylinder_magT(mag_cylinder,forwardinfo,type);
	//调用立方体正演函数
	_3DCube_magT(mag_cube,forwardinfo,type);
	//调用球体正演函数
	_3DGlobal_magT(mag_global,forwardinfo,type);
	//将以上几种不同类型的模型正演值相加
	for (int i = 0; i < number_y; i++)
	{
		for (int j = 0; j < number_x; j++)
		{
			mag[i][j]=mag_cube[i][j]+mag_cylinder[i][j]+mag_global[i][j];
		}
	}
	GetGrdMinMax(mag,forwardinfo.model.grddatainfo);
	//销毁二维数组
	DeleteArray2(mag_cylinder,number_y,number_x);
	DeleteArray2(mag_cube,number_y,number_x);
	DeleteArray2(mag_global,number_y,number_x);
	return 1;
}
int _3DRegularModel_mag(double** mag,RegularGeometry3DForward& forwardinfo,int type)
{
	int number_x=forwardinfo.model.grddatainfo.GetNumber_x();
	int number_y=forwardinfo.model.grddatainfo.GetNumber_y();
	double** mag_cylinder=CreateArray2(number_y,number_x);
	double** mag_cube=CreateArray2(number_y,number_x);
	double** mag_global=CreateArray2(number_y,number_x);
	Assign_Array2(mag_cube,number_y,number_x,0);
	Assign_Array2(mag_cylinder,number_y,number_x,0);
	Assign_Array2(mag_global,number_y,number_x,0);
	Assign_Array2(mag,number_y,number_x,0);
	
	//调用水平圆柱体正演函数
	_3DHorizontalCylinder_mag(mag_cylinder,forwardinfo,type);
	//调用立方体正演函数
	_3DCube_mag(mag_cube,forwardinfo,type);
	//调用球体正演函数
	_3DGlobal_mag(mag_global,forwardinfo,type);
	//将以上几种不同类型的模型正演值相加
	for (int i = 0; i < number_y; i++)
	{
		for (int j = 0; j < number_x; j++)
		{
			mag[i][j]=mag_cube[i][j]+mag_cylinder[i][j]+mag_global[i][j];
		}
	}
	GetGrdMinMax(mag,forwardinfo.model.grddatainfo);
	//销毁二维数组
	DeleteArray2(mag_cylinder,number_y,number_x);
	DeleteArray2(mag_cube,number_y,number_x);
	DeleteArray2(mag_global,number_y,number_x);
	return 1;
}
int _3DRegularModel_grav(double** grav,RegularGeometry3DForward& forwardinfo,int type)
{
	int number_x=forwardinfo.model.grddatainfo.GetNumber_x();
	int number_y=forwardinfo.model.grddatainfo.GetNumber_y();
	double** grav_cylinder=CreateArray2(number_y,number_x);
	double** grav_cube=CreateArray2(number_y,number_x);
	double** grav_global=CreateArray2(number_y,number_x);
	Assign_Array2(grav_cube,number_y,number_x,0);
	Assign_Array2(grav_cylinder,number_y,number_x,0);
	Assign_Array2(grav_global,number_y,number_x,0);
	Assign_Array2(grav,number_y,number_x,0);
	
	//调用水平圆柱体正演函数
	//_3DHorizontalCylinder_grav(grav_cylinder,forwardinfo,type);//线模型
	_3DFiniteCylinder_Grav(grav_cylinder, forwardinfo, type);
	//调用立方体正演函数
	_3DCube_grav(grav_cube,forwardinfo,type);
	//调用球体正演函数
	_3DGlobal_grav(grav_global,forwardinfo,type);
	//将以上几种不同类型的模型正演值相加
	for (int i = 0; i < number_y; i++)
	{
		for (int j = 0; j < number_x; j++)
		{
			grav[i][j]=grav_cube[i][j]+grav_cylinder[i][j]+grav_global[i][j];
		}
	}
	//获取数据最大值最小值
	GetGrdMinMax(grav,forwardinfo.model.grddatainfo);
	//销毁二维数组
	DeleteArray2(grav_cylinder,number_y,number_x);
	DeleteArray2(grav_cube,number_y,number_x);
	DeleteArray2(grav_global,number_y,number_x);
	return 1;
}

int _3DGlobal_grav(double** grav,RegularGeometry3DForward& forwardinfo,int type)
{
	//首先将grav赋值为0，后面累加
	int number_x=forwardinfo.model.grddatainfo.GetNumber_x();
	int number_y=forwardinfo.model.grddatainfo.GetNumber_y();
	Assign_Array2(grav,number_y,number_x,0);
	
	//坐标范围
	double xmin=forwardinfo.model.grddatainfo.m_AxisBounds[0];
	double xmax=forwardinfo.model.grddatainfo.m_AxisBounds[1];
	double ymin=forwardinfo.model.grddatainfo.m_AxisBounds[2];
	double ymax=forwardinfo.model.grddatainfo.m_AxisBounds[3];
	double dx=forwardinfo.model.grddatainfo.m_Dx;
	double dy=forwardinfo.model.grddatainfo.m_Dy;
	double forwardheight=forwardinfo.model.grddatainfo.m_Height_data;
	switch (type)
	{
	case FORWARD_V:
		{
			for (int k = 0; k < (int)forwardinfo.model.global_vec.size(); k++)
			{
				
				double R=forwardinfo.model.global_vec[k].Radius;
				double x0=forwardinfo.model.global_vec[k].Center[0];
				double y0=forwardinfo.model.global_vec[k].Center[1];
				double D=fabs(forwardinfo.model.global_vec[k].Center[2]-forwardheight);
				double density=forwardinfo.model.global_vec[k].Density;
				double M=4.0/3.0*PI*pow(R,3.0)*density;
				double x,y;
				double temp;
				for (int i=0;i< number_y;i++)
				{
					y=ymin+i*dy-y0;
					for (int j=0;j<number_x;j++)
					{
						x=xmin+j*dx-x0;
						temp=G*M*pow(x*x+y*y+D*D,-1.0/2.0);
						grav[i][j]+=temp*1E8;					//mGal`m
					}
				}
			}
		}
		break;
	case FORWARD_Vz:
		{
			for (int k = 0; k < (int)forwardinfo.model.global_vec.size(); k++)
			{
				double R=forwardinfo.model.global_vec[k].Radius;
				double x0=forwardinfo.model.global_vec[k].Center[0];
				double y0=forwardinfo.model.global_vec[k].Center[1];
				double D=fabs(forwardinfo.model.global_vec[k].Center[2]-forwardheight);
				double density=forwardinfo.model.global_vec[k].Density;
				double M=4.0/3.0*PI*pow(R,3.0)*density;
				double x,y;
				double temp;
				for (int i=0;i< number_y;i++)
				{
					y=ymin+i*dy-y0;
					for (int j=0;j<number_x;j++)
					{
						x=xmin+j*dx-x0;
						temp=G*M*D*pow(x*x+y*y+D*D,-3.0/2.0);
						grav[i][j]+=temp*1E8;					//mGal
					}
				}
			}
		}
		break;
	case FORWARD_Vxx:
		{
			for (int k = 0; k < (int)forwardinfo.model.global_vec.size(); k++)
			{
				double R=forwardinfo.model.global_vec[k].Radius;
				double x0=forwardinfo.model.global_vec[k].Center[0];
				double y0=forwardinfo.model.global_vec[k].Center[1];
				double D=fabs(forwardinfo.model.global_vec[k].Center[2]-forwardheight);
				double density=forwardinfo.model.global_vec[k].Density;
				double M=4.0/3.0*PI*pow(R,3.0)*density;
				double x,y;
				double temp,temp2;
				for (int i=0;i< number_y;i++)
				{
					y=ymin+i*dy-y0;
					for (int j=0;j<number_x;j++)
					{
						x=xmin+j*dx-x0;
						temp2=x*x+y*y+D*D;
						temp=(2.0*x*x-y*y-D*D)/pow(temp2,2.5);
						grav[i][j]+=temp*G*M*1E12;					//E
					}
				}
			}
		}
		break;
	case FORWARD_Vyy:
		{
			for (int k = 0; k < (int)forwardinfo.model.global_vec.size(); k++)
			{
				double R=forwardinfo.model.global_vec[k].Radius;
				double x0=forwardinfo.model.global_vec[k].Center[0];
				double y0=forwardinfo.model.global_vec[k].Center[1];
				double D=fabs(forwardinfo.model.global_vec[k].Center[2]-forwardheight);
				double density=forwardinfo.model.global_vec[k].Density;
				double M=4.0/3.0*PI*pow(R,3.0)*density;
				double x,y;
				double temp,temp2;
				for (int i=0;i< number_y;i++)
				{
					y=ymin+i*dy-y0;
					for (int j=0;j<number_x;j++)
					{
						x=xmin+j*dx-x0;
						temp2=x*x+y*y+D*D;
						temp=(2.0*y*y-x*x-D*D)/pow(temp2,2.5);
						grav[i][j]+=temp*G*M*1E12;					//E
					}
				}
			}
		}
		break;
	case FORWARD_Vzz:
		{
			for (int k = 0; k < (int)forwardinfo.model.global_vec.size(); k++)
			{
				double R=forwardinfo.model.global_vec[k].Radius;
				double x0=forwardinfo.model.global_vec[k].Center[0];
				double y0=forwardinfo.model.global_vec[k].Center[1];
				double D=fabs(forwardinfo.model.global_vec[k].Center[2]-forwardheight);
				double density=forwardinfo.model.global_vec[k].Density;
				double M=4.0/3.0*PI*pow(R,3.0)*density;
				double x,y;
				double temp,temp2;
				for (int i=0;i< number_y;i++)
				{
					y=ymin+i*dy-y0;
					for (int j=0;j<number_x;j++)
					{
						x=xmin+j*dx-x0;
						temp2=x*x+y*y+D*D;
						temp=(2.0*D*D-x*x-y*y)/pow(temp2,2.5);
						grav[i][j]+=temp*G*M*1E12;					//E
					}
				}
			}
		}
		break;
	case FORWARD_Vxy:
		{
			for (int k = 0; k < (int)forwardinfo.model.global_vec.size(); k++)
			{
				double R=forwardinfo.model.global_vec[k].Radius;
				double x0=forwardinfo.model.global_vec[k].Center[0];
				double y0=forwardinfo.model.global_vec[k].Center[1];
				double D=fabs(forwardinfo.model.global_vec[k].Center[2]-forwardheight);
				double density=forwardinfo.model.global_vec[k].Density;
				double M=4.0/3.0*PI*pow(R,3.0)*density;
				double x,y;
				double temp,temp2;
				for (int i=0;i< number_y;i++)
				{
					y=ymin+i*dy-y0;
					for (int j=0;j<number_x;j++)
					{
						x=xmin+j*dx-x0;
						temp2=x*x+y*y+D*D;
						temp=(3.0*x*y)/pow(temp2,2.5);
						grav[i][j]+=temp*G*M*1E12;					//E
					}
				}
			}
		}
		break;
	case FORWARD_Vxz:
		{
			for (int k = 0; k < (int)forwardinfo.model.global_vec.size(); k++)
			{
				double R=forwardinfo.model.global_vec[k].Radius;
				double x0=forwardinfo.model.global_vec[k].Center[0];
				double y0=forwardinfo.model.global_vec[k].Center[1];
				double D=fabs(forwardinfo.model.global_vec[k].Center[2]-forwardheight);
				double density=forwardinfo.model.global_vec[k].Density;
				double M=4.0/3.0*PI*pow(R,3.0)*density;
				double x,y;
				double temp,temp2;
				for (int i=0;i< number_y;i++)
				{
					y=ymin+i*dy-y0;
					for (int j=0;j<number_x;j++)
					{
						x=xmin+j*dx-x0;
						temp2=x*x+y*y+D*D;
						temp=(-3.0*x*D)/pow(temp2,2.5);
						grav[i][j]+=temp*G*M*1E12;					//E
					}
				}
			}
		}
		break;
	case FORWARD_Vyz:
		{
			for (int k = 0; k < (int)forwardinfo.model.global_vec.size(); k++)
			{
				double R=forwardinfo.model.global_vec[k].Radius;
				double x0=forwardinfo.model.global_vec[k].Center[0];
				double y0=forwardinfo.model.global_vec[k].Center[1];
				double D=fabs(forwardinfo.model.global_vec[k].Center[2]-forwardheight);
				double density=forwardinfo.model.global_vec[k].Density;
				double M=4.0/3.0*PI*pow(R,3.0)*density;
				double x,y;
				double temp,temp2;
				for (int i=0;i< number_y;i++)
				{
					y=ymin+i*dy-y0;
					for (int j=0;j<number_x;j++)
					{
						x=xmin+j*dx-x0;
						temp2=x*x+y*y+D*D;
						temp=(-3.0*D*y)/pow(temp2,2.5);
						grav[i][j]+=temp*G*M*1E12;					//E
					}
				}
			}
		}
		break;
	default:
		MessageBox(NULL,_T("这是球体重力正演函数，请输入正确的type"),_T("错误提示"),MB_OK);
		return 0;
	}

	//获取数据最大值最小值
	GetGrdMinMax(grav,forwardinfo.model.grddatainfo);
	return 1;
}
int _3DGlobal_mag(double** mag,RegularGeometry3DForward& forwardinfo,int type)
{
	//首先将grav赋值为0，后面累加
	int number_x=forwardinfo.model.grddatainfo.GetNumber_x();
	int number_y=forwardinfo.model.grddatainfo.GetNumber_y();
	
	Assign_Array2(mag,number_y,number_x,0);

	//坐标范围
	double xmin=forwardinfo.model.grddatainfo.m_AxisBounds[0];
	double xmax=forwardinfo.model.grddatainfo.m_AxisBounds[1];
	double ymin=forwardinfo.model.grddatainfo.m_AxisBounds[2];
	double ymax=forwardinfo.model.grddatainfo.m_AxisBounds[3];
	double dx=forwardinfo.model.grddatainfo.m_Dx;
	double dy=forwardinfo.model.grddatainfo.m_Dy;
	double forwardheight=forwardinfo.model.grddatainfo.m_Height_data;
	switch (type)
	{
	case FORWARD_Za:
		{
			for (int k = 0; k < (int)forwardinfo.model.global_vec.size(); k++)
			{
				double R=forwardinfo.model.global_vec[k].Radius;
				double x0=forwardinfo.model.global_vec[k].Center[0];
				double y0=forwardinfo.model.global_vec[k].Center[1];
				double D=fabs(forwardinfo.model.global_vec[k].Center[2]-forwardheight);

				double angle_I=forwardinfo.model.global_vec[k].Angle_I;						//磁倾角（度）
				double angle_Apie=forwardinfo.model.global_vec[k].Angle_D;					//磁偏角（度）
				double CHL=forwardinfo.model.global_vec[k].CiHuaLv;							//磁化率
				double EarthMag=forwardinfo.EarthMag;
				double M=CHL*EarthMag;
				double m=4.0/3.0*PI*pow(R,3.0)*M;											//磁化强度（A/m）//这里没有除以U0是跟系数抵消掉了，

				//参数换算
				angle_I=angle_I/180.0*PI;												//换算为弧度
				angle_Apie=angle_Apie/180.0*PI;											//换算为弧度
				double x,y, temp,temp2,temp11,temp22,temp33;
				for (int i=0;i< number_y;i++)
				{
					x=ymin+i*dy-y0;
					for (int j=0;j<number_x;j++)
					{
						y=xmin+j*dx-x0;

						temp2=pow(x*x+y*y+D*D,-2.5);
						temp11=(2.0*D*D-x*x-y*y)*sin(angle_I);
						temp22=-3.0*D*x*cos(angle_I)*cos(angle_Apie);
						temp33=-3.0*D*y*cos(angle_I)*sin(angle_Apie);
						temp=temp2*(temp11+temp22+temp33);
						mag[i][j]+=m*temp/4.0/PI;					//nT
					}
				}
			}
		}
		break;
	case FORWARD_Hax:
		{
			for (int k = 0; k < (int)forwardinfo.model.global_vec.size(); k++)
			{
				double R=forwardinfo.model.global_vec[k].Radius;
				double x0=forwardinfo.model.global_vec[k].Center[0];
				double y0=forwardinfo.model.global_vec[k].Center[1];
				double D=fabs(forwardinfo.model.global_vec[k].Center[2]-forwardheight);

				double angle_I=forwardinfo.model.global_vec[k].Angle_I;						//磁倾角（度）
				double angle_Apie=forwardinfo.model.global_vec[k].Angle_D;					//磁偏角（度）
				double CHL=forwardinfo.model.global_vec[k].CiHuaLv;							//磁化率
				double EarthMag=forwardinfo.EarthMag;
				double M=CHL*EarthMag;
				double m=4.0/3.0*PI*pow(R,3.0)*M;											//磁化强度（A/m）//这里没有除以U0是跟系数抵消掉了，
				angle_I=angle_I/180.0*PI;												//换算为弧度
				angle_Apie=angle_Apie/180.0*PI;											//换算为弧度
				double x,y, temp,temp2,temp11,temp22,temp33;
				for (int i=0;i< number_y;i++)
				{
					x=ymin+i*dy-y0;
					for (int j=0;j<number_x;j++)
					{
						y=xmin+j*dx-x0;

						temp2=pow(x*x+y*y+D*D,-2.5);
						temp11=(2.0*x*x-y*y-D*D)*cos(angle_I)*cos(angle_Apie);
						temp22=-3.0*D*x*sin(angle_I);
						temp33=3.0*x*y*cos(angle_I)*sin(angle_Apie);
						temp=temp2*(temp11+temp22+temp33);
						mag[i][j]+=m*temp/4.0/PI;					//nT
					}
				}
			}
		}
		break;
	case FORWARD_Hay:
		{
			for (int k = 0; k < (int)forwardinfo.model.global_vec.size(); k++)
			{
				double R=forwardinfo.model.global_vec[k].Radius;
				double x0=forwardinfo.model.global_vec[k].Center[0];
				double y0=forwardinfo.model.global_vec[k].Center[1];
				double D=fabs(forwardinfo.model.global_vec[k].Center[2]-forwardheight);

				double angle_I=forwardinfo.model.global_vec[k].Angle_I;						//磁倾角（度）
				double angle_Apie=forwardinfo.model.global_vec[k].Angle_D;					//磁偏角（度）
				double CHL=forwardinfo.model.global_vec[k].CiHuaLv;							//磁化率
				double EarthMag=forwardinfo.EarthMag;
				double M=CHL*EarthMag;
				double m=4.0/3.0*PI*pow(R,3.0)*M;											//磁化强度（A/m）//这里没有除以U0是跟系数抵消掉了，
				//参数换算
				angle_I=angle_I/180.0*PI;												//换算为弧度
				angle_Apie=angle_Apie/180.0*PI;											//换算为弧度
				double x,y, temp,temp2,temp11,temp22,temp33;
				for (int i=0;i< number_y;i++)
				{
					x=ymin+i*dy-y0;
					for (int j=0;j<number_x;j++)
					{
						y=xmin+j*dx-x0;

						temp2=pow(x*x+y*y+D*D,-2.5);
						temp11=(2.0*y*y-x*x-D*D)*cos(angle_I)*sin(angle_Apie);
						temp22=-3.0*D*y*sin(angle_I);
						temp33=3.0*x*y*cos(angle_I)*cos(angle_Apie);
						temp=temp2*(temp11+temp22+temp33);
						mag[i][j]+=m*temp/4.0/PI;					//nT
					}
				}
			}
		}
		break;
	default:
		MessageBox(NULL,_T("这是球体磁力正演函数，请输入正确的type"),_T("错误提示"),MB_OK);
		return 0;
	}

	//获取数据最大值最小值
	GetGrdMinMax(mag,forwardinfo.model.grddatainfo);
	return 1;
}
int _3DGlobal_magT(double** mag,RegularGeometry3DForward& forwardinfo,int type)
{
	//首先将grav赋值为0，后面累加
	int number_x=forwardinfo.model.grddatainfo.GetNumber_x();
	int number_y=forwardinfo.model.grddatainfo.GetNumber_y();
	double** Hax=CreateArray2(number_y,number_x);
	double** Hay=CreateArray2(number_y,number_x);
	double** Za=CreateArray2(number_y,number_x);
	Assign_Array2(Hax,number_y,number_x,0);
	Assign_Array2(Hay,number_y,number_x,0);
	Assign_Array2(Za,number_y,number_x,0);
	Assign_Array2(mag,number_y,number_x,0);
	double I=forwardinfo.EarthAngle_I/180.0*PI;
	double D=forwardinfo.EarthAngle_D/180.0*PI;

	//计算Hax
	_3DGlobal_mag(Hax,forwardinfo,FORWARD_Hax);
	//计算Hay
	_3DGlobal_mag(Hay,forwardinfo,FORWARD_Hay);
	//计算Za
	_3DGlobal_mag(Za,forwardinfo,FORWARD_Za);
	switch (type)
	{
	case FORWARD_Ta:
		{
			//计算Ta
			for (int i = 0; i < number_y; i++)
			{
				for (int j = 0; j < number_x; j++)
				{
					mag[i][j]=Hax[i][j]*cos(I)*cos(D)+Hay[i][j]*cos(I)*sin(D)+Za[i][j]*sin(I);
				}
			}
		}
		break;
	case FORWARD_Module:
		{
			//计算模量
			for (int i = 0; i < number_y; i++)
			{
				for (int j = 0; j < number_x; j++)
				{
					mag[i][j]=sqrt(Hax[i][j]*Hax[i][j]+Hay[i][j]*Hay[i][j]+Za[i][j]*Za[i][j]);
				}
			}
		}
		break;
	default:
		MessageBox(NULL,_T("这是球体Ta和模量正演，输入正确的类型参数"),_T("出错提示"),MB_OK);
		return 0;
	}
	
	GetGrdMinMax(mag,forwardinfo.model.grddatainfo);
	//销毁二维数组
	DeleteArray2(Hax,number_y,number_x);
	DeleteArray2(Hay,number_y,number_x);
	DeleteArray2(Za,number_y,number_x);
	return 1;
}

int _3DHorizontalCylinder_magT(double** mag,RegularGeometry3DForward& forwardinfo,int type)
{
	//首先将grav赋值为0，后面累加
	int number_x=forwardinfo.model.grddatainfo.GetNumber_x();
	int number_y=forwardinfo.model.grddatainfo.GetNumber_y();
	double** Hax=CreateArray2(number_y,number_x);
	double** Hay=CreateArray2(number_y,number_x);
	double** Za=CreateArray2(number_y,number_x);
	Assign_Array2(Hax,number_y,number_x,0);
	Assign_Array2(Hay,number_y,number_x,0);
	Assign_Array2(Za,number_y,number_x,0);
	Assign_Array2(mag,number_y,number_x,0);
	double I=forwardinfo.EarthAngle_I/180.0*PI;
	double D=forwardinfo.EarthAngle_D/180.0*PI;

	//计算Hax
	_3DHorizontalCylinder_mag(Hax,forwardinfo,FORWARD_Hax);
	//计算Hay
	_3DHorizontalCylinder_mag(Hay,forwardinfo,FORWARD_Hay);
	//计算Za
	_3DHorizontalCylinder_mag(Za,forwardinfo,FORWARD_Za);
	switch (type)
	{
	case FORWARD_Ta:
		{
			//计算Ta
			for (int i = 0; i < number_y; i++)
			{
				for (int j = 0; j < number_x; j++)
				{
					mag[i][j]=Hax[i][j]*cos(I)*cos(D)+Hay[i][j]*cos(I)*sin(D)+Za[i][j]*sin(I);
				}
			}
		}
		break;
	case FORWARD_Module:
		{
			//计算模量
			for (int i = 0; i < number_y; i++)
			{
				for (int j = 0; j < number_x; j++)
				{
					mag[i][j]=sqrt(Hax[i][j]*Hax[i][j]+Hay[i][j]*Hay[i][j]+Za[i][j]*Za[i][j]);
				}
			}
		}
		break;
	default:
		MessageBox(NULL, _T("这是Ta和模量正演，输入正确的类型参数"), _T("出错提示"), MB_OK);
		return 0;
	}
	
	GetGrdMinMax(mag,forwardinfo.model.grddatainfo);
	//销毁二维数组
	DeleteArray2(Hax,number_y,number_x);
	DeleteArray2(Hay,number_y,number_x);
	DeleteArray2(Za,number_y,number_x);

	return 0;
}

int _3DHorizontalCylinder_mag(double** mag,RegularGeometry3DForward& forwardinfo,int type)
{
	//首先将mag赋值为0，后面累加
	int number_x=forwardinfo.model.grddatainfo.GetNumber_x();
	int number_y=forwardinfo.model.grddatainfo.GetNumber_y();
	Assign_Array2(mag,number_y,number_x,0);
	
	//坐标范围
	double xmin=forwardinfo.model.grddatainfo.m_AxisBounds[0];
	double xmax=forwardinfo.model.grddatainfo.m_AxisBounds[1];
	double ymin=forwardinfo.model.grddatainfo.m_AxisBounds[2];
	double ymax=forwardinfo.model.grddatainfo.m_AxisBounds[3];
	double dx=forwardinfo.model.grddatainfo.m_Dx;
	double dy=forwardinfo.model.grddatainfo.m_Dy;
	double xorigin,yorigin,x,y,L1,L2,temp1,temp2;
	double Center[3],Length;
	switch (type)
	{
	case FORWARD_Hay:
		{
			double TX,TY,TZ;
			double xx,DD;
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double angle_cylinder=atan2((forwardinfo.model.cylinder_vec[i].Pos1[1]-forwardinfo.model.cylinder_vec[i].Pos2[1]),(forwardinfo.model.cylinder_vec[i].Pos1[0]-forwardinfo.model.cylinder_vec[i].Pos2[0]));
				if (angle_cylinder<0)//atan2得出的是一个向量与x轴正方向的夹角，我们只需要这个直线与x轴正方向的夹角即可，因此需要把复制转换为正的，比如-135换算为45度
				{
					angle_cylinder+=PI;
				}
				angle_cylinder=PI/2.0-fabs(angle_cylinder);
				double Trans[2][2];//坐标变换矩阵
				Trans[0][0]=cos(angle_cylinder);Trans[0][1]=-sin(angle_cylinder);
				Trans[1][0]=sin(angle_cylinder);Trans[1][1]=cos(angle_cylinder);
				double xyOrigin[2],xyTrans[2];//前者是原始坐标的列向量，后者是坐标变换后的列向量
				//计算圆柱体中心点坐标和长度
				Length=Distance(forwardinfo.model.cylinder_vec[i].Pos1,forwardinfo.model.cylinder_vec[i].Pos2,2);
				Center[0]=(forwardinfo.model.cylinder_vec[i].Pos1[0]+forwardinfo.model.cylinder_vec[i].Pos2[0])/2.0;//圆柱体原始中心点位置
				Center[1]=(forwardinfo.model.cylinder_vec[i].Pos1[1]+forwardinfo.model.cylinder_vec[i].Pos2[1])/2.0;
				Center[2]=(forwardinfo.model.cylinder_vec[i].Pos1[2]+forwardinfo.model.cylinder_vec[i].Pos2[2])/2.0;
				//将中心点坐标变换到圆柱体坐标系
				Mat_Multiply(Trans[0],Center,xyTrans,2,2);
				double x0=xyTrans[0];
				double y0=xyTrans[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r=forwardinfo.model.cylinder_vec[i].Radius;
				double L=Length/2.0;
				double D=fabs(Center[2])+forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				/*double x0=Center[0];
				double y0=Center[1];*/
				double S=PI*r*r;
				double CHL=forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag=forwardinfo.EarthMag;
				double M=S*CHL*EarthMag;											//与系数消去了U0
				double angle_I=forwardinfo.model.cylinder_vec[i].Angle_I/180.0*PI;//化为弧度
				double angle_D=forwardinfo.model.cylinder_vec[i].Angle_D/180.0*PI;
				double Mx=M*cos(angle_I)*sin(angle_D);
				double My=M*cos(angle_I)*cos(angle_D);
				double Mz=M*sin(angle_I);
				for (int j = 0; j < number_y; j++)
				{
					//y=ymin+dy*j-y0;
					yorigin=ymin+dy*j;
					for (int k = 0; k < number_x; k++)
					{
						//x=xmin+k*dx-x0;
						xorigin=xmin+k*dx;
						//将计算点变换到圆柱体坐标系中
						xyOrigin[0]=xorigin;xyOrigin[1]=yorigin;
						Mat_Multiply(Trans[0],xyOrigin,xyTrans,2,2);
						x=xyTrans[0]-x0;y=xyTrans[1]-y0;

						L1=L-y;L2=L+y;
						temp1=sqrt(x*x+D*D+L1*L1);
						temp2=sqrt(x*x+D*D+L2*L2);
						xx=x*x;DD=D*D;
						TX=1.0/(xx+DD)*(Sign(L1)*abs(L1)/temp1*((xx-DD)/(xx+DD)-xx/pow(temp1,2.0))+Sign(L2)*abs(L2)/temp2*((xx-DD)/(xx+DD)-xx/pow(temp2,2.0)));
						TY=x*(1.0/pow(temp1,3.0)-1.0/pow(temp2,3.0));
						TZ=-D*x/(xx+DD)*(L1/temp1*((2.0/(xx+DD))+1.0/pow(temp1,2.0))+L2/temp2*(2.0/(xx+DD)+1.0/pow(temp2,2.0)));
						mag[j][k]+=1.0/4.0/PI*(Mx*TX+My*TY+Mz*TZ);
					}
				}
			}
		}
		break;
	case FORWARD_Hax:
		{
			double TX,TY,TZ;
			double temp13,temp23;
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double angle_cylinder=atan2((forwardinfo.model.cylinder_vec[i].Pos1[1]-forwardinfo.model.cylinder_vec[i].Pos2[1]),(forwardinfo.model.cylinder_vec[i].Pos1[0]-forwardinfo.model.cylinder_vec[i].Pos2[0]));
				if (angle_cylinder<0)//atan2得出的是一个向量与x轴正方向的夹角，我们只需要这个直线与x轴正方向的夹角即可，因此需要把复制转换为正的，比如-135换算为45度
				{
					angle_cylinder+=PI;
				}
				angle_cylinder=PI/2.0-fabs(angle_cylinder);
				double Trans[2][2];//坐标变换矩阵
				Trans[0][0]=cos(angle_cylinder);Trans[0][1]=-sin(angle_cylinder);
				Trans[1][0]=sin(angle_cylinder);Trans[1][1]=cos(angle_cylinder);
				double xyOrigin[2],xyTrans[2];//前者是原始坐标的列向量，后者是坐标变换后的列向量
				//计算圆柱体中心点坐标和长度
				Length=Distance(forwardinfo.model.cylinder_vec[i].Pos1,forwardinfo.model.cylinder_vec[i].Pos2,2);
				Center[0]=(forwardinfo.model.cylinder_vec[i].Pos1[0]+forwardinfo.model.cylinder_vec[i].Pos2[0])/2.0;//圆柱体原始中心点位置
				Center[1]=(forwardinfo.model.cylinder_vec[i].Pos1[1]+forwardinfo.model.cylinder_vec[i].Pos2[1])/2.0;
				Center[2]=(forwardinfo.model.cylinder_vec[i].Pos1[2]+forwardinfo.model.cylinder_vec[i].Pos2[2])/2.0;
				//将中心点坐标变换到圆柱体坐标系
				Mat_Multiply(Trans[0],Center,xyTrans,2,2);
				double x0=xyTrans[0];
				double y0=xyTrans[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r=forwardinfo.model.cylinder_vec[i].Radius;
				double L=Length/2.0;//L这里的L是半长度
				double D=fabs(Center[2])+forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				//double x0=Center[0];
				//double y0=Center[1];
				double S=PI*r*r;
				double CHL=forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag=forwardinfo.EarthMag;
				double M=S*CHL*EarthMag;
				double angle_I=forwardinfo.model.cylinder_vec[i].Angle_I/180.0*PI;//化为弧度
				double angle_D=forwardinfo.model.cylinder_vec[i].Angle_D/180.0*PI;
				double Mx=M*cos(angle_I)*sin(angle_D);
				double My=M*cos(angle_I)*cos(angle_D);
				double Mz=M*sin(angle_I);
				for (int j = 0; j < number_y; j++)
				{
					//y=ymin+dy*j-y0;
					yorigin=ymin+dy*j;
					for (int k = 0; k < number_x; k++)
					{
						//x=xmin+k*dx-x0;
						xorigin=xmin+k*dx;
						//将计算点变换到圆柱体坐标系中
						xyOrigin[0]=xorigin;xyOrigin[1]=yorigin;
						Mat_Multiply(Trans[0],xyOrigin,xyTrans,2,2);
						x=xyTrans[0]-x0;y=xyTrans[1]-y0;

						L1=L-y;L2=L+y;
						temp1=sqrt(x*x+D*D+L1*L1);
						temp2=sqrt(x*x+D*D+L2*L2);
						temp13=pow(temp1,3.0);
						temp23=pow(temp2,3.0);
						TX=x*(1.0/temp13-1.0/temp23);
						TY=-(L1/temp13+L2/temp23);
						TZ=-D*(1.0/temp13-1.0/temp23);
						mag[j][k]+=1.0/4.0/PI*(Mx*TX+My*TY+Mz*TZ);
					}
				}
			}
		}
		break;
	case FORWARD_Za:
		{
			double TX,TY,TZ;
			double temp12,temp22,xx,DD,temp13,temp23;
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double angle_cylinder=atan2((forwardinfo.model.cylinder_vec[i].Pos1[1]-forwardinfo.model.cylinder_vec[i].Pos2[1]),(forwardinfo.model.cylinder_vec[i].Pos1[0]-forwardinfo.model.cylinder_vec[i].Pos2[0]));
				if (angle_cylinder<0)//atan2得出的是一个向量与x轴正方向的夹角，我们只需要这个直线与x轴正方向的夹角即可，因此需要把复制转换为正的，比如-135换算为45度
				{
					angle_cylinder+=PI;
				}
				angle_cylinder=PI/2.0-fabs(angle_cylinder);
				double Trans[2][2];//坐标变换矩阵
				Trans[0][0]=cos(angle_cylinder);Trans[0][1]=-sin(angle_cylinder);
				Trans[1][0]=sin(angle_cylinder);Trans[1][1]=cos(angle_cylinder);
				double xyOrigin[2],xyTrans[2];//前者是原始坐标的列向量，后者是坐标变换后的列向量
				//计算圆柱体中心点坐标
				Length=Distance(forwardinfo.model.cylinder_vec[i].Pos1,forwardinfo.model.cylinder_vec[i].Pos2,2);
				Center[0]=(forwardinfo.model.cylinder_vec[i].Pos1[0]+forwardinfo.model.cylinder_vec[i].Pos2[0])/2.0;//圆柱体原始中心点位置
				Center[1]=(forwardinfo.model.cylinder_vec[i].Pos1[1]+forwardinfo.model.cylinder_vec[i].Pos2[1])/2.0;
				Center[2]=(forwardinfo.model.cylinder_vec[i].Pos1[2]+forwardinfo.model.cylinder_vec[i].Pos2[2])/2.0;
				//将中心点坐标变换到圆柱体坐标系
				Mat_Multiply(Trans[0],Center,xyTrans,2,2);
				double x0=xyTrans[0];
				double y0=xyTrans[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r=forwardinfo.model.cylinder_vec[i].Radius;
				double L=Length/2.0;//L这里的L是半长度
				double D=fabs(Center[2])+forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
//				double x0=Center[0];
//				double y0=Center[1];
				double S=PI*r*r;
				double CHL=forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag=forwardinfo.EarthMag;
				double M=S*CHL*EarthMag;
				double angle_I=forwardinfo.model.cylinder_vec[i].Angle_I/180.0*PI;//化为弧度
				double angle_D=forwardinfo.model.cylinder_vec[i].Angle_D/180.0*PI;
				double Mx=M*cos(angle_I)*sin(angle_D);
				double My=M*cos(angle_I)*cos(angle_D);
				double Mz=M*sin(angle_I);
				for (int j = 0; j < number_y; j++)
				{
					//y=ymin+dy*j-y0;
					yorigin=ymin+dy*j;
					for (int k = 0; k < number_x; k++)
					{
						//x=xmin+k*dx-x0;
						xorigin=xmin+k*dx;
						//将计算点变换到圆柱体坐标系中
						xyOrigin[0]=xorigin;xyOrigin[1]=yorigin;
						Mat_Multiply(Trans[0],xyOrigin,xyTrans,2,2);
						x=xyTrans[0]-x0;y=xyTrans[1]-y0;

						L1=L-y;L2=L+y;
						temp1=sqrt(x*x+D*D+L1*L1);
						temp2=sqrt(x*x+D*D+L2*L2);
						temp12=pow(temp1,2.0);
						temp22=pow(temp2,2.0);
						temp13=temp12*temp1;
						temp23=temp22*temp2;
						xx=x*x;DD=D*D;
						TX=-D*x/(xx+DD)*(L1/temp1*(2.0/(xx+DD)+1.0/temp12)+L2/temp2*(2.0/(xx+DD)+1.0/temp22));
						TY=-D*(1.0/temp13-1.0/temp23);
						TZ=1/(xx+DD)*(L1/temp1*((DD-xx)/(xx+DD)+DD/temp12)+L2/temp2*((DD-xx)/(xx+DD)+DD/temp22));
						mag[j][k]+=1.0/4.0/PI*(Mx*TX+My*TY+Mz*TZ);
					}
				}
			}
		}
		break;
	default:
		MessageBox(NULL, _T("这是磁三分量正演函数，请输入正确的type"), _T("错误提示"), MB_OK);
		return 0;

	}
	GetGrdMinMax(mag,forwardinfo.model.grddatainfo);
	return 1;
}

int _3DHorizontalCylinder_grav(double** grav,RegularGeometry3DForward& forwardinfo,int type)
{
	//首先将grav赋值为0，后面累加
	int number_x=forwardinfo.model.grddatainfo.GetNumber_x();
	int number_y=forwardinfo.model.grddatainfo.GetNumber_y();
	
	Assign_Array2(grav,number_y,number_x,0);

	//坐标范围
	double xmin=forwardinfo.model.grddatainfo.m_AxisBounds[0];
	double xmax=forwardinfo.model.grddatainfo.m_AxisBounds[1];
	double ymin=forwardinfo.model.grddatainfo.m_AxisBounds[2];
	double ymax=forwardinfo.model.grddatainfo.m_AxisBounds[3];
	double dx=forwardinfo.model.grddatainfo.m_Dx;
	double dy=forwardinfo.model.grddatainfo.m_Dy;
	double yorigin,xorigin,x,y,L1,L2,B1,B2,temp1,temp2,v1,v2;
	double Center[3],Length;
	switch (type)
	{
	case FORWARD_V:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double angle_cylinder=atan2((forwardinfo.model.cylinder_vec[i].Pos1[1]-forwardinfo.model.cylinder_vec[i].Pos2[1]),(forwardinfo.model.cylinder_vec[i].Pos1[0]-forwardinfo.model.cylinder_vec[i].Pos2[0]));
				if (angle_cylinder<0)//atan2得出的是一个向量与x轴正方向的夹角，我们只需要这个直线与x轴正方向的夹角即可，因此需要把复制转换为正的，比如-135换算为45度
				{
					angle_cylinder+=PI;
				}
				angle_cylinder=PI/2.0-fabs(angle_cylinder);
				double Trans[2][2];//坐标变换矩阵
				Trans[0][0]=cos(angle_cylinder);Trans[0][1]=-sin(angle_cylinder);
				Trans[1][0]=sin(angle_cylinder);Trans[1][1]=cos(angle_cylinder);
				double xyOrigin[2],xyTrans[2];//前者是原始坐标的列向量，后者是坐标变换后的列向量
				//计算圆柱体中心点坐标
				Length=Distance(forwardinfo.model.cylinder_vec[i].Pos1,forwardinfo.model.cylinder_vec[i].Pos2,2);
				Center[0]=(forwardinfo.model.cylinder_vec[i].Pos1[0]+forwardinfo.model.cylinder_vec[i].Pos2[0])/2.0;//圆柱体原始中心点位置
				Center[1]=(forwardinfo.model.cylinder_vec[i].Pos1[1]+forwardinfo.model.cylinder_vec[i].Pos2[1])/2.0;
				Center[2]=(forwardinfo.model.cylinder_vec[i].Pos1[2]+forwardinfo.model.cylinder_vec[i].Pos2[2])/2.0;
				//将中心点坐标变换到圆柱体坐标系
				Mat_Multiply(Trans[0],Center,xyTrans,2,2);
				double x0=xyTrans[0];
				double y0=xyTrans[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r=forwardinfo.model.cylinder_vec[i].Radius;
				double L=Length/2.0;//L这里的L是半长度
				double D=fabs(Center[2])+forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double den=forwardinfo.model.cylinder_vec[i].Density;
				double S=PI*r*r;
				double LineDen=den*S;
				//double x0=Center[0];
				//double y0=Center[1];

				for (int j = 0; j < number_y; j++)
				{
					//y=ymin+dy*j-y0;
					yorigin=ymin+dy*j;
					for (int k = 0; k < number_x; k++)
					{
						//x=xmin+k*dx-x0;
						xorigin=xmin+k*dx;
						//将计算点变换到圆柱体坐标系中
						xyOrigin[0]=xorigin;xyOrigin[1]=yorigin;
						Mat_Multiply(Trans[0],xyOrigin,xyTrans,2,2);
						x=xyTrans[0]-x0;y=xyTrans[1]-y0;

						L1=L-y;L2=L+y;
						B1=L1*L1;B2=L2*L2;
						temp1=sqrt(x*x+D*D+L1*L1);
						temp2=sqrt(x*x+D*D+L2*L2);
						v1=Sign(L1)*log((temp1-sqrt(B1))/(temp1+sqrt(B1)));
						v2=Sign(L2)*log((temp2-sqrt(B2))/(temp2+sqrt(B2)));
						grav[j][k]+=-1E8*G*LineDen*(v1+v2)/2.0;
				
					}
				}
			}
		}
		break;
	case FORWARD_Vx:
		{
			double vx1,vx2;
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double angle_cylinder=atan2((forwardinfo.model.cylinder_vec[i].Pos1[1]-forwardinfo.model.cylinder_vec[i].Pos2[1]),(forwardinfo.model.cylinder_vec[i].Pos1[0]-forwardinfo.model.cylinder_vec[i].Pos2[0]));
				if (angle_cylinder<0)//atan2得出的是一个向量与x轴正方向的夹角，我们只需要这个直线与x轴正方向的夹角即可，因此需要把复制转换为正的，比如-135换算为45度
				{
					angle_cylinder+=PI;
				}
				angle_cylinder=PI/2.0-fabs(angle_cylinder);
				double Trans[2][2];//坐标变换矩阵
				Trans[0][0]=cos(angle_cylinder);Trans[0][1]=-sin(angle_cylinder);
				Trans[1][0]=sin(angle_cylinder);Trans[1][1]=cos(angle_cylinder);
				double xyOrigin[2],xyTrans[2];//前者是原始坐标的列向量，后者是坐标变换后的列向量
				//计算圆柱体中心点坐标
				Length=Distance(forwardinfo.model.cylinder_vec[i].Pos1,forwardinfo.model.cylinder_vec[i].Pos2,2);
				Center[0]=(forwardinfo.model.cylinder_vec[i].Pos1[0]+forwardinfo.model.cylinder_vec[i].Pos2[0])/2.0;//圆柱体原始中心点位置
				Center[1]=(forwardinfo.model.cylinder_vec[i].Pos1[1]+forwardinfo.model.cylinder_vec[i].Pos2[1])/2.0;
				Center[2]=(forwardinfo.model.cylinder_vec[i].Pos1[2]+forwardinfo.model.cylinder_vec[i].Pos2[2])/2.0;
				//将中心点坐标变换到圆柱体坐标系
				Mat_Multiply(Trans[0],Center,xyTrans,2,2);
				double x0=xyTrans[0];
				double y0=xyTrans[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r=forwardinfo.model.cylinder_vec[i].Radius;
				double L=Length/2.0;//L这里的L是半长度
				double den=forwardinfo.model.cylinder_vec[i].Density;
				double S=PI*r*r;
				double LineDen=den*S;
				double D=fabs(Center[2])+forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
//				double x0=Center[0];
//				double y0=Center[1];

				for (int j = 0; j < number_y; j++)
				{
					//y=ymin+dy*j-y0;
					yorigin=ymin+dy*j;
					for (int k = 0; k < number_x; k++)
					{
						//x=xmin+k*dx-x0;
						xorigin=xmin+k*dx;
						//将计算点变换到圆柱体坐标系中
						xyOrigin[0]=xorigin;xyOrigin[1]=yorigin;
						Mat_Multiply(Trans[0],xyOrigin,xyTrans,2,2);
						x=xyTrans[0]-x0;y=xyTrans[1]-y0;

						L1=L-y;L2=L+y;
						B1=L1*L1;B2=L2*L2;
						temp1=sqrt(x*x+D*D+L1*L1);
						temp2=sqrt(x*x+D*D+L2*L2);
						vx1=Sign(L1)*fabs(L1)/temp1;
						vx2=Sign(L2)*fabs(L2)/temp2;
						grav[j][k]+=-1E8*G*LineDen*x/(x*x+D*D)*(vx1+vx2);
				
					}
				}
			}
		}
		break;
	case FORWARD_Vy:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double angle_cylinder=atan2((forwardinfo.model.cylinder_vec[i].Pos1[1]-forwardinfo.model.cylinder_vec[i].Pos2[1]),(forwardinfo.model.cylinder_vec[i].Pos1[0]-forwardinfo.model.cylinder_vec[i].Pos2[0]));
				if (angle_cylinder<0)//atan2得出的是一个向量与x轴正方向的夹角，我们只需要这个直线与x轴正方向的夹角即可，因此需要把复制转换为正的，比如-135换算为45度
				{
					angle_cylinder+=PI;
				}
				angle_cylinder=PI/2.0-fabs(angle_cylinder);
				double Trans[2][2];//坐标变换矩阵
				Trans[0][0]=cos(angle_cylinder);Trans[0][1]=-sin(angle_cylinder);
				Trans[1][0]=sin(angle_cylinder);Trans[1][1]=cos(angle_cylinder);
				double xyOrigin[2],xyTrans[2];//前者是原始坐标的列向量，后者是坐标变换后的列向量
				//计算圆柱体中心点坐标
				Length=Distance(forwardinfo.model.cylinder_vec[i].Pos1,forwardinfo.model.cylinder_vec[i].Pos2,2);
				Center[0]=(forwardinfo.model.cylinder_vec[i].Pos1[0]+forwardinfo.model.cylinder_vec[i].Pos2[0])/2.0;//圆柱体原始中心点位置
				Center[1]=(forwardinfo.model.cylinder_vec[i].Pos1[1]+forwardinfo.model.cylinder_vec[i].Pos2[1])/2.0;
				Center[2]=(forwardinfo.model.cylinder_vec[i].Pos1[2]+forwardinfo.model.cylinder_vec[i].Pos2[2])/2.0;
				//将中心点坐标变换到圆柱体坐标系
				Mat_Multiply(Trans[0],Center,xyTrans,2,2);
				double x0=xyTrans[0];
				double y0=xyTrans[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r=forwardinfo.model.cylinder_vec[i].Radius;
				double L=Length/2.0;//L这里的L是半长度
				double den=forwardinfo.model.cylinder_vec[i].Density;
				double S=PI*r*r;
				double LineDen=den*S;
				double D=fabs(Center[2])+forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
//				double x0=Center[0];
//				double y0=Center[1];

				for (int j = 0; j < number_y; j++)
				{
					//y=ymin+dy*j-y0;
					yorigin=ymin+dy*j;
					for (int k = 0; k < number_x; k++)
					{
						//x=xmin+k*dx-x0;
						xorigin=xmin+k*dx;
						//将计算点变换到圆柱体坐标系中
						xyOrigin[0]=xorigin;xyOrigin[1]=yorigin;
						Mat_Multiply(Trans[0],xyOrigin,xyTrans,2,2);
						x=xyTrans[0]-x0;y=xyTrans[1]-y0;

						L1=L-y;L2=L+y;
						B1=L1*L1;B2=L2*L2;
						temp1=sqrt(x*x+D*D+L1*L1);
						temp2=sqrt(x*x+D*D+L2*L2);
						grav[j][k]+=-1E8*G*LineDen*(1.0/temp1-1.0/temp2);
				
					}
				}
			}
		}
		break;
	case FORWARD_Vz:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double angle_cylinder=atan2((forwardinfo.model.cylinder_vec[i].Pos1[1]-forwardinfo.model.cylinder_vec[i].Pos2[1]),(forwardinfo.model.cylinder_vec[i].Pos1[0]-forwardinfo.model.cylinder_vec[i].Pos2[0]));
				if (angle_cylinder<0)//atan2得出的是一个向量与x轴正方向的夹角，我们只需要这个直线与x轴正方向的夹角即可，因此需要把复制转换为正的，比如-135换算为45度
				{
					angle_cylinder+=PI;
				}
				angle_cylinder=PI/2.0-fabs(angle_cylinder);
				double Trans[2][2];//坐标变换矩阵
				Trans[0][0]=cos(angle_cylinder);Trans[0][1]=-sin(angle_cylinder);
				Trans[1][0]=sin(angle_cylinder);Trans[1][1]=cos(angle_cylinder);
				double xyOrigin[2],xyTrans[2];//前者是原始坐标的列向量，后者是坐标变换后的列向量
				//计算圆柱体中心点坐标
				Length=Distance(forwardinfo.model.cylinder_vec[i].Pos1,forwardinfo.model.cylinder_vec[i].Pos2,2);
				Center[0]=(forwardinfo.model.cylinder_vec[i].Pos1[0]+forwardinfo.model.cylinder_vec[i].Pos2[0])/2.0;//圆柱体原始中心点位置
				Center[1]=(forwardinfo.model.cylinder_vec[i].Pos1[1]+forwardinfo.model.cylinder_vec[i].Pos2[1])/2.0;
				Center[2]=(forwardinfo.model.cylinder_vec[i].Pos1[2]+forwardinfo.model.cylinder_vec[i].Pos2[2])/2.0;
				//将中心点坐标变换到圆柱体坐标系
				Mat_Multiply(Trans[0],Center,xyTrans,2,2);
				double x0=xyTrans[0];
				double y0=xyTrans[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r=forwardinfo.model.cylinder_vec[i].Radius;
				double L=Length/2.0;//L这里的L是半长度
				double den=forwardinfo.model.cylinder_vec[i].Density;
				double S=PI*r*r;
				double LineDen=den*S;
				double D=fabs(Center[2])+forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
//				double x0=Center[0];
//				double y0=Center[1];

				for (int j = 0; j < number_y; j++)
				{
					//y=ymin+dy*j-y0;
					yorigin=ymin+dy*j;
					for (int k = 0; k < number_x; k++)
					{
						//x=xmin+k*dx-x0;
						xorigin=xmin+k*dx;
						//将计算点变换到圆柱体坐标系中
						xyOrigin[0]=xorigin;xyOrigin[1]=yorigin;
						Mat_Multiply(Trans[0],xyOrigin,xyTrans,2,2);
						x=xyTrans[0]-x0;y=xyTrans[1]-y0;

						L1=L-y;L2=L+y;
						B1=L1*L1;B2=L2*L2;
						temp1=sqrt(x*x+D*D+B1);
						temp2=sqrt(x*x+D*D+B2);
						
						grav[j][k]+=1E8*G*LineDen*D*(L1/temp1+L2/temp2)/(x*x+D*D);
				
					}
				}
			}
		}
		break;
	case FORWARD_Vxx:
		{
			double vxx1,vxx2,vx1,vx2;
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double angle_cylinder=atan2((forwardinfo.model.cylinder_vec[i].Pos1[1]-forwardinfo.model.cylinder_vec[i].Pos2[1]),(forwardinfo.model.cylinder_vec[i].Pos1[0]-forwardinfo.model.cylinder_vec[i].Pos2[0]));
				if (angle_cylinder<0)//atan2得出的是一个向量与x轴正方向的夹角，我们只需要这个直线与x轴正方向的夹角即可，因此需要把复制转换为正的，比如-135换算为45度
				{
					angle_cylinder+=PI;
				}
				angle_cylinder=PI/2.0-fabs(angle_cylinder);
				double Trans[2][2];//坐标变换矩阵
				Trans[0][0]=cos(angle_cylinder);Trans[0][1]=-sin(angle_cylinder);
				Trans[1][0]=sin(angle_cylinder);Trans[1][1]=cos(angle_cylinder);
				double xyOrigin[2],xyTrans[2];//前者是原始坐标的列向量，后者是坐标变换后的列向量
				//计算圆柱体中心点坐标
				Length=Distance(forwardinfo.model.cylinder_vec[i].Pos1,forwardinfo.model.cylinder_vec[i].Pos2,2);
				Center[0]=(forwardinfo.model.cylinder_vec[i].Pos1[0]+forwardinfo.model.cylinder_vec[i].Pos2[0])/2.0;//圆柱体原始中心点位置
				Center[1]=(forwardinfo.model.cylinder_vec[i].Pos1[1]+forwardinfo.model.cylinder_vec[i].Pos2[1])/2.0;
				Center[2]=(forwardinfo.model.cylinder_vec[i].Pos1[2]+forwardinfo.model.cylinder_vec[i].Pos2[2])/2.0;
				//将中心点坐标变换到圆柱体坐标系
				Mat_Multiply(Trans[0],Center,xyTrans,2,2);
				double x0=xyTrans[0];
				double y0=xyTrans[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r=forwardinfo.model.cylinder_vec[i].Radius;
				double L=Length/2.0;//L这里的L是半长度
				double den=forwardinfo.model.cylinder_vec[i].Density;
				double S=PI*r*r;
				double LineDen=den*S;
				double D=fabs(Center[2])+forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
//				double x0=Center[0];
//				double y0=Center[1];

				for (int j = 0; j < number_y; j++)
				{
					//y=ymin+dy*j-y0;
					yorigin=ymin+dy*j;
					for (int k = 0; k < number_x; k++)
					{
						//x=xmin+k*dx-x0;
						xorigin=xmin+k*dx;
						//将计算点变换到圆柱体坐标系中
						xyOrigin[0]=xorigin;xyOrigin[1]=yorigin;
						Mat_Multiply(Trans[0],xyOrigin,xyTrans,2,2);
						x=xyTrans[0]-x0;y=xyTrans[1]-y0;

						L1=L-y;L2=L+y;
						B1=L1*L1;B2=L2*L2;
						temp1=sqrt(x*x+D*D+B1);
						temp2=sqrt(x*x+D*D+B2);
						vx1=Sign(L1)*fabs(L1)/temp1;
						vx2=Sign(L2)*fabs(L2)/temp2;
						vxx1=(L1)/pow(temp1,3.0);
						vxx2=(L2)/pow(temp2,3.0);
						grav[j][k]+=1E12*G*LineDen*((x*x-D*D)/pow((x*x+D*D),2.0)*(vx1+vx2)-x*x/(x*x+D*D)*(vxx1+vxx2));
				
					}
				}
			}
		}
		break;
	case FORWARD_Vxy:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double angle_cylinder=atan2((forwardinfo.model.cylinder_vec[i].Pos1[1]-forwardinfo.model.cylinder_vec[i].Pos2[1]),(forwardinfo.model.cylinder_vec[i].Pos1[0]-forwardinfo.model.cylinder_vec[i].Pos2[0]));
				if (angle_cylinder<0)//atan2得出的是一个向量与x轴正方向的夹角，我们只需要这个直线与x轴正方向的夹角即可，因此需要把复制转换为正的，比如-135换算为45度
				{
					angle_cylinder+=PI;
				}
				angle_cylinder=PI/2.0-fabs(angle_cylinder);
				double Trans[2][2];//坐标变换矩阵
				Trans[0][0]=cos(angle_cylinder);Trans[0][1]=-sin(angle_cylinder);
				Trans[1][0]=sin(angle_cylinder);Trans[1][1]=cos(angle_cylinder);
				double xyOrigin[2],xyTrans[2];//前者是原始坐标的列向量，后者是坐标变换后的列向量
				//计算圆柱体中心点坐标
				Length=Distance(forwardinfo.model.cylinder_vec[i].Pos1,forwardinfo.model.cylinder_vec[i].Pos2,2);
				Center[0]=(forwardinfo.model.cylinder_vec[i].Pos1[0]+forwardinfo.model.cylinder_vec[i].Pos2[0])/2.0;//圆柱体原始中心点位置
				Center[1]=(forwardinfo.model.cylinder_vec[i].Pos1[1]+forwardinfo.model.cylinder_vec[i].Pos2[1])/2.0;
				Center[2]=(forwardinfo.model.cylinder_vec[i].Pos1[2]+forwardinfo.model.cylinder_vec[i].Pos2[2])/2.0;
				//将中心点坐标变换到圆柱体坐标系
				Mat_Multiply(Trans[0],Center,xyTrans,2,2);
				double x0=xyTrans[0];
				double y0=xyTrans[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r=forwardinfo.model.cylinder_vec[i].Radius;
				double L=Length/2.0;//L这里的L是半长度
				double den=forwardinfo.model.cylinder_vec[i].Density;
				double S=PI*r*r;
				double LineDen=den*S;
				double D=fabs(Center[2])+forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
//				double x0=Center[0];
//				double y0=Center[1];

				for (int j = 0; j < number_y; j++)
				{
					//y=ymin+dy*j-y0;
					yorigin=ymin+dy*j;
					for (int k = 0; k < number_x; k++)
					{
						//x=xmin+k*dx-x0;
						xorigin=xmin+k*dx;
						//将计算点变换到圆柱体坐标系中
						xyOrigin[0]=xorigin;xyOrigin[1]=yorigin;
						Mat_Multiply(Trans[0],xyOrigin,xyTrans,2,2);
						x=xyTrans[0]-x0;y=xyTrans[1]-y0;

						L1=L-y;L2=L+y;
						B1=L1*L1;B2=L2*L2;
						temp1=sqrt(x*x+D*D+B1);
						temp2=sqrt(x*x+D*D+B2);
						grav[j][k]+=1E12*G*LineDen*(x/pow(temp1,3.0)-x/pow(temp2,3.0));
				
					}
				}
			}
		}
		break;
	case FORWARD_Vxz:
		{
			double vxz1,vxz2;
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double angle_cylinder=atan2((forwardinfo.model.cylinder_vec[i].Pos1[1]-forwardinfo.model.cylinder_vec[i].Pos2[1]),(forwardinfo.model.cylinder_vec[i].Pos1[0]-forwardinfo.model.cylinder_vec[i].Pos2[0]));
				if (angle_cylinder<0)//atan2得出的是一个向量与x轴正方向的夹角，我们只需要这个直线与x轴正方向的夹角即可，因此需要把复制转换为正的，比如-135换算为45度
				{
					angle_cylinder+=PI;
				}
				angle_cylinder=PI/2.0-fabs(angle_cylinder);
				double Trans[2][2];//坐标变换矩阵
				Trans[0][0]=cos(angle_cylinder);Trans[0][1]=-sin(angle_cylinder);
				Trans[1][0]=sin(angle_cylinder);Trans[1][1]=cos(angle_cylinder);
				double xyOrigin[2],xyTrans[2];//前者是原始坐标的列向量，后者是坐标变换后的列向量
				//计算圆柱体中心点坐标
				Length=Distance(forwardinfo.model.cylinder_vec[i].Pos1,forwardinfo.model.cylinder_vec[i].Pos2,2);
				Center[0]=(forwardinfo.model.cylinder_vec[i].Pos1[0]+forwardinfo.model.cylinder_vec[i].Pos2[0])/2.0;//圆柱体原始中心点位置
				Center[1]=(forwardinfo.model.cylinder_vec[i].Pos1[1]+forwardinfo.model.cylinder_vec[i].Pos2[1])/2.0;
				Center[2]=(forwardinfo.model.cylinder_vec[i].Pos1[2]+forwardinfo.model.cylinder_vec[i].Pos2[2])/2.0;
				//将中心点坐标变换到圆柱体坐标系
				Mat_Multiply(Trans[0],Center,xyTrans,2,2);
				double x0=xyTrans[0];
				double y0=xyTrans[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r=forwardinfo.model.cylinder_vec[i].Radius;
				double L=Length/2.0;//L这里的L是半长度
				double den=forwardinfo.model.cylinder_vec[i].Density;
				double S=PI*r*r;
				double LineDen=den*S;
				double D=fabs(Center[2])+forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
//				double x0=Center[0];
//				double y0=Center[1];

				for (int j = 0; j < number_y; j++)
				{
					//y=ymin+dy*j-y0;
					yorigin=ymin+dy*j;
					for (int k = 0; k < number_x; k++)
					{
						//x=xmin+k*dx-x0;
						xorigin=xmin+k*dx;
						//将计算点变换到圆柱体坐标系中
						xyOrigin[0]=xorigin;xyOrigin[1]=yorigin;
						Mat_Multiply(Trans[0],xyOrigin,xyTrans,2,2);
						x=xyTrans[0]-x0;y=xyTrans[1]-y0;

						L1=L-y;L2=L+y;
						B1=L1*L1;B2=L2*L2;
						temp1=sqrt(x*x+D*D+B1);
						temp2=sqrt(x*x+D*D+B2);
						vxz1=(L1/(temp1)+L2/(temp2))*2*x/pow((x*x+D*D),2.0);
						vxz2=(L1/pow(temp1,3.0)+L2/pow(temp2,3.0))*x/(x*x+D*D);
						grav[j][k]+=-1E12*G*LineDen*D*(vxz1+vxz2);
				
					}
				}
			}
		}
		break;
	case FORWARD_Vyy:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double angle_cylinder=atan2((forwardinfo.model.cylinder_vec[i].Pos1[1]-forwardinfo.model.cylinder_vec[i].Pos2[1]),(forwardinfo.model.cylinder_vec[i].Pos1[0]-forwardinfo.model.cylinder_vec[i].Pos2[0]));
				if (angle_cylinder<0)//atan2得出的是一个向量与x轴正方向的夹角，我们只需要这个直线与x轴正方向的夹角即可，因此需要把复制转换为正的，比如-135换算为45度
				{
					angle_cylinder+=PI;
				}
				angle_cylinder=PI/2.0-fabs(angle_cylinder);
				double Trans[2][2];//坐标变换矩阵
				Trans[0][0]=cos(angle_cylinder);Trans[0][1]=-sin(angle_cylinder);
				Trans[1][0]=sin(angle_cylinder);Trans[1][1]=cos(angle_cylinder);
				double xyOrigin[2],xyTrans[2];//前者是原始坐标的列向量，后者是坐标变换后的列向量
				//计算圆柱体中心点坐标
				Length=Distance(forwardinfo.model.cylinder_vec[i].Pos1,forwardinfo.model.cylinder_vec[i].Pos2,2);
				Center[0]=(forwardinfo.model.cylinder_vec[i].Pos1[0]+forwardinfo.model.cylinder_vec[i].Pos2[0])/2.0;//圆柱体原始中心点位置
				Center[1]=(forwardinfo.model.cylinder_vec[i].Pos1[1]+forwardinfo.model.cylinder_vec[i].Pos2[1])/2.0;
				Center[2]=(forwardinfo.model.cylinder_vec[i].Pos1[2]+forwardinfo.model.cylinder_vec[i].Pos2[2])/2.0;
				//将中心点坐标变换到圆柱体坐标系
				Mat_Multiply(Trans[0],Center,xyTrans,2,2);
				double x0=xyTrans[0];
				double y0=xyTrans[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r=forwardinfo.model.cylinder_vec[i].Radius;
				double L=Length/2.0;//L这里的L是半长度
				double den=forwardinfo.model.cylinder_vec[i].Density;
				double S=PI*r*r;
				double LineDen=den*S;
				double D=fabs(Center[2])+forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
//				double x0=Center[0];
//				double y0=Center[1];

				for (int j = 0; j < number_y; j++)
				{
					//y=ymin+dy*j-y0;
					yorigin=ymin+dy*j;
					for (int k = 0; k < number_x; k++)
					{
						//x=xmin+k*dx-x0;
						xorigin=xmin+k*dx;
						//将计算点变换到圆柱体坐标系中
						xyOrigin[0]=xorigin;xyOrigin[1]=yorigin;
						Mat_Multiply(Trans[0],xyOrigin,xyTrans,2,2);
						x=xyTrans[0]-x0;y=xyTrans[1]-y0;

						L1=L-y;L2=L+y;
						B1=L1*L1;B2=L2*L2;
						temp1=sqrt(x*x+D*D+B1);
						temp2=sqrt(x*x+D*D+B2);
						grav[j][k]+=1E12*G*LineDen*(-L1/pow(temp1,3.0)-L2/pow(temp2,3.0));
				
					}
				}
			}
		}
		break;
	case FORWARD_Vyz:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double angle_cylinder=atan2((forwardinfo.model.cylinder_vec[i].Pos1[1]-forwardinfo.model.cylinder_vec[i].Pos2[1]),(forwardinfo.model.cylinder_vec[i].Pos1[0]-forwardinfo.model.cylinder_vec[i].Pos2[0]));
				if (angle_cylinder<0)//atan2得出的是一个向量与x轴正方向的夹角，我们只需要这个直线与x轴正方向的夹角即可，因此需要把复制转换为正的，比如-135换算为45度
				{
					angle_cylinder+=PI;
				}
				angle_cylinder=PI/2.0-fabs(angle_cylinder);
				double Trans[2][2];//坐标变换矩阵
				Trans[0][0]=cos(angle_cylinder);Trans[0][1]=-sin(angle_cylinder);
				Trans[1][0]=sin(angle_cylinder);Trans[1][1]=cos(angle_cylinder);
				double xyOrigin[2],xyTrans[2];//前者是原始坐标的列向量，后者是坐标变换后的列向量
				//计算圆柱体中心点坐标
				Length=Distance(forwardinfo.model.cylinder_vec[i].Pos1,forwardinfo.model.cylinder_vec[i].Pos2,2);
				Center[0]=(forwardinfo.model.cylinder_vec[i].Pos1[0]+forwardinfo.model.cylinder_vec[i].Pos2[0])/2.0;//圆柱体原始中心点位置
				Center[1]=(forwardinfo.model.cylinder_vec[i].Pos1[1]+forwardinfo.model.cylinder_vec[i].Pos2[1])/2.0;
				Center[2]=(forwardinfo.model.cylinder_vec[i].Pos1[2]+forwardinfo.model.cylinder_vec[i].Pos2[2])/2.0;
				//将中心点坐标变换到圆柱体坐标系
				Mat_Multiply(Trans[0],Center,xyTrans,2,2);
				double x0=xyTrans[0];
				double y0=xyTrans[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r=forwardinfo.model.cylinder_vec[i].Radius;
				double L=Length/2.0;//L这里的L是半长度
				double den=forwardinfo.model.cylinder_vec[i].Density;
				double S=PI*r*r;
				double LineDen=den*S;
				double D=fabs(Center[2])+forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
//				double x0=Center[0];
//				double y0=Center[1];

				for (int j = 0; j < number_y; j++)
				{
					//y=ymin+dy*j-y0;
					yorigin=ymin+dy*j;
					for (int k = 0; k < number_x; k++)
					{
						//x=xmin+k*dx-x0;
						xorigin=xmin+k*dx;
						//将计算点变换到圆柱体坐标系中
						xyOrigin[0]=xorigin;xyOrigin[1]=yorigin;
						Mat_Multiply(Trans[0],xyOrigin,xyTrans,2,2);
						x=xyTrans[0]-x0;y=xyTrans[1]-y0;

						L1=L-y;L2=L+y;
						B1=L1*L1;B2=L2*L2;
						temp1=sqrt(x*x+D*D+B1);
						temp2=sqrt(x*x+D*D+B2);
						grav[j][k]+=-1E12*G*LineDen*D*(1.0/pow(temp1,3.0)-1.0/pow(temp2,3.0));
				
					}
				}
			}
		}
		break;
	case FORWARD_Vzz:
		{
			double vzz1,vzz2;
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double angle_cylinder=atan2((forwardinfo.model.cylinder_vec[i].Pos1[1]-forwardinfo.model.cylinder_vec[i].Pos2[1]),(forwardinfo.model.cylinder_vec[i].Pos1[0]-forwardinfo.model.cylinder_vec[i].Pos2[0]));
				if (angle_cylinder<0)//atan2得出的是一个向量与x轴正方向的夹角，我们只需要这个直线与x轴正方向的夹角即可，因此需要把复制转换为正的，比如-135换算为45度
				{
					angle_cylinder+=PI;
				}
				angle_cylinder=PI/2.0-fabs(angle_cylinder);
				double Trans[2][2];//坐标变换矩阵
				Trans[0][0]=cos(angle_cylinder);Trans[0][1]=-sin(angle_cylinder);
				Trans[1][0]=sin(angle_cylinder);Trans[1][1]=cos(angle_cylinder);
				double xyOrigin[2],xyTrans[2];//前者是原始坐标的列向量，后者是坐标变换后的列向量
				//计算圆柱体中心点坐标
				Length=Distance(forwardinfo.model.cylinder_vec[i].Pos1,forwardinfo.model.cylinder_vec[i].Pos2,2);
				Center[0]=(forwardinfo.model.cylinder_vec[i].Pos1[0]+forwardinfo.model.cylinder_vec[i].Pos2[0])/2.0;//圆柱体原始中心点位置
				Center[1]=(forwardinfo.model.cylinder_vec[i].Pos1[1]+forwardinfo.model.cylinder_vec[i].Pos2[1])/2.0;
				Center[2]=(forwardinfo.model.cylinder_vec[i].Pos1[2]+forwardinfo.model.cylinder_vec[i].Pos2[2])/2.0;
				//将中心点坐标变换到圆柱体坐标系
				Mat_Multiply(Trans[0],Center,xyTrans,2,2);
				double x0=xyTrans[0];
				double y0=xyTrans[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r=forwardinfo.model.cylinder_vec[i].Radius;
				double L=Length/2.0;//L这里的L是半长度
				double den=forwardinfo.model.cylinder_vec[i].Density;
				double S=PI*r*r;
				double LineDen=den*S;
				double D=fabs(Center[2])+forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
//				double x0=Center[0];
//				double y0=Center[1];

				for (int j = 0; j < number_y; j++)
				{
					//y=ymin+dy*j-y0;
					yorigin=ymin+dy*j;
					for (int k = 0; k < number_x; k++)
					{
						//x=xmin+k*dx-x0;
						xorigin=xmin+k*dx;
						//将计算点变换到圆柱体坐标系中
						xyOrigin[0]=xorigin;xyOrigin[1]=yorigin;
						Mat_Multiply(Trans[0],xyOrigin,xyTrans,2,2);
						x=xyTrans[0]-x0;y=xyTrans[1]-y0;

						L1=L-y;L2=L+y;
						B1=L1*L1;B2=L2*L2;
						temp1=sqrt(x*x+D*D+B1);
						temp2=sqrt(x*x+D*D+B2);
						vzz1=(x*x-D*D)/pow((x*x+D*D),2.0)*(L1/temp1+L2/temp2); 
						vzz2=(D*D)/(x*x+D*D)*(L1/pow(temp1,3.0)+L2/pow(temp2,3.0));
						grav[j][k]+=-1E12*G*LineDen*(vzz1-vzz2);				
					}
				}
			}
		}
		break;
	default:
		MessageBox(NULL, _T("正演分量类型不在范围之内"), _T("出错提示"), MB_OK);
		break;
	}
	
	GetGrdMinMax(grav,forwardinfo.model.grddatainfo);
	return 0;
}

int _3DCube_grav(double** grav,RegularGeometry3DForward& forwardinfo,int type)
{
	//首先将grav赋值为0，后面累加
	int number_x=forwardinfo.model.grddatainfo.AutoGetNumber_x();
	int number_y=forwardinfo.model.grddatainfo.AutoGetNumber_y();
	
	Assign_Array2(grav,number_y,number_x,0);
	
	//坐标范围
	double xmin=forwardinfo.model.grddatainfo.m_AxisBounds[0];
	double xmax=forwardinfo.model.grddatainfo.m_AxisBounds[1];
	double ymin=forwardinfo.model.grddatainfo.m_AxisBounds[2];
	double ymax=forwardinfo.model.grddatainfo.m_AxisBounds[3];
	double dx=forwardinfo.model.grddatainfo.m_Dx;
	double dy=forwardinfo.model.grddatainfo.m_Dy;

	
	switch (type)
	{
	case FORWARD_Vz:
	{
			for (int k = 0; k < (int)forwardinfo.model.cube_vec.size(); k++)
			{
				double* bounds=forwardinfo.model.cube_vec[k].bounds;
				double forwardheight=forwardinfo.model.grddatainfo.m_Height_data;
				double x1=bounds[0],x2=bounds[1],y1=bounds[2],y2=bounds[3];
				double z1=fabs(bounds[5]-forwardheight),z2=fabs(bounds[4]-forwardheight);//公式中的z1和z2分别表示顶面和地面埋深，因此这里需要转换一下
				double density=forwardinfo.model.cube_vec[k].Density;
				//中间变量
				double tempx1,tempx2,tempy1,tempy2,G122,G121,G112,G111,G222,G221,G212,G211;
				double p222,p221,p212,p211,p122,p121,p112,p111;		//按x，y，z的顺序排列
				double p2,p1;
				p2=z2*z2;	//X方向计算式这是常数
				p1=z1*z1;
				double x,y;
				for (int i=0;i< number_y;i++)
				{
					for (int j=0;j<number_x;j++)
					{
						y=ymin+i*dy;
						x=xmin+j*dx;;
						tempx1=x1-x;		//计算积分限
						tempx2=x2-x;
		
						tempy1=y1-y;
						tempy2=y2-y;
						//tempy1=-20;tempy2=20;
			
						p222=sqrt(pow(tempx2,2)+pow(tempy2,2)+p2);		//距离
						p221=sqrt(pow(tempx2,2)+pow(tempy2,2)+p1);
						p212=sqrt(pow(tempx2,2)+pow(tempy1,2)+p2);
						p211=sqrt(pow(tempx2,2)+pow(tempy1,2)+p1);
			
						p122=sqrt(pow(tempx1,2)+pow(tempy2,2)+p2);		//距离
						p121=sqrt(pow(tempx1,2)+pow(tempy2,2)+p1);
						p112=sqrt(pow(tempx1,2)+pow(tempy1,2)+p2);
						p111=sqrt(pow(tempx1,2)+pow(tempy1,2)+p1);
			
						G222=tempx2*log(tempy2+p222)+tempy2*log(tempx2+p222)+z2*atan2((z2*p222),(tempx2*tempy2));		//积分项
						G221=tempx2*log(tempy2+p221)+tempy2*log(tempx2+p221)+z1*atan2((z1*p221),(tempx2*tempy2));
						G212=tempx2*log(tempy1+p212)+tempy1*log(tempx2+p212)+z2*atan2((z2*p212),(tempx2*tempy1));
						G211=tempx2*log(tempy1+p211)+tempy1*log(tempx2+p211)+z1*atan2((z1*p211),(tempx2*tempy1));
			
						G122=tempx1*log(tempy2+p122)+tempy2*log(tempx1+p122)+z2*atan2((z2*p122),(tempx1*tempy2));
						G121=tempx1*log(tempy2+p121)+tempy2*log(tempx1+p121)+z1*atan2((z1*p121),(tempx1*tempy2));
						G112=tempx1*log(tempy1+p112)+tempy1*log(tempx1+p112)+z2*atan2((z2*p112),(tempx1*tempy1));
						G111=tempx1*log(tempy1+p111)+tempy1*log(tempx1+p111)+z1*atan2((z1*p111),(tempx1*tempy1));
			
						grav[i][j]+=-G*1E8*density*(G222+G211+G121+G112-G221-G212-G122-G111);	//mGal
						
					}
				}
			}
			
		}
		break;
	case FORWARD_Vxx:
		{
			for (int k = 0; k < (int)forwardinfo.model.cube_vec.size(); k++)
			{
				double* bounds=forwardinfo.model.cube_vec[k].bounds;
				double forwardheight=forwardinfo.model.grddatainfo.m_Height_data;
				double x1=bounds[0],x2=bounds[1],y1=bounds[2],y2=bounds[3];
				double z1=fabs(bounds[5]-forwardheight),z2=fabs(bounds[4]-forwardheight);//公式中的z1和z2分别表示顶面和地面埋深，因此这里需要转换一下
				double density=forwardinfo.model.cube_vec[k].Density;
				//中间变量
				double tempx1,tempx2,tempy1,tempy2,G122,G121,G112,G111,G222,G221,G212,G211;
				double p222,p221,p212,p211,p122,p121,p112,p111;		//按x，y，z的顺序排列
				double p2,p1;
				p2=z2*z2;	//X方向计算式这是常数
				p1=z1*z1;
				double x,y;
				for (int i=0;i< number_y;i++)
				{
					for (int j=0;j<number_x;j++)
					{
						y=ymin+i*dy;
						x=xmin+j*dx;;
						tempx1=x1-x;		//计算积分限
						tempx2=x2-x;
		
						tempy1=y1-y;
						tempy2=y2-y;
						//tempy1=-20;tempy2=20;
			
						p222=sqrt(pow(tempx2,2)+pow(tempy2,2)+p2);		//距离
						p221=sqrt(pow(tempx2,2)+pow(tempy2,2)+p1);
						p212=sqrt(pow(tempx2,2)+pow(tempy1,2)+p2);
						p211=sqrt(pow(tempx2,2)+pow(tempy1,2)+p1);
			
						p122=sqrt(pow(tempx1,2)+pow(tempy2,2)+p2);		//距离
						p121=sqrt(pow(tempx1,2)+pow(tempy2,2)+p1);
						p112=sqrt(pow(tempx1,2)+pow(tempy1,2)+p2);
						p111=sqrt(pow(tempx1,2)+pow(tempy1,2)+p1);
			
						G222=atan2(tempx2*tempy2,pow(tempx2,2.0)+p222*z2+p2);		//积分项
						G221=atan2(tempx2*tempy2,pow(tempx2,2.0)+p221*z1+p1);
						G212=atan2(tempx2*tempy1,pow(tempx2,2.0)+p212*z2+p2);
						G211=atan2(tempx2*tempy1,pow(tempx2,2.0)+p211*z1+p1);

						G122=atan2(tempx1*tempy2,pow(tempx1,2.0)+p122*z2+p2);		//积分项
						G121=atan2(tempx1*tempy2,pow(tempx1,2.0)+p121*z1+p1);
						G112=atan2(tempx1*tempy1,pow(tempx1,2.0)+p112*z2+p2);
						G111=atan2(tempx1*tempy1,pow(tempx1,2.0)+p111*z1+p1);

						grav[i][j]+=G*1E12*density*(G222+G211+G121+G112-G221-G212-G122-G111);	//E
			
					}
				}
			}
		}
		break;
	case FORWARD_Vxy:
		{
			for (int k = 0; k < (int)forwardinfo.model.cube_vec.size(); k++)
			{
				double* bounds=forwardinfo.model.cube_vec[k].bounds;
				double forwardheight=forwardinfo.model.grddatainfo.m_Height_data;
				double x1=bounds[0],x2=bounds[1],y1=bounds[2],y2=bounds[3];
				double z1=fabs(bounds[5]-forwardheight),z2=fabs(bounds[4]-forwardheight);//公式中的z1和z2分别表示顶面和地面埋深，因此这里需要转换一下
				double density=forwardinfo.model.cube_vec[k].Density;
				//中间变量
				double tempx1,tempx2,tempy1,tempy2,G122,G121,G112,G111,G222,G221,G212,G211;
				double p222,p221,p212,p211,p122,p121,p112,p111;		//按x，y，z的顺序排列
				double p2,p1;
				p2=z2*z2;	//X方向计算式这是常数
				p1=z1*z1;
				double x,y;
				for (int i=0;i< number_y;i++)
				{
					for (int j=0;j<number_x;j++)
					{
						y=ymin+i*dy;
						x=xmin+j*dx;;
						tempx1=x1-x;		//计算积分限
						tempx2=x2-x;
		
						tempy1=y1-y;
						tempy2=y2-y;
						//tempy1=-20;tempy2=20;
			
						p222=sqrt(pow(tempx2,2)+pow(tempy2,2)+p2);		//距离
						p221=sqrt(pow(tempx2,2)+pow(tempy2,2)+p1);
						p212=sqrt(pow(tempx2,2)+pow(tempy1,2)+p2);
						p211=sqrt(pow(tempx2,2)+pow(tempy1,2)+p1);
			
						p122=sqrt(pow(tempx1,2)+pow(tempy2,2)+p2);		//距离
						p121=sqrt(pow(tempx1,2)+pow(tempy2,2)+p1);
						p112=sqrt(pow(tempx1,2)+pow(tempy1,2)+p2);
						p111=sqrt(pow(tempx1,2)+pow(tempy1,2)+p1);
			
						G222=log(z2+p222);		//积分项
						G221=log(z1+p221);
						G212=log(z2+p212);
						G211=log(z1+p211);

						G122=log(z2+p122);		//积分项
						G121=log(z1+p121);
						G112=log(z2+p112);
						G111=log(z1+p111);

						grav[i][j]+=G*1E12*density*(G222+G211+G121+G112-G221-G212-G122-G111);	//E
					}
				}
			}
		}
		break;
	case FORWARD_Vxz:
		{
			for (int k = 0; k < (int)forwardinfo.model.cube_vec.size(); k++)
			{
				double* bounds=forwardinfo.model.cube_vec[k].bounds;
				double forwardheight=forwardinfo.model.grddatainfo.m_Height_data;
				double x1=bounds[0],x2=bounds[1],y1=bounds[2],y2=bounds[3];
				double z1=fabs(bounds[5]-forwardheight),z2=fabs(bounds[4]-forwardheight);//公式中的z1和z2分别表示顶面和地面埋深，因此这里需要转换一下
				double density=forwardinfo.model.cube_vec[k].Density;
				//中间变量
				double tempx1,tempx2,tempy1,tempy2,G122,G121,G112,G111,G222,G221,G212,G211;
				double p222,p221,p212,p211,p122,p121,p112,p111;		//按x，y，z的顺序排列
				double p2,p1;
				p2=z2*z2;	//X方向计算式这是常数
				p1=z1*z1;
				double x,y;
				for (int i=0;i< number_y;i++)
				{
					for (int j=0;j<number_x;j++)
					{
						y=ymin+i*dy;
						x=xmin+j*dx;;
						tempx1=x1-x;		//计算积分限
						tempx2=x2-x;
		
						tempy1=y1-y;
						tempy2=y2-y;
						//tempy1=-20;tempy2=20;
			
						p222=sqrt(pow(tempx2,2)+pow(tempy2,2)+p2);		//距离
						p221=sqrt(pow(tempx2,2)+pow(tempy2,2)+p1);
						p212=sqrt(pow(tempx2,2)+pow(tempy1,2)+p2);
						p211=sqrt(pow(tempx2,2)+pow(tempy1,2)+p1);
			
						p122=sqrt(pow(tempx1,2)+pow(tempy2,2)+p2);		//距离
						p121=sqrt(pow(tempx1,2)+pow(tempy2,2)+p1);
						p112=sqrt(pow(tempx1,2)+pow(tempy1,2)+p2);
						p111=sqrt(pow(tempx1,2)+pow(tempy1,2)+p1);
			
						G222=log(tempy2+p222);		//积分项
						G221=log(tempy2+p221);
						G212=log(tempy1+p212);
						G211=log(tempy1+p211);

						G122=log(tempy2+p122);		//积分项
						G121=log(tempy2+p121);
						G112=log(tempy1+p112);
						G111=log(tempy1+p111);

						grav[i][j]+=G*1E12*density*(G222+G211+G121+G112-G221-G212-G122-G111);	//E
					}
				}
			}
		}
		break;
	case FORWARD_Vyy:
		{
			for (int k = 0; k < (int)forwardinfo.model.cube_vec.size(); k++)
			{
				double* bounds=forwardinfo.model.cube_vec[k].bounds;
				double forwardheight=forwardinfo.model.grddatainfo.m_Height_data;
				double x1=bounds[0],x2=bounds[1],y1=bounds[2],y2=bounds[3];
				double z1=fabs(bounds[5]-forwardheight),z2=fabs(bounds[4]-forwardheight);//公式中的z1和z2分别表示顶面和地面埋深，因此这里需要转换一下
				double density=forwardinfo.model.cube_vec[k].Density;
				//中间变量
				double tempx1,tempx2,tempy1,tempy2,G122,G121,G112,G111,G222,G221,G212,G211;
				double p222,p221,p212,p211,p122,p121,p112,p111;		//按x，y，z的顺序排列
				double p2,p1;
				p2=z2*z2;	//X方向计算式这是常数
				p1=z1*z1;
				double x,y;
				for (int i=0;i< number_y;i++)
				{
					for (int j=0;j<number_x;j++)
					{
						y=ymin+i*dy;
						x=xmin+j*dx;;
						tempx1=x1-x;		//计算积分限
						tempx2=x2-x;
		
						tempy1=y1-y;
						tempy2=y2-y;
						//tempy1=-20;tempy2=20;
			
						p222=sqrt(pow(tempx2,2)+pow(tempy2,2)+p2);		//距离
						p221=sqrt(pow(tempx2,2)+pow(tempy2,2)+p1);
						p212=sqrt(pow(tempx2,2)+pow(tempy1,2)+p2);
						p211=sqrt(pow(tempx2,2)+pow(tempy1,2)+p1);
			
						p122=sqrt(pow(tempx1,2)+pow(tempy2,2)+p2);		//距离
						p121=sqrt(pow(tempx1,2)+pow(tempy2,2)+p1);
						p112=sqrt(pow(tempx1,2)+pow(tempy1,2)+p2);
						p111=sqrt(pow(tempx1,2)+pow(tempy1,2)+p1);
			
						G222=atan2(tempx2*tempy2,pow(tempy2,2.0)+p222*z2+p2);		//积分项
						G221=atan2(tempx2*tempy2,pow(tempy2,2.0)+p221*z1+p1);
						G212=atan2(tempx2*tempy1,pow(tempy1,2.0)+p212*z2+p2);
						G211=atan2(tempx2*tempy1,pow(tempy1,2.0)+p211*z1+p1);

						G122=atan2(tempx1*tempy2,pow(tempy2,2.0)+p122*z2+p2);		//积分项
						G121=atan2(tempx1*tempy2,pow(tempy2,2.0)+p121*z1+p1);
						G112=atan2(tempx1*tempy1,pow(tempy1,2.0)+p112*z2+p2);
						G111=atan2(tempx1*tempy1,pow(tempy1,2.0)+p111*z1+p1);

						grav[i][j]+=G*1E12*density*(G222+G211+G121+G112-G221-G212-G122-G111);	//E
					}
				}
			}
		}
		break;
	case FORWARD_Vyz:
		{
			for (int k = 0; k < (int)forwardinfo.model.cube_vec.size(); k++)
			{
				double* bounds=forwardinfo.model.cube_vec[k].bounds;
				double forwardheight=forwardinfo.model.grddatainfo.m_Height_data;
				double x1=bounds[0],x2=bounds[1],y1=bounds[2],y2=bounds[3];
				double z1=fabs(bounds[5]-forwardheight),z2=fabs(bounds[4]-forwardheight);//公式中的z1和z2分别表示顶面和地面埋深，因此这里需要转换一下
				double density=forwardinfo.model.cube_vec[k].Density;
				//中间变量
				double tempx1,tempx2,tempy1,tempy2,G122,G121,G112,G111,G222,G221,G212,G211;
				double p222,p221,p212,p211,p122,p121,p112,p111;		//按x，y，z的顺序排列
				double p2,p1;
				p2=z2*z2;	//X方向计算式这是常数
				p1=z1*z1;
				double x,y;
				for (int i=0;i< number_y;i++)
				{
					for (int j=0;j<number_x;j++)
					{
						y=ymin+i*dy;
						x=xmin+j*dx;;
						tempx1=x1-x;		//计算积分限
						tempx2=x2-x;
		
						tempy1=y1-y;
						tempy2=y2-y;
						//tempy1=-20;tempy2=20;
			
						p222=sqrt(pow(tempx2,2)+pow(tempy2,2)+p2);		//距离
						p221=sqrt(pow(tempx2,2)+pow(tempy2,2)+p1);
						p212=sqrt(pow(tempx2,2)+pow(tempy1,2)+p2);
						p211=sqrt(pow(tempx2,2)+pow(tempy1,2)+p1);
			
						p122=sqrt(pow(tempx1,2)+pow(tempy2,2)+p2);		//距离
						p121=sqrt(pow(tempx1,2)+pow(tempy2,2)+p1);
						p112=sqrt(pow(tempx1,2)+pow(tempy1,2)+p2);
						p111=sqrt(pow(tempx1,2)+pow(tempy1,2)+p1);
			
						G222=log(tempx2+p222);		//积分项
						G221=log(tempx2+p221);
						G212=log(tempx2+p212);
						G211=log(tempx2+p211);

						G122=log(tempx1+p122);		//积分项
						G121=log(tempx1+p121);
						G112=log(tempx1+p112);
						G111=log(tempx1+p111);

						grav[i][j]+=G*1E12*density*(G222+G211+G121+G112-G221-G212-G122-G111);	//E
					}
				}
			}
		}
		break;
	case FORWARD_Vzz:
		{
			for (int k = 0; k < (int)forwardinfo.model.cube_vec.size(); k++)
			{
				double* bounds=forwardinfo.model.cube_vec[k].bounds;
				double forwardheight=forwardinfo.model.grddatainfo.m_Height_data;
				double x1=bounds[0],x2=bounds[1],y1=bounds[2],y2=bounds[3];
				double z1=fabs(bounds[5]-forwardheight),z2=fabs(bounds[4]-forwardheight);//公式中的z1和z2分别表示顶面和地面埋深，因此这里需要转换一下
				double density=forwardinfo.model.cube_vec[k].Density;
				//中间变量
				double tempx1,tempx2,tempy1,tempy2,G122,G121,G112,G111,G222,G221,G212,G211;
				double p222,p221,p212,p211,p122,p121,p112,p111;		//按x，y，z的顺序排列
				double p2,p1;
				p2=z2*z2;	//X方向计算式这是常数
				p1=z1*z1;
				double x,y;
				for (int i=0;i< number_y;i++)
				{
					for (int j=0;j<number_x;j++)
					{
						y=ymin+i*dy;
						x=xmin+j*dx;;
						tempx1=x1-x;		//计算积分限
						tempx2=x2-x;
		
						tempy1=y1-y;
						tempy2=y2-y;
						//tempy1=-20;tempy2=20;
			
						p222=sqrt(pow(tempx2,2)+pow(tempy2,2)+p2);		//距离
						p221=sqrt(pow(tempx2,2)+pow(tempy2,2)+p1);
						p212=sqrt(pow(tempx2,2)+pow(tempy1,2)+p2);
						p211=sqrt(pow(tempx2,2)+pow(tempy1,2)+p1);
			
						p122=sqrt(pow(tempx1,2)+pow(tempy2,2)+p2);		//距离
						p121=sqrt(pow(tempx1,2)+pow(tempy2,2)+p1);
						p112=sqrt(pow(tempx1,2)+pow(tempy1,2)+p2);
						p111=sqrt(pow(tempx1,2)+pow(tempy1,2)+p1);
			
						G222=atan2(tempx2*tempy2,p222*z2);		//积分项
						G221=atan2(tempx2*tempy2,p221*z1);
						G212=atan2(tempx2*tempy1,p212*z2);
						G211=atan2(tempx2*tempy1,p211*z1);

						G122=atan2(tempx1*tempy2,p122*z2);		//积分项
						G121=atan2(tempx1*tempy2,p121*z1);
						G112=atan2(tempx1*tempy1,p112*z2);
						G111=atan2(tempx1*tempy1,p111*z1);

						grav[i][j]+=-G*1E12*density*(G222+G211+G121+G112-G221-G212-G122-G111);	//E
					}
				}
			}
		}
		break;
	default:
		MessageBox(NULL, _T("正演分量类型不在范围之内"), _T("出错提示"), MB_OK);
		return 0;
	}
	
	GetGrdMinMax(grav,forwardinfo.model.grddatainfo);
	return 1;
}

int _3DCube_mag(double** mag,RegularGeometry3DForward& forwardinfo,int type)
{
	//首先将grav赋值为0，后面累加
	int number_x=forwardinfo.model.grddatainfo.GetNumber_x();
	int number_y=forwardinfo.model.grddatainfo.GetNumber_y();
	
	Assign_Array2(mag,number_y,number_x,0);

	//坐标范围
	double xmin=forwardinfo.model.grddatainfo.m_AxisBounds[0];
	double xmax=forwardinfo.model.grddatainfo.m_AxisBounds[1];
	double ymin=forwardinfo.model.grddatainfo.m_AxisBounds[2];
	double ymax=forwardinfo.model.grddatainfo.m_AxisBounds[3];
	double dx=forwardinfo.model.grddatainfo.m_Dx;
	double dy=forwardinfo.model.grddatainfo.m_Dy;
	switch (type)
	{
	case FORWARD_Hax:
		{
			for (int k = 0; k < (int)forwardinfo.model.cube_vec.size(); k++)
			{
				double* bounds=forwardinfo.model.cube_vec[k].bounds;
				double forwardheight=forwardinfo.model.grddatainfo.m_Height_data;
				double x1=bounds[2],x2=bounds[3],y1=bounds[0],y2=bounds[1];				//为了则surfer中绘制模型投影图不出错，模型参数里的y其实才是磁的北，应赋值给x变量
				double z1=fabs(bounds[5]-forwardheight),z2=fabs(bounds[4]-forwardheight);//公式中的z1和z2分别表示顶面和地面埋深，因此这里需要转换一下
				double angle_I=forwardinfo.model.cube_vec[k].Angle_I;						//磁倾角（度）
				double angle_Apie=forwardinfo.model.cube_vec[k].Angle_D;					//磁偏角（度）
				double CHL=forwardinfo.model.cube_vec[k].CiHuaLv;							//磁化率
				double EarthMag=forwardinfo.EarthMag;
				double m=CHL*EarthMag;;													//磁化强度（A/m）//这里没有除以U0是跟系数抵消掉了，
				//参数换算
				angle_I=angle_I/180.0*PI;												//换算为弧度
				angle_Apie=angle_Apie/180.0*PI;											//换算为弧度
				double mx=m*cos(angle_I)*cos(angle_Apie);
				double my=m*cos(angle_I)*sin(angle_Apie);
				double mz=m*sin(angle_I);
				//中间变量
				double tempx1,tempx2,tempy1,tempy2,G122,G121,G112,G111,G222,G221,G212,G211;
				double p222,p221,p212,p211,p122,p121,p112,p111;							//按x，y，z的顺序排列
				double p2,p1;
				p2=pow(z2,2);															//X方向计算时这是常数
				p1=pow(z1,2);
				double x,y;
				for (int i=0;i< number_y;i++)
				{
					x=ymin+i*dy;
					for (int j=0;j<number_x;j++)
					{
						y=xmin+j*dx;
				
						tempx1=x1-x;													//计算积分限
						tempx2=x2-x;
			
						tempy1=y1-y;
						tempy2=y2-y;
						//tempy1=-20;tempy2=20;
			
						p222=sqrt(pow(tempx2,2)+pow(tempy2,2)+p2);						//距离
						p221=sqrt(pow(tempx2,2)+pow(tempy2,2)+p1);
						p212=sqrt(pow(tempx2,2)+pow(tempy1,2)+p2);
						p211=sqrt(pow(tempx2,2)+pow(tempy1,2)+p1);
			
						p122=sqrt(pow(tempx1,2)+pow(tempy2,2)+p2);						//距离
						p121=sqrt(pow(tempx1,2)+pow(tempy2,2)+p1);
						p112=sqrt(pow(tempx1,2)+pow(tempy1,2)+p2);
						p111=sqrt(pow(tempx1,2)+pow(tempy1,2)+p1);
			
						G222=mx*atan2((x2-x)*(y2-y),pow((x2-x),2)+p222*(z2)+pow(z2,2))+my*log(p222+z2)+mz*log(p222+y2-y);		//积分项
						G221=mx*atan2((x2-x)*(y2-y),pow((x2-x),2)+p221*(z1)+pow(z1,2))+my*log(p221+z1)+mz*log(p221+y2-y);
						G212=mx*atan2((x2-x)*(y1-y),pow((x2-x),2)+p212*(z2)+pow(z2,2))+my*log(p212+z2)+mz*log(p212+y1-y);
						G211=mx*atan2((x2-x)*(y1-y),pow((x2-x),2)+p211*(z1)+pow(z1,2))+my*log(p211+z1)+mz*log(p211+y1-y);
			
						G122=mx*atan2((x1-x)*(y2-y),pow((x1-x),2)+p122*(z2)+pow(z2,2))+my*log(p122+z2)+mz*log(p122+y2-y);		//积分项
						G121=mx*atan2((x1-x)*(y2-y),pow((x1-x),2)+p121*(z1)+pow(z1,2))+my*log(p121+z1)+mz*log(p121+y2-y);
						G112=mx*atan2((x1-x)*(y1-y),pow((x1-x),2)+p112*(z2)+pow(z2,2))+my*log(p112+z2)+mz*log(p112+y1-y);
						G111=mx*atan2((x1-x)*(y1-y),pow((x1-x),2)+p111*(z1)+pow(z1,2))+my*log(p111+z1)+mz*log(p111+y1-y);
			
						mag[i][j]+=(G222+G211+G121+G112-G221-G212-G122-G111)/(4.0*PI);			//nT
					}
				}
			}
		}
		break;
	case FORWARD_Hay:
		{
			for (int k = 0; k < (int)forwardinfo.model.cube_vec.size(); k++)
			{
				double* bounds=forwardinfo.model.cube_vec[k].bounds;
				double forwardheight=forwardinfo.model.grddatainfo.m_Height_data;
				double x1=bounds[2],x2=bounds[3],y1=bounds[0],y2=bounds[1];				//为了则surfer中绘制模型投影图不出错，模型参数里的y其实才是磁的北，应赋值给x变量
				double z1=fabs(bounds[5]-forwardheight),z2=fabs(bounds[4]-forwardheight);//公式中的z1和z2分别表示顶面和地面埋深，因此这里需要转换一下
				double angle_I=forwardinfo.model.cube_vec[k].Angle_I;						//磁倾角（度）
				double angle_Apie=forwardinfo.model.cube_vec[k].Angle_D;					//磁偏角（度）
				double CHL=forwardinfo.model.cube_vec[k].CiHuaLv;							//磁化率
				double EarthMag=forwardinfo.EarthMag;
				double m=CHL*EarthMag;;													//磁化强度（A/m）
				//参数换算
				angle_I=angle_I/180.0*PI;												//换算为弧度
				angle_Apie=angle_Apie/180.0*PI;											//换算为弧度
				double mx=m*cos(angle_I)*cos(angle_Apie);
				double my=m*cos(angle_I)*sin(angle_Apie);
				double mz=m*sin(angle_I);
				//中间变量
				double tempx1,tempx2,tempy1,tempy2,G122,G121,G112,G111,G222,G221,G212,G211;
				double p222,p221,p212,p211,p122,p121,p112,p111;							//按x，y，z的顺序排列
				double p2,p1;
				p2=pow(z2,2);															//X方向计算时这是常数
				p1=pow(z1,2);
				double x,y;
				for (int i=0;i< number_y;i++)
				{
					x=ymin+i*dy;
					for (int j=0;j<number_x;j++)
					{
						y=xmin+j*dx;
				
						tempx1=x1-x;													//计算积分限
						tempx2=x2-x;
			
						tempy1=y1-y;
						tempy2=y2-y;
						//tempy1=-20;tempy2=20;
			
						p222=sqrt(pow(tempx2,2)+pow(tempy2,2)+p2);		//距离
						p221=sqrt(pow(tempx2,2)+pow(tempy2,2)+p1);
						p212=sqrt(pow(tempx2,2)+pow(tempy1,2)+p2);
						p211=sqrt(pow(tempx2,2)+pow(tempy1,2)+p1);
			
						p122=sqrt(pow(tempx1,2)+pow(tempy2,2)+p2);		//距离
						p121=sqrt(pow(tempx1,2)+pow(tempy2,2)+p1);
						p112=sqrt(pow(tempx1,2)+pow(tempy1,2)+p2);
						p111=sqrt(pow(tempx1,2)+pow(tempy1,2)+p1);
			
						G222=my*atan2((x2-x)*(y2-y),pow((y2-y),2)+p222*(z2)+pow(z2,2))+mx*log(p222+z2)+mz*log(p222+x2-x);		//积分项
						G221=my*atan2((x2-x)*(y2-y),pow((y2-y),2)+p221*(z1)+pow(z1,2))+mx*log(p221+z1)+mz*log(p221+x2-x);
						G212=my*atan2((x2-x)*(y1-y),pow((y1-y),2)+p212*(z2)+pow(z2,2))+mx*log(p212+z2)+mz*log(p212+x2-x);
						G211=my*atan2((x2-x)*(y1-y),pow((y1-y),2)+p211*(z1)+pow(z1,2))+mx*log(p211+z1)+mz*log(p211+x2-x);
			
						G122=my*atan2((x1-x)*(y2-y),pow((y2-y),2)+p122*(z2)+pow(z2,2))+mx*log(p122+z2)+mz*log(p122+x1-x);		//积分项
						G121=my*atan2((x1-x)*(y2-y),pow((y2-y),2)+p121*(z1)+pow(z1,2))+mx*log(p121+z1)+mz*log(p121+x1-x);
						G112=my*atan2((x1-x)*(y1-y),pow((y1-y),2)+p112*(z2)+pow(z2,2))+mx*log(p112+z2)+mz*log(p112+x1-x);
						G111=my*atan2((x1-x)*(y1-y),pow((y1-y),2)+p111*(z1)+pow(z1,2))+mx*log(p111+z1)+mz*log(p111+x1-x);
			
						mag[i][j]+=(G222+G211+G121+G112-G221-G212-G122-G111)/(4.0*PI);			//nT
					}
				}
			}
		}
		break;
	case FORWARD_Za:
		{
			for (int k = 0; k < (int)forwardinfo.model.cube_vec.size(); k++)
			{
				double* bounds=forwardinfo.model.cube_vec[k].bounds;
				double forwardheight=forwardinfo.model.grddatainfo.m_Height_data;
				double x1=bounds[2],x2=bounds[3],y1=bounds[0],y2=bounds[1];				//为了则surfer中绘制模型投影图不出错，模型参数里的y其实才是磁的北，应赋值给x变量
				double z1=fabs(bounds[5]-forwardheight),z2=fabs(bounds[4]-forwardheight);//公式中的z1和z2分别表示顶面和地面埋深，因此这里需要转换一下
				double angle_I=forwardinfo.model.cube_vec[k].Angle_I;						//磁倾角（度）
				double angle_Apie=forwardinfo.model.cube_vec[k].Angle_D;					//磁偏角（度）
				double CHL=forwardinfo.model.cube_vec[k].CiHuaLv;							//磁化率
				double EarthMag=forwardinfo.EarthMag;
				double m=CHL*EarthMag;;													//磁化强度（A/m）
				//参数换算
				angle_I=angle_I/180.0*PI;												//换算为弧度
				angle_Apie=angle_Apie/180.0*PI;											//换算为弧度
				double mx=m*cos(angle_I)*cos(angle_Apie);
				double my=m*cos(angle_I)*sin(angle_Apie);
				double mz=m*sin(angle_I);
				//中间变量
				double tempx1,tempx2,tempy1,tempy2,G122,G121,G112,G111,G222,G221,G212,G211;
				double p222,p221,p212,p211,p122,p121,p112,p111;							//按x，y，z的顺序排列
				double p2,p1;
				p2=pow(z2,2);															//X方向计算时这是常数
				p1=pow(z1,2);
				double x,y;
				for (int i=0;i< number_y;i++)
				{
					x=ymin+i*dy;
					for (int j=0;j<number_x;j++)
					{
						y=xmin+j*dx;
				
						tempx1=x1-x;		//计算积分限
						tempx2=x2-x;
			
						tempy1=y1-y;
						tempy2=y2-y;
						//tempy1=-20;tempy2=20;
			
						p222=sqrt(pow(tempx2,2)+pow(tempy2,2)+p2);		//距离
						p221=sqrt(pow(tempx2,2)+pow(tempy2,2)+p1);
						p212=sqrt(pow(tempx2,2)+pow(tempy1,2)+p2);
						p211=sqrt(pow(tempx2,2)+pow(tempy1,2)+p1);
			
						p122=sqrt(pow(tempx1,2)+pow(tempy2,2)+p2);		//距离
						p121=sqrt(pow(tempx1,2)+pow(tempy2,2)+p1);
						p112=sqrt(pow(tempx1,2)+pow(tempy1,2)+p2);
						p111=sqrt(pow(tempx1,2)+pow(tempy1,2)+p1);
			
						G222=mx*log(p222+y2-y)+my*log(p222+x2-x)-mz*atan2((x2-x)*(y2-y),p222*(z2));		//积分项
						G221=mx*log(p221+y2-y)+my*log(p221+x2-x)-mz*atan2((x2-x)*(y2-y),p221*(z1));
						G212=mx*log(p212+y1-y)+my*log(p212+x2-x)-mz*atan2((x2-x)*(y1-y),p212*(z2));
						G211=mx*log(p211+y1-y)+my*log(p211+x2-x)-mz*atan2((x2-x)*(y1-y),p211*(z1));
			
						G122=mx*log(p122+y2-y)+my*log(p122+x1-x)-mz*atan2((x1-x)*(y2-y),p122*(z2));		//积分项
						G121=mx*log(p121+y2-y)+my*log(p121+x1-x)-mz*atan2((x1-x)*(y2-y),p121*(z1));
						G112=mx*log(p112+y1-y)+my*log(p112+x1-x)-mz*atan2((x1-x)*(y1-y),p112*(z2));
						G111=mx*log(p111+y1-y)+my*log(p111+x1-x)-mz*atan2((x1-x)*(y1-y),p111*(z1));
			
						mag[i][j]+=(G222+G211+G121+G112-G221-G212-G122-G111)/(4.0*PI);			//nT
					}
				}
			}
		}
		break;
	default:
		MessageBox(NULL, _T("立方体磁正演分量类型不在范围之内"), _T("出错提示"), MB_OK);
		return 0;
	}
	GetGrdMinMax(mag,forwardinfo.model.grddatainfo);
	return 1;
}

int _3DCube_magT(double** mag,RegularGeometry3DForward& forwardinfo,int type)
{
	//首先将grav赋值为0，后面累加
	int number_x=forwardinfo.model.grddatainfo.GetNumber_x();
	int number_y=forwardinfo.model.grddatainfo.GetNumber_y();
	double** Hax=CreateArray2(number_y,number_x);
	double** Hay=CreateArray2(number_y,number_x);
	double** Za=CreateArray2(number_y,number_x);
	Assign_Array2(Hax,number_y,number_x,0);
	Assign_Array2(Hay,number_y,number_x,0);
	Assign_Array2(Za,number_y,number_x,0);
	Assign_Array2(mag,number_y,number_x,0);
	double I=forwardinfo.EarthAngle_I/180.0*PI;
	double D=forwardinfo.EarthAngle_D/180.0*PI;

	//计算Hax
	_3DCube_mag(Hax,forwardinfo,FORWARD_Hax);
	//计算Hay
	_3DCube_mag(Hay,forwardinfo,FORWARD_Hay);
	//计算Za
	_3DCube_mag(Za,forwardinfo,FORWARD_Za);
	switch (type)
	{
	case FORWARD_Ta:
		{
			//计算Ta
			for (int i = 0; i < number_y; i++)
			{
				for (int j = 0; j < number_x; j++)
				{
					mag[i][j]=Hax[i][j]*cos(I)*cos(D)+Hay[i][j]*cos(I)*sin(D)+Za[i][j]*sin(I);
				}
			}
		}
		break;
	case FORWARD_Module:
		{
			//计算模量
			for (int i = 0; i < number_y; i++)
			{
				for (int j = 0; j < number_x; j++)
				{
					mag[i][j]=sqrt(Hax[i][j]*Hax[i][j]+Hay[i][j]*Hay[i][j]+Za[i][j]*Za[i][j]);
				}
			}
		}
		break;
	default:
		MessageBox(NULL, _T("这是立方体Ta和模量正演，输入正确的类型参数"), _T("出错提示"), MB_OK);
		return 0;
	}
	
	GetGrdMinMax(mag,forwardinfo.model.grddatainfo);
	//销毁二维数组
	DeleteArray2(Hax,number_y,number_x);
	DeleteArray2(Hay,number_y,number_x);
	DeleteArray2(Za,number_y,number_x);

	return 0;
}

int _3DFiniteVercicalLine_mag(double** mag, RegularGeometry3DForward& forwardinfo, int type)
{
	//首先将mag赋值为0，后面累加
	int number_x = forwardinfo.model.grddatainfo.AutoGetNumber_x();
	int number_y = forwardinfo.model.grddatainfo.AutoGetNumber_y();
	Assign_Array2(mag, number_y, number_x, 0);
	//坐标范围
	double xmin = forwardinfo.model.grddatainfo.m_AxisBounds[0];
	double xmax = forwardinfo.model.grddatainfo.m_AxisBounds[1];
	double ymin = forwardinfo.model.grddatainfo.m_AxisBounds[2];
	double ymax = forwardinfo.model.grddatainfo.m_AxisBounds[3];
	double dx = forwardinfo.model.grddatainfo.m_Dx;
	double dy = forwardinfo.model.grddatainfo.m_Dy;
	switch (type)
	{
		case FORWARD_Hax:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
				
				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double S = PI*r*r;
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y,x,y;
				double r1, r2,xx0,yy0,zz1,zz2,Vxz,Vxy,Vxx,xx02,yy02,rxy2,rxy;
				//north_x = 190;//东向切面
				//east_y = 140;//北向切面
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x - y0 ; yy0 = y - x0 ; zz1 = z1-z; zz2 = z2-z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						rxy = (xx02 + yy02);
						rxy2 = rxy*rxy;
						r1 = sqrt(xx02+yy02+zz1*zz1);
						r2 = sqrt(xx02 + yy02 + zz2*zz2);
						//--------------------------------------------------------------------
						//Vxx = (yy02 - xx02) / rxy2*(zz1 / r1 - zz2 / r2) - xx02 / rxy*(zz1/pow(r1,3.0)-zz2/pow(r2,3.0));//化简的公式
						Vxx = (-xx02 / pow(zz2 + r2, 2.0) / r2 / r2 + (yy02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0)) 
							- (-xx02 / pow(zz1 + r1, 2.0) / r1 / r1 + (yy02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
						//Vxx = Vxx*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

						//Vxy = -2 * xx0*yy0 / rxy2*(zz1 / r1-zz2 / r2) - xx0*yy0 / rxy*(zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));//化简的公式
						Vxy = (yy0*xx0 / pow(zz1 + r1, 2.0) / r1 / r1 + xx0*yy0 / (zz1 + r1) / pow(r1, 3.0)) 
							- (yy0*xx0 / pow(zz2 + r2, 2.0) / r2 / r2 + xx0*yy0 / (zz2 + r2) / pow(r2, 3.0));
						//Vxy = Vxy*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E					//有问题

						Vxz = xx0*(1.0 / pow(r2, 3.0)-1.0 / pow(r1, 3.0)); //Vxz = Vxz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

						mag[j][k] += 1.0 / 4.0 / PI*(Mx*Vxx + My*Vxy + Mz*Vxz);
					}
				}
			}
		}
		break;
		case FORWARD_Hay:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];

				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double S = PI*r*r;
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double r1, r2, xx0, yy0, zz1, zz2,Vxy, Vyy, Vyz,xx02, yy02, rxy2, rxy;
				//north_x = 190;//东向切面
				//east_y = 140;//北向切面
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x - y0; yy0 = y - x0; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						rxy = (xx02 + yy02);
						rxy2 = rxy*rxy;
						r1 = sqrt(xx02 + yy02 + zz1*zz1);
						r2 = sqrt(xx02 + yy02 + zz2*zz2);
						//--------------------------------------------------------------------

						//Vyy = (xx02-yy02) / rxy2*(zz1 / r1 - zz2 / r2) - yy02 / rxy*(zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));//化简的公式
						Vyy = (-yy02 / pow(zz2 + r2, 2.0) / r2 / r2 + (xx02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
							- (-yy02 / pow(zz1 + r1, 2.0) / r1 / r1 + (xx02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));

						//Vxy = -2 * xx0*yy0 / rxy2*(zz1 / r1-zz2 / r2) - xx0*yy0 / rxy*(zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));//化简的公式
						Vxy = (yy0*xx0 / pow(zz1 + r1, 2.0) / r1 / r1 + xx0*yy0 / (zz1 + r1) / pow(r1, 3.0))
							- (yy0*xx0 / pow(zz2 + r2, 2.0) / r2 / r2 + xx0*yy0 / (zz2 + r2) / pow(r2, 3.0));

						Vyz = yy0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0)); //Vyz = Vyz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

						mag[j][k] += 1.0 / 4.0 / PI*(Mx*Vxy + My*Vyy + Mz*Vyz);
					}
				}
			}
		}
			break;
		case FORWARD_Za:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];

				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double S = PI*r*r;
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double r1, r2, xx0, yy0, zz1, zz2, Vxz, Vyz, Vzz, xx02, yy02, rxy2, rxy;
				//north_x = 190;//东向切面
				//east_y = 140;//北向切面
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x - y0; yy0 = y - x0; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						rxy = (xx02 + yy02);
						rxy2 = rxy*rxy;
						r1 = sqrt(xx02 + yy02 + zz1*zz1);
						r2 = sqrt(xx02 + yy02 + zz2*zz2);
						//--------------------------------------------------------------------

						Vxz = xx0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0)); //Vxz = Vxz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

						Vyz = yy0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));// Vyz = Vyz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

						Vzz = (zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0)); //Vzz = Vzz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

						mag[j][k] += 1.0 / 4.0 / PI*(Mx*Vxz + My*Vyz + Mz*Vzz);
					}
				}
			}
		}
			break;

		default:
			MessageBox(NULL, _T("这是磁三分量正演函数，请输入正确的type"), _T("有限长直立线模型正演-错误提示"), MB_OK);
			return 0;

	}
	GetGrdMinMax(mag, forwardinfo.model.grddatainfo);
	return 1;
}

int _3DFiniteVercicalLine_magT(double** mag, RegularGeometry3DForward& forwardinfo, int type)
{
	//首先将grav赋值为0，后面累加
	int number_x = forwardinfo.model.grddatainfo.GetNumber_x();
	int number_y = forwardinfo.model.grddatainfo.GetNumber_y();
	double** Hax = CreateArray2(number_y, number_x);
	double** Hay = CreateArray2(number_y, number_x);
	double** Za = CreateArray2(number_y, number_x);
	Assign_Array2(Hax, number_y, number_x, 0);
	Assign_Array2(Hay, number_y, number_x, 0);
	Assign_Array2(Za, number_y, number_x, 0);
	Assign_Array2(mag, number_y, number_x, 0);
	double I = forwardinfo.EarthAngle_I / 180.0*PI;
	double D = forwardinfo.EarthAngle_D / 180.0*PI;

	//计算Hax
	_3DFiniteVercicalLine_mag(Hax, forwardinfo, FORWARD_Hax);
	//计算Hay
	_3DFiniteVercicalLine_mag(Hay, forwardinfo, FORWARD_Hay);
	//计算Za
	_3DFiniteVercicalLine_mag(Za, forwardinfo, FORWARD_Za);
	switch (type)
	{
	case FORWARD_Ta:
	{
		//计算Ta
		for (int i = 0; i < number_y; i++)
		{
			for (int j = 0; j < number_x; j++)
			{
				mag[i][j] = Hax[i][j] * cos(I)*cos(D) + Hay[i][j] * cos(I)*sin(D) + Za[i][j] * sin(I);
			}
		}
	}
		break;
	case FORWARD_Module:
	{
		//计算模量
		for (int i = 0; i < number_y; i++)
		{
			for (int j = 0; j < number_x; j++)
			{
				mag[i][j] = sqrt(Hax[i][j] * Hax[i][j] + Hay[i][j] * Hay[i][j] + Za[i][j] * Za[i][j]);
			}
		}
	}
		break;
	default:
		MessageBox(NULL, _T("这是Ta和模量正演，输入正确的类型参数"), _T("出错提示"), MB_OK);
		return 0;
	}

	GetGrdMinMax(mag, forwardinfo.model.grddatainfo);
	//销毁二维数组
	DeleteArray2(Hax, number_y, number_x);
	DeleteArray2(Hay, number_y, number_x);
	DeleteArray2(Za, number_y, number_x);

	return 0;
}
int _3DFiniteVercicalLine_grav(double** grav, RegularGeometry3DForward& forwardinfo, int type)
{
	//首先将mag赋值为0，后面累加
	int number_x = forwardinfo.model.grddatainfo.AutoGetNumber_x();
	int number_y = forwardinfo.model.grddatainfo.AutoGetNumber_y();
	Assign_Array2(grav, number_y, number_x, 0);
	//坐标范围
	double xmin = forwardinfo.model.grddatainfo.m_AxisBounds[0];
	double xmax = forwardinfo.model.grddatainfo.m_AxisBounds[1];
	double ymin = forwardinfo.model.grddatainfo.m_AxisBounds[2];
	double ymax = forwardinfo.model.grddatainfo.m_AxisBounds[3];
	double dx = forwardinfo.model.grddatainfo.m_Dx;
	double dy = forwardinfo.model.grddatainfo.m_Dy;
	switch (type)
	{
		case FORWARD_V:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];

				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double S = PI*r*r;
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
				double V;
				//north_x = 190;//东向切面
				//east_y = 140;//北向切面
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x - y0; yy0 = y - x0; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						rxy = (xx02 + yy02);
						rxy2 = rxy*rxy;
						r1 = sqrt(xx02 + yy02 + zz1*zz1);
						r2 = sqrt(xx02 + yy02 + zz2*zz2);
						//--------------------------------------------------------------------
						V = log((zz2)+r2) - log((zz1)+r1); V = V*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;
						grav[j][k] += V;			//
					}
				}
			}
		}
		break;
		case FORWARD_Vx:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];

				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double S = PI*r*r;
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
				double Vx;
				//north_x = 190;//东向切面
				//east_y = 140;//北向切面
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x - y0; yy0 = y - x0; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						rxy = (xx02 + yy02);
						rxy2 = rxy*rxy;
						r1 = sqrt(xx02 + yy02 + zz1*zz1);
						r2 = sqrt(xx02 + yy02 + zz2*zz2);
						//--------------------------------------------------------------------
						//Vx = xx0 / rxy*(zz1 / r1 - zz2 / r2); Vx = Vx*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;//化简的结果在中心投影点存在分母等于零的情况
						Vx = xx0 / (zz2 + r2) / r2 - xx0 / (zz1 + r1) / r1; Vx = Vx*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;

						grav[j][k] += Vx;			
					}
				}
			}
		}
			break;
		case FORWARD_Vy:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];

				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double S = PI*r*r;
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
				double Vy;
				//north_x = 190;//东向切面
				//east_y = 140;//北向切面
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x - y0; yy0 = y - x0; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						rxy = (xx02 + yy02);
						rxy2 = rxy*rxy;
						r1 = sqrt(xx02 + yy02 + zz1*zz1);
						r2 = sqrt(xx02 + yy02 + zz2*zz2);
						//--------------------------------------------------------------------

						//Vy = yy0 / rxy*(zz1 / r1 - zz2 / r2); Vy = Vy*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;//化简公式在中心投影点存在分母等于零的情况
						Vy = yy0 / (zz2 + r2) / r2 - yy0 / (zz1 + r1) / r1; Vy = Vy*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;

						grav[j][k] += Vy;			
					}
				}
			}
		}
			break;
		case FORWARD_Vz:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];

				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double S = PI*r*r;
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
				double Vz;
				//north_x = 190;//东向切面
				//east_y = 140;//北向切面
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x - y0; yy0 = y - x0; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						rxy = (xx02 + yy02);
						rxy2 = rxy*rxy;
						r1 = sqrt(xx02 + yy02 + zz1*zz1);
						r2 = sqrt(xx02 + yy02 + zz2*zz2);
						//--------------------------------------------------------------------

						Vz = 1.0 / r1 - 1.0 / r2; Vz = Vz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;				//mGal

						grav[j][k] += Vz;
					}
				}
			}
		}
			break;
		case FORWARD_Vxx:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];

				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double S = PI*r*r;
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double r1, r2, xx0, yy0, zz1, zz2, Vxx, xx02, yy02, rxy2, rxy;
				//north_x = 190;//东向切面
				//east_y = 140;//北向切面
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x - y0; yy0 = y - x0; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						rxy = (xx02 + yy02);
						rxy2 = rxy*rxy;
						r1 = sqrt(xx02 + yy02 + zz1*zz1);
						r2 = sqrt(xx02 + yy02 + zz2*zz2);
						//--------------------------------------------------------------------

						//Vxx = (yy02 - xx02) / rxy2*(zz1 / r1 - zz2 / r2) - xx02 / rxy*(zz1/pow(r1,3.0)-zz2/pow(r2,3.0));//化简的公式
						Vxx = (-xx02 / pow(zz2 + r2, 2.0) / r2 / r2 + (yy02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
							- (-xx02 / pow(zz1 + r1, 2.0) / r1 / r1 + (yy02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
						Vxx = Vxx*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

						grav[j][k] += Vxx;			
					}
				}
			}
		}
			break;
		case FORWARD_Vxy:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];

				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double S = PI*r*r;
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double r1, r2, xx0, yy0, zz1, zz2, Vxy, xx02, yy02, rxy2, rxy;
				//north_x = 190;//东向切面
				//east_y = 140;//北向切面
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x - y0; yy0 = y - x0; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						rxy = (xx02 + yy02);
						rxy2 = rxy*rxy;
						r1 = sqrt(xx02 + yy02 + zz1*zz1);
						r2 = sqrt(xx02 + yy02 + zz2*zz2);
						//--------------------------------------------------------------------
						//Vxy = -2 * xx0*yy0 / rxy2*(zz1 / r1-zz2 / r2) - xx0*yy0 / rxy*(zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));//化简的公式
						Vxy = (yy0*xx0 / pow(zz1 + r1, 2.0) / r1 / r1 + xx0*yy0 / (zz1 + r1) / pow(r1, 3.0))
							- (yy0*xx0 / pow(zz2 + r2, 2.0) / r2 / r2 + xx0*yy0 / (zz2 + r2) / pow(r2, 3.0));
						Vxy = Vxy*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

						grav[j][k] += Vxy;
					}
				}
			}
		}
			break;
		case FORWARD_Vxz:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];

				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double S = PI*r*r;
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double r1, r2, xx0, yy0, zz1, zz2, Vxz, xx02, yy02, rxy2, rxy;
				//north_x = 190;//东向切面
				//east_y = 140;//北向切面
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x - y0; yy0 = y - x0; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						rxy = (xx02 + yy02);
						rxy2 = rxy*rxy;
						r1 = sqrt(xx02 + yy02 + zz1*zz1);
						r2 = sqrt(xx02 + yy02 + zz2*zz2);
						//--------------------------------------------------------------------

						Vxz = xx0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0)); Vxz = Vxz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

						grav[j][k] += Vxz;			//Vxy有点问题
					}
				}
			}
		}
			break;
		case FORWARD_Vyy:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];

				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double S = PI*r*r;
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double r1, r2, xx0, yy0, zz1, zz2, Vyy, xx02, yy02, rxy2, rxy;
				//north_x = 190;//东向切面
				//east_y = 140;//北向切面
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x - y0; yy0 = y - x0; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						rxy = (xx02 + yy02);
						rxy2 = rxy*rxy;
						r1 = sqrt(xx02 + yy02 + zz1*zz1);
						r2 = sqrt(xx02 + yy02 + zz2*zz2);
						//--------------------------------------------------------------------

						////Vyy = (xx02-yy02) / rxy2*(zz1 / r1 - zz2 / r2) - yy02 / rxy*(zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));//化简的公式
						Vyy = (-yy02 / pow(zz2 + r2, 2.0) / r2 / r2 + (xx02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
							- (-yy02 / pow(zz1 + r1, 2.0) / r1 / r1 + (xx02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
						Vyy = Vyy*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

						grav[j][k] += Vyy;			//Vxy有点问题
					}
				}
			}
		}
			break;
		case FORWARD_Vyz:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];

				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double S = PI*r*r;
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double r1, r2, xx0, yy0, zz1, zz2, Vyz,xx02, yy02, rxy2, rxy;
				//north_x = 190;//东向切面
				//east_y = 140;//北向切面
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x - y0; yy0 = y - x0; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						rxy = (xx02 + yy02);
						rxy2 = rxy*rxy;
						r1 = sqrt(xx02 + yy02 + zz1*zz1);
						r2 = sqrt(xx02 + yy02 + zz2*zz2);
						//--------------------------------------------------------------------
						Vyz = yy0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0)); Vyz = Vyz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

						grav[j][k] += Vyz;			//Vxy有点问题
					}
				}
			}
		}
			break;
		case FORWARD_Vzz:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				//计算圆柱体水平偏离y轴的夹角:用这两个端点确定一个角度，使用atan2函数表示这个方向与x轴正向的夹角，我们只需要0~180度即可，因此这里对atan2的结果取绝对值
				//然后用90度减去这个角度，如果是正值说明观测坐标系顺时针转可以使得y轴与圆柱体平行；如果为负值则需要沿逆时针方向旋转，两种情况的坐标变换矩阵不同而已
				double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];

				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				//==================================================上面是平面旋转的代码=====================================================================
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double S = PI*r*r;
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double r1, r2, xx0, yy0, zz1, zz2,Vzz, xx02, yy02, rxy2, rxy;
				//north_x = 190;//东向切面
				//east_y = 140;//北向切面
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x - y0; yy0 = y - x0; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						rxy = (xx02 + yy02);
						rxy2 = rxy*rxy;
						r1 = sqrt(xx02 + yy02 + zz1*zz1);
						r2 = sqrt(xx02 + yy02 + zz2*zz2);
						//--------------------------------------------------------------------
					
						Vzz = (zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0)); Vzz = Vzz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

						grav[j][k] += Vzz;			//Vxy有点问题
					}
				}
			}
		}
			break;
	default:
		MessageBox(NULL, _T("这是重力位及其导数正演函数，请输入正确的type"), _T("有限长直立线模型正演-错误提示"), MB_OK);
		return 0;

	}
	GetGrdMinMax(grav, forwardinfo.model.grddatainfo);
	return 1;
}

int _3DFiniteLine_Grav(double** grav, RegularGeometry3DForward& forwardinfo0, int type)
{
	RegularGeometry3DForward forwardinfo = forwardinfo0;
	//首先将mag赋值为0，后面累加
	int number_x = forwardinfo.model.grddatainfo.AutoGetNumber_x();
	int number_y = forwardinfo.model.grddatainfo.AutoGetNumber_y();
	Assign_Array2(grav, number_y, number_x, 0);
	//坐标范围
	double xmin = forwardinfo.model.grddatainfo.m_AxisBounds[0];
	double xmax = forwardinfo.model.grddatainfo.m_AxisBounds[1];
	double ymin = forwardinfo.model.grddatainfo.m_AxisBounds[2];
	double ymax = forwardinfo.model.grddatainfo.m_AxisBounds[3];
	double dx = forwardinfo.model.grddatainfo.m_Dx;
	double dy = forwardinfo.model.grddatainfo.m_Dy;
	switch (type)
	{
		case FORWARD_V:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
				{
					for (int ii = 0; ii < 3; ii++)
					{
						double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
						forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
						forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
					}
				}
				//0. 与坐标系无关的量
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double S = PI*r*r;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				//1. 根据线模型的两个端点计算转换矩阵
				double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
				forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
				temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
				forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
				double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
				
				if (x1==x2 && y1==y2)		//直立模型
				{
					double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
					double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double V;
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						//z = ymin+dy*j;//计算沿东向的切片
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//z = xmin +dx*k;//沿北向的切片
							//将计算点变换到圆柱体坐标系中
							x = north_x; y = east_y;
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							V = log((zz2)+r2) - log((zz1)+r1); V = V*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;
							grav[j][k] += V;			//Vxy有点问题
						}
					}
				}
				else			//非直立模型
				{
					//1. 圆柱体的计算倾角和偏角
					double direct_vector_horizontal[3],direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
					double L, L_H;
					L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
					L = sqrt(L_H*L_H+(z1-z2)*(z1-z2));	//圆柱体长度
					direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
					direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
					double *cross_horizontal, dot_horizontal, D,I;//x与水平投影方向的叉乘及点乘，偏角,倾角
					cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
					//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
					dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
					//printf("%lf\n",dot_horizontal);
					D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
					if (dot_horizontal < 0)
					{
						if (D < 0)	//第三象限
						{
							D = -(D + PI);
						}
						else		//第四象限
						{
							D = PI - D;
						}
					}
					//printf("%lf\n", D*180/PI);
					I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角
					
					//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
					double TransMat[9];
					double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
					TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
					TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
					TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;
					
					//3. 根据变换矩阵计算变换后的点
					double newPos1[3], newPos2[3];
					Mat_Multiply(TransMat,forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
					
					//4. 将观测坐标转换到模型坐标系进行计算
					double x0 = newPos2[0];
					double y0 = newPos2[1];
					z1 = newPos1[2], z2 = newPos2[2];
					//printf("%lf %lf \n%lf \n%lf \n",x0,y0,z1,z2);
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
					angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
					if (angle_I < 0)				//当磁化方向为向上时
					{
						angle_I = PI / 2.0 - (fabs(angle_I)+PI/2.0 - I);		
					}
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D-D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
					//********************
					//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
					//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double V;
					double SurveyPoint[3], TransPoint[3];
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//将计算点变换到圆柱体坐标系中
							SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
							Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
							x = TransPoint[0]; y = TransPoint[1];
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							V = log((zz2)+r2) - log((zz1)+r1); V = V*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;
							grav[j][k] += V;			//Vxy有点问题
						}
					}
				}
			}
		}
		break;
		case FORWARD_Vx:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
				{
					for (int ii = 0; ii < 3; ii++)
					{
						double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
						forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
						forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
					}
				}
				//0. 与坐标系无关的量
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double S = PI*r*r;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				//1. 根据线模型的两个端点计算转换矩阵
				double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
				forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
				temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
				forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
				double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
				
				if (x1 == x2 && y1 == y2)		//直立模型
				{
					double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
					double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vx;
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						//z = ymin+dy*j;//计算沿东向的切片
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//z = xmin +dx*k;//沿北向的切片
							//将计算点变换到圆柱体坐标系中
							x = north_x; y = east_y;
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							Vx = xx0 / (zz2 + r2) / r2 - xx0 / (zz1 + r1) / r1; Vx = Vx*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;
							grav[j][k] += Vx;			//Vxy有点问题
						}
					}
				}
				else			//非直立模型
				{
					//1. 圆柱体的计算倾角和偏角
					double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
					double L, L_H;
					L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
					L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
					direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
					direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
					double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
					cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
					//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
					dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
					//printf("%lf\n",dot_horizontal);
					D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
					if (dot_horizontal < 0)
					{
						if (D < 0)	//第三象限
						{
							D = -(D + PI);
						}
						else		//第四象限
						{
							D = PI - D;
						}
					}
					//printf("%lf\n", D*180/PI);
					I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

					//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
					double TransMat[9];
					double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
					TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
					TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
					TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

					//3. 根据变换矩阵计算变换后的点
					double newPos1[3], newPos2[3];
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
					//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
					double point_x[3]; point_x[0] = 1; point_x[1] = 0; point_x[2] = 0;
					double DirectCos[3];
					Mat_Multiply(TransMat, point_x, DirectCos, 3, 3);//北向在模型坐标系下的方向余弦
					//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
					//4. 将观测坐标转换到模型坐标系进行计算
					double x0 = newPos2[0];
					double y0 = newPos2[1];
					z1 = newPos1[2], z2 = newPos2[2];
					//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
					angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
					if (angle_I < 0)				//当磁化方向为向上时
					{
						angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
					}
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
					//********************
					//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
					//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vx,Vy,Vz;
					double SurveyPoint[3], TransPoint[3];
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//将计算点变换到圆柱体坐标系中
							SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
							Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
							x = TransPoint[0]; y = TransPoint[1];
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							//Vx = xx0 / rxy*(zz1 / r1 - zz2 / r2); Vx = Vx*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;//化简的结果在中心投影点存在分母等于零的情况
							Vx = xx0 / (zz2 + r2) / r2 - xx0 / (zz1 + r1) / r1; 
							Vy = yy0 / (zz2 + r2) / r2 - yy0 / (zz1 + r1) / r1; 
							Vz = 1.0 / r1 - 1.0 / r2; 
							grav[j][k] += (Vx*DirectCos[0] + Vy*DirectCos[1] + Vz*DirectCos[2])*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;//mGal/m
						}
					}
				}
			}
		}
		break;
		case FORWARD_Vy:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
				{
					for (int ii = 0; ii < 3; ii++)
					{
						double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
						forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
						forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
					}
				}
				//0. 与坐标系无关的量
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double S = PI*r*r;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				//1. 根据线模型的两个端点计算转换矩阵
				double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
				forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
				temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
				forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
				double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
				
				if (x1 == x2 && y1 == y2)		//直立模型
				{
					double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
					double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vy;
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						//z = ymin+dy*j;//计算沿东向的切片
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//z = xmin +dx*k;//沿北向的切片
							//将计算点变换到圆柱体坐标系中
							x = north_x; y = east_y;
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							Vy = yy0 / (zz2 + r2) / r2 - yy0 / (zz1 + r1) / r1; Vy = Vy*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;
							grav[j][k] += Vy;			//Vxy有点问题
						}
					}
				}
				else			//非直立模型
				{
					//1. 圆柱体的计算倾角和偏角
					double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
					double L, L_H;
					L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
					L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
					direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
					direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
					double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
					cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
					//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
					dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
					//printf("%lf\n",dot_horizontal);
					D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
					if (dot_horizontal < 0)
					{
						if (D < 0)	//第三象限
						{
							D = -(D + PI);
						}
						else		//第四象限
						{
							D = PI - D;
						}
					}
					//printf("%lf\n", D*180/PI);
					I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

					//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
					double TransMat[9];
					double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
					TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
					TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
					TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

					//3. 根据变换矩阵计算变换后的点
					double newPos1[3], newPos2[3];
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
					//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
					double point_y[3]; point_y[0] = 0; point_y[1] = 1; point_y[2] = 0;
					double DirectCos[3];
					Mat_Multiply(TransMat, point_y, DirectCos, 3, 3);//北向在模型坐标系下的方向余弦
					//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
					//4. 将观测坐标转换到模型坐标系进行计算
					double x0 = newPos2[0];
					double y0 = newPos2[1];
					z1 = newPos1[2], z2 = newPos2[2];
					//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
					angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
					if (angle_I < 0)				//当磁化方向为向上时
					{
						angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
					}
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
					//********************
					//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
					//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vx, Vy, Vz;
					double SurveyPoint[3], TransPoint[3];
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//将计算点变换到圆柱体坐标系中
							SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
							Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
							x = TransPoint[0]; y = TransPoint[1];
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							//Vx = xx0 / rxy*(zz1 / r1 - zz2 / r2); Vx = Vx*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;//化简的结果在中心投影点存在分母等于零的情况
							Vx = xx0 / (zz2 + r2) / r2 - xx0 / (zz1 + r1) / r1;
							Vy = yy0 / (zz2 + r2) / r2 - yy0 / (zz1 + r1) / r1;
							Vz = 1.0 / r1 - 1.0 / r2;
							grav[j][k] += (Vx*DirectCos[0] + Vy*DirectCos[1] + Vz*DirectCos[2])*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;//mGal/m
						}
					}
				}
			}
		}
		break;
		case FORWARD_Vz:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
				{
					for (int ii = 0; ii < 3; ii++)
					{
						double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
						forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
						forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
					}
				}
				//0. 与坐标系无关的量
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double S = PI*r*r;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				//1. 根据线模型的两个端点计算转换矩阵
				double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
				forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
				temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
				forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
				double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
				
				if (x1 == x2 && y1 == y2)		//直立模型
				{
					double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
					double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vz;
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						//z = ymin+dy*j;//计算沿东向的切片
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//z = xmin +dx*k;//沿北向的切片
							//将计算点变换到圆柱体坐标系中
							x = north_x; y = east_y;
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							Vz = 1.0 / r1 - 1.0 / r2; Vz = Vz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;
							grav[j][k] += Vz;			//Vxy有点问题
						}
					}
				}
				else			//非直立模型
				{
					//1. 圆柱体的计算倾角和偏角
					double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
					double L, L_H;
					L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
					L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
					direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
					direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
					double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
					cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
					//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
					dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
					//printf("%lf\n",dot_horizontal);
					D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
					if (dot_horizontal < 0)
					{
						if (D < 0)	//第三象限
						{
							D = -(D + PI);
						}
						else		//第四象限
						{
							D = PI - D;
						}
					}
					//printf("%lf\n", D*180/PI);
					I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

					//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
					double TransMat[9];
					double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
					TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
					TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
					TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

					//3. 根据变换矩阵计算变换后的点
					double newPos1[3], newPos2[3];
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
					//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
					double point_z[3]; point_z[0] = 0; point_z[1] = 0; point_z[2] = 1;
					double DirectCos[3];
					Mat_Multiply(TransMat, point_z, DirectCos, 3, 3);//北向在模型坐标系下的方向余弦
					//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
					//4. 将观测坐标转换到模型坐标系进行计算
					double x0 = newPos2[0];
					double y0 = newPos2[1];
					z1 = newPos1[2], z2 = newPos2[2];
					//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
					angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
					if (angle_I < 0)				//当磁化方向为向上时
					{
						angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
					}
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
					//********************
					//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
					//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vx, Vy, Vz;
					double SurveyPoint[3], TransPoint[3];
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//将计算点变换到圆柱体坐标系中
							SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
							Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
							x = TransPoint[0]; y = TransPoint[1];
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							//Vx = xx0 / rxy*(zz1 / r1 - zz2 / r2); Vx = Vx*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;//化简的结果在中心投影点存在分母等于零的情况
							Vx = xx0 / (zz2 + r2) / r2 - xx0 / (zz1 + r1) / r1;
							Vy = yy0 / (zz2 + r2) / r2 - yy0 / (zz1 + r1) / r1;
							Vz = 1.0 / r1 - 1.0 / r2;
							grav[j][k] += (Vx*DirectCos[0] + Vy*DirectCos[1] + Vz*DirectCos[2])*G*forwardinfo.model.cylinder_vec[i].Density*S*1E8;//mGal/m
						}
					}
				}
			}
		}
		break;
		case FORWARD_Vxx:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
				{
					for (int ii = 0; ii < 3; ii++)
					{
						double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
						forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
						forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
					}
				}
				//0. 与坐标系无关的量
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double S = PI*r*r;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				//1. 根据线模型的两个端点计算转换矩阵
				double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
				forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
				temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
				forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
				double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
				
				if (x1 == x2 && y1 == y2)		//直立模型
				{
					double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
					double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vxx;
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						//z = ymin+dy*j;//计算沿东向的切片
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//z = xmin +dx*k;//沿北向的切片
							//将计算点变换到圆柱体坐标系中
							x = north_x; y = east_y;
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							Vxx = (-xx02 / pow(zz2 + r2, 2.0) / r2 / r2 + (yy02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
								- (-xx02 / pow(zz1 + r1, 2.0) / r1 / r1 + (yy02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
							Vxx = Vxx*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
							grav[j][k] += Vxx;
						}
					}
				}
				else			//非直立模型
				{
					//1. 圆柱体的计算倾角和偏角
					double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
					double L, L_H;
					L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
					L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
					direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
					direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
					double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
					cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
					//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
					dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
					//printf("%lf\n",dot_horizontal);
					D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
					if (dot_horizontal < 0)
					{
						if (D < 0)	//第三象限
						{
							D = -(D + PI);
						}
						else		//第四象限
						{
							D = PI - D;
						}
					}
					//printf("%lf\n", D*180/PI);
					I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

					//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
					double TransMat[9];
					double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
					TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
					TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
					TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

					//3. 根据变换矩阵计算变换后的点
					double newPos1[3], newPos2[3];
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
					//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
					double point_x[3]; point_x[0] = 1; point_x[1] = 0; point_x[2] = 0;
					double DirectCos[3];
					Mat_Multiply(TransMat, point_x, DirectCos, 3, 3);//北向在模型坐标系下的方向余弦
					//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
					//4. 将观测坐标转换到模型坐标系进行计算
					double x0 = newPos2[0];
					double y0 = newPos2[1];
					z1 = newPos1[2], z2 = newPos2[2];
					//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
					angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
					if (angle_I < 0)				//当磁化方向为向上时
					{
						angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
					}
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
					//********************
					//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
					//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vxx,Vxy,Vxz,Vyy,Vyz,Vzz;
					double SurveyPoint[3], TransPoint[3];
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//将计算点变换到圆柱体坐标系中
							SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
							Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
							x = TransPoint[0]; y = TransPoint[1];
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							Vxx = (-xx02 / pow(zz2 + r2, 2.0) / r2 / r2 + (yy02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
								- (-xx02 / pow(zz1 + r1, 2.0) / r1 / r1 + (yy02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
							Vxy = (yy0*xx0 / pow(zz1 + r1, 2.0) / r1 / r1 + xx0*yy0 / (zz1 + r1) / pow(r1, 3.0))
								- (yy0*xx0 / pow(zz2 + r2, 2.0) / r2 / r2 + xx0*yy0 / (zz2 + r2) / pow(r2, 3.0));//*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
							Vxz = xx0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
							Vyy = (-yy02 / pow(zz2 + r2, 2.0) / r2 / r2 + (xx02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
								- (-yy02 / pow(zz1 + r1, 2.0) / r1 / r1 + (xx02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
							Vyz = yy0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
							Vzz = (zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));
							grav[j][k] += (DirectCos[0] * (Vxx*DirectCos[0] + Vxy*DirectCos[1] + Vxz*DirectCos[2]) + 
								DirectCos[1] *	(Vxy*DirectCos[0] + Vyy*DirectCos[1] + Vyz*DirectCos[2])+
								DirectCos[2] * (Vxz*DirectCos[0] + Vyz*DirectCos[1] + Vzz*DirectCos[2]))*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
						}
					}
				}
			}
		}
		break;
		case FORWARD_Vxy:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
				{
					for (int ii = 0; ii < 3; ii++)
					{
						double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
						forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
						forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
					}
				}
				//0. 与坐标系无关的量
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double S = PI*r*r;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				//1. 根据线模型的两个端点计算转换矩阵
				double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
				forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
				temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
				forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
				double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
				
				if (x1 == x2 && y1 == y2)		//直立模型
				{
					double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
					double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vxy;
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						//z = ymin+dy*j;//计算沿东向的切片
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//z = xmin +dx*k;//沿北向的切片
							//将计算点变换到圆柱体坐标系中
							x = north_x; y = east_y;
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							Vxy = (yy0*xx0 / pow(zz1 + r1, 2.0) / r1 / r1 + xx0*yy0 / (zz1 + r1) / pow(r1, 3.0))
								- (yy0*xx0 / pow(zz2 + r2, 2.0) / r2 / r2 + xx0*yy0 / (zz2 + r2) / pow(r2, 3.0));
							Vxy = Vxy*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
							grav[j][k] += Vxy;
						}
					}
				}
				else			//非直立模型
				{
					//1. 圆柱体的计算倾角和偏角
					double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
					double L, L_H;
					L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
					L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
					direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
					direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
					double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
					cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
					//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
					dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
					//printf("%lf\n",dot_horizontal);
					D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
					if (dot_horizontal < 0)
					{
						if (D < 0)	//第三象限
						{
							D = -(D + PI);
						}
						else		//第四象限
						{
							D = PI - D;
						}
					}
					//printf("%lf\n", D*180/PI);
					I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

					//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
					double TransMat[9];
					double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
					TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
					TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
					TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

					//3. 根据变换矩阵计算变换后的点
					double newPos1[3], newPos2[3];
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
					//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
					double point_x[3]; point_x[0] = 1; point_x[1] = 0; point_x[2] = 0;
					double point_y[3]; point_y[0] = 0; point_y[1] = 1; point_y[2] = 0;
					double XDirectCos[3],YDirectCos[3];
					Mat_Multiply(TransMat, point_x, XDirectCos, 3, 3);//北向在模型坐标系下的方向余弦
					Mat_Multiply(TransMat, point_y, YDirectCos, 3, 3);
					//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
					//4. 将观测坐标转换到模型坐标系进行计算
					double x0 = newPos2[0];
					double y0 = newPos2[1];
					z1 = newPos1[2], z2 = newPos2[2];
					//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
					angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
					if (angle_I < 0)				//当磁化方向为向上时
					{
						angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
					}
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
					//********************
					//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
					//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
					double SurveyPoint[3], TransPoint[3];
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//将计算点变换到圆柱体坐标系中
							SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
							Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
							x = TransPoint[0]; y = TransPoint[1];
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							Vxx = (-xx02 / pow(zz2 + r2, 2.0) / r2 / r2 + (yy02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
								- (-xx02 / pow(zz1 + r1, 2.0) / r1 / r1 + (yy02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
							Vxy = (yy0*xx0 / pow(zz1 + r1, 2.0) / r1 / r1 + xx0*yy0 / (zz1 + r1) / pow(r1, 3.0))
								- (yy0*xx0 / pow(zz2 + r2, 2.0) / r2 / r2 + xx0*yy0 / (zz2 + r2) / pow(r2, 3.0));//*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
							Vxz = xx0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
							Vyy = (-yy02 / pow(zz2 + r2, 2.0) / r2 / r2 + (xx02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
								- (-yy02 / pow(zz1 + r1, 2.0) / r1 / r1 + (xx02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
							Vyz = yy0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
							Vzz = (zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));
							grav[j][k] += (YDirectCos[0] * (Vxx*XDirectCos[0] + Vxy*XDirectCos[1] + Vxz*XDirectCos[2]) +
								YDirectCos[1] * (Vxy*XDirectCos[0] + Vyy*XDirectCos[1] + Vyz*XDirectCos[2]) +
								YDirectCos[2] * (Vxz*XDirectCos[0] + Vyz*XDirectCos[1] + Vzz*XDirectCos[2]))*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
						}
					}
				}
			}
		}
			break;
			case FORWARD_Vxz:
			{
				for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
				{
					if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
					{
						for (int ii = 0; ii < 3; ii++)
						{
							double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
							forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
							forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
						}
					}
					//0. 与坐标系无关的量
					double r = forwardinfo.model.cylinder_vec[i].Radius;
					double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
					double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
					double S = PI*r*r;
					double EarthMag = forwardinfo.EarthMag;
					double M = S*CHL*EarthMag;											//与系数消去了U0
					//1. 根据线模型的两个端点计算转换矩阵
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
					forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
					temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
					forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
					double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
					double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
					
					if (x1 == x2 && y1 == y2)		//直立模型
					{
						double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
						double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
						double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
						double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
						double Mx = M*cos(angle_I)*sin(angle_D);
						double My = M*cos(angle_I)*cos(angle_D);
						double Mz = M*sin(angle_I);
						double north_x, east_y, x, y;
						double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
						double Vxz;
						for (int j = 0; j < number_y; j++)
						{
							north_x = ymin + dy*j;
							//z = ymin+dy*j;//计算沿东向的切片
							for (int k = 0; k < number_x; k++)
							{
								east_y = xmin + k*dx;
								//z = xmin +dx*k;//沿北向的切片
								//将计算点变换到圆柱体坐标系中
								x = north_x; y = east_y;
								xx0 = x0-x; yy0 = y0-y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
								xx02 = xx0*xx0; yy02 = yy0*yy0;
								rxy = (xx02 + yy02);
								rxy2 = rxy*rxy;
								r1 = sqrt(xx02 + yy02 + zz1*zz1);
								r2 = sqrt(xx02 + yy02 + zz2*zz2);
								//--------------------------------------------------------------------
								Vxz = xx0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
								Vxz = Vxz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
								grav[j][k] += Vxz;
							}
						}
					}
					else			//非直立模型
					{
						//1. 圆柱体的计算倾角和偏角
						double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
						double L, L_H;
						L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
						L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
						direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
						direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
						double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
						cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
						//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
						dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
						//printf("%lf\n",dot_horizontal);
						D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
						if (dot_horizontal < 0)
						{
							if (D < 0)	//第三象限
							{
								D = -(D + PI);
							}
							else		//第四象限
							{
								D = PI - D;
							}
						}
						//printf("%lf\n", D*180/PI);
						I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

						//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
						double TransMat[9];
						double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
						TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
						TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
						TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

						//3. 根据变换矩阵计算变换后的点
						double newPos1[3], newPos2[3];
						Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
						Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
						//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
						double point_x[3]; point_x[0] = 1; point_x[1] = 0; point_x[2] = 0;
						double point_z[3]; point_z[0] = 0; point_z[1] = 0; point_z[2] = 1;
						double XDirectCos[3], ZDirectCos[3];
						Mat_Multiply(TransMat, point_x, XDirectCos, 3, 3);//北向在模型坐标系下的方向余弦
						Mat_Multiply(TransMat, point_z, ZDirectCos, 3, 3);
						//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
						//4. 将观测坐标转换到模型坐标系进行计算
						double x0 = newPos2[0];
						double y0 = newPos2[1];
						z1 = newPos1[2], z2 = newPos2[2];
						//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
						double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
						angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
						if (angle_I < 0)				//当磁化方向为向上时
						{
							angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
						}
						double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
						//********************
						//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
						//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
						double Mx = M*cos(angle_I)*sin(angle_D);
						double My = M*cos(angle_I)*cos(angle_D);
						double Mz = M*sin(angle_I);
						double north_x, east_y, x, y;
						double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
						double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
						double SurveyPoint[3], TransPoint[3];
						for (int j = 0; j < number_y; j++)
						{
							north_x = ymin + dy*j;
							for (int k = 0; k < number_x; k++)
							{
								east_y = xmin + k*dx;
								//将计算点变换到圆柱体坐标系中
								SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
								Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
								x = TransPoint[0]; y = TransPoint[1];
								xx0 = x0-x; yy0 = y0-y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
								xx02 = xx0*xx0; yy02 = yy0*yy0;
								rxy = (xx02 + yy02);
								rxy2 = rxy*rxy;
								r1 = sqrt(xx02 + yy02 + zz1*zz1);
								r2 = sqrt(xx02 + yy02 + zz2*zz2);
								//--------------------------------------------------------------------
								Vxx = (-xx02 / pow(zz2 + r2, 2.0) / r2 / r2 + (yy02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
									- (-xx02 / pow(zz1 + r1, 2.0) / r1 / r1 + (yy02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
								Vxy = (yy0*xx0 / pow(zz1 + r1, 2.0) / r1 / r1 + xx0*yy0 / (zz1 + r1) / pow(r1, 3.0))
									- (yy0*xx0 / pow(zz2 + r2, 2.0) / r2 / r2 + xx0*yy0 / (zz2 + r2) / pow(r2, 3.0));//*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
								Vxz = xx0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
								Vyy = (-yy02 / pow(zz2 + r2, 2.0) / r2 / r2 + (xx02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
									- (-yy02 / pow(zz1 + r1, 2.0) / r1 / r1 + (xx02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
								Vyz = yy0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
								Vzz = (zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));
								grav[j][k] += (ZDirectCos[0] * (Vxx*XDirectCos[0] + Vxy*XDirectCos[1] + Vxz*XDirectCos[2]) +
									ZDirectCos[1] * (Vxy*XDirectCos[0] + Vyy*XDirectCos[1] + Vyz*XDirectCos[2]) +
									ZDirectCos[2] * (Vxz*XDirectCos[0] + Vyz*XDirectCos[1] + Vzz*XDirectCos[2]))*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
							}
						}
					}
				}
			}
			break;
			case FORWARD_Vyy:
			{
				for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
				{
					if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
					{
						for (int ii = 0; ii < 3; ii++)
						{
							double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
							forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
							forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
						}
					}
					//0. 与坐标系无关的量
					double r = forwardinfo.model.cylinder_vec[i].Radius;
					double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
					double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
					double S = PI*r*r;
					double EarthMag = forwardinfo.EarthMag;
					double M = S*CHL*EarthMag;											//与系数消去了U0
					//1. 根据线模型的两个端点计算转换矩阵
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
					forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
					temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
					forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
					double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
					double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
					
					if (x1 == x2 && y1 == y2)		//直立模型
					{
						double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
						double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
						double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
						double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
						double Mx = M*cos(angle_I)*sin(angle_D);
						double My = M*cos(angle_I)*cos(angle_D);
						double Mz = M*sin(angle_I);
						double north_x, east_y, x, y;
						double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
						double Vyy;
						for (int j = 0; j < number_y; j++)
						{
							north_x = ymin + dy*j;
							//z = ymin+dy*j;//计算沿东向的切片
							for (int k = 0; k < number_x; k++)
							{
								east_y = xmin + k*dx;
								//z = xmin +dx*k;//沿北向的切片
								//将计算点变换到圆柱体坐标系中
								x = north_x; y = east_y;
								xx0 = x0-x; yy0 = y0-y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
								xx02 = xx0*xx0; yy02 = yy0*yy0;
								rxy = (xx02 + yy02);
								rxy2 = rxy*rxy;
								r1 = sqrt(xx02 + yy02 + zz1*zz1);
								r2 = sqrt(xx02 + yy02 + zz2*zz2);
								//--------------------------------------------------------------------
								Vyy = (-yy02 / pow(zz2 + r2, 2.0) / r2 / r2 + (xx02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
									- (-yy02 / pow(zz1 + r1, 2.0) / r1 / r1 + (xx02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
								Vyy = Vyy*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
								grav[j][k] += Vyy;
							}
						}
					}
					else			//非直立模型
					{
						//1. 圆柱体的计算倾角和偏角
						double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
						double L, L_H;
						L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
						L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
						direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
						direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
						double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
						cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
						//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
						dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
						//printf("%lf\n",dot_horizontal);
						D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
						if (dot_horizontal < 0)
						{
							if (D < 0)	//第三象限
							{
								D = -(D + PI);
							}
							else		//第四象限
							{
								D = PI - D;
							}
						}
						//printf("%lf\n", D*180/PI);
						I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

						//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
						double TransMat[9];
						double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
						TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
						TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
						TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

						//3. 根据变换矩阵计算变换后的点
						double newPos1[3], newPos2[3];
						Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
						Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
						//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
						double point_y[3]; point_y[0] = 0; point_y[1] = 1; point_y[2] = 0;
						double DirectCos[3];
						Mat_Multiply(TransMat, point_y, DirectCos, 3, 3);//北向在模型坐标系下的方向余弦
						//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
						//4. 将观测坐标转换到模型坐标系进行计算
						double x0 = newPos2[0];
						double y0 = newPos2[1];
						z1 = newPos1[2], z2 = newPos2[2];
						//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
						double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
						angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
						if (angle_I < 0)				//当磁化方向为向上时
						{
							angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
						}
						double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
						//********************
						//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
						//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
						double Mx = M*cos(angle_I)*sin(angle_D);
						double My = M*cos(angle_I)*cos(angle_D);
						double Mz = M*sin(angle_I);
						double north_x, east_y, x, y;
						double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
						double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
						double SurveyPoint[3], TransPoint[3];
						for (int j = 0; j < number_y; j++)
						{
							north_x = ymin + dy*j;
							for (int k = 0; k < number_x; k++)
							{
								east_y = xmin + k*dx;
								//将计算点变换到圆柱体坐标系中
								SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
								Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
								x = TransPoint[0]; y = TransPoint[1];
								xx0 = x0-x; yy0 = y0-y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
								xx02 = xx0*xx0; yy02 = yy0*yy0;
								rxy = (xx02 + yy02);
								rxy2 = rxy*rxy;
								r1 = sqrt(xx02 + yy02 + zz1*zz1);
								r2 = sqrt(xx02 + yy02 + zz2*zz2);
								//--------------------------------------------------------------------
								Vxx = (-xx02 / pow(zz2 + r2, 2.0) / r2 / r2 + (yy02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
									- (-xx02 / pow(zz1 + r1, 2.0) / r1 / r1 + (yy02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
								Vxy = (yy0*xx0 / pow(zz1 + r1, 2.0) / r1 / r1 + xx0*yy0 / (zz1 + r1) / pow(r1, 3.0))
									- (yy0*xx0 / pow(zz2 + r2, 2.0) / r2 / r2 + xx0*yy0 / (zz2 + r2) / pow(r2, 3.0));//*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
								Vxz = xx0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
								Vyy = (-yy02 / pow(zz2 + r2, 2.0) / r2 / r2 + (xx02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
									- (-yy02 / pow(zz1 + r1, 2.0) / r1 / r1 + (xx02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
								Vyz = yy0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
								Vzz = (zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));
								grav[j][k] += (DirectCos[0] * (Vxx*DirectCos[0] + Vxy*DirectCos[1] + Vxz*DirectCos[2]) +
									DirectCos[1] * (Vxy*DirectCos[0] + Vyy*DirectCos[1] + Vyz*DirectCos[2]) +
									DirectCos[2] * (Vxz*DirectCos[0] + Vyz*DirectCos[1] + Vzz*DirectCos[2]))*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
							}
						}
					}
				}
			}
			break;
			case FORWARD_Vyz:
			{
				for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
				{
					if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
					{
						for (int ii = 0; ii < 3; ii++)
						{
							double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
							forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
							forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
						}
					}
					//0. 与坐标系无关的量
					double r = forwardinfo.model.cylinder_vec[i].Radius;
					double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
					double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
					double S = PI*r*r;
					double EarthMag = forwardinfo.EarthMag;
					double M = S*CHL*EarthMag;											//与系数消去了U0
					//1. 根据线模型的两个端点计算转换矩阵
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
					forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
					temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
					forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
					double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
					double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
					
					if (x1 == x2 && y1 == y2)		//直立模型
					{
						double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
						double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
						double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
						double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
						double Mx = M*cos(angle_I)*sin(angle_D);
						double My = M*cos(angle_I)*cos(angle_D);
						double Mz = M*sin(angle_I);
						double north_x, east_y, x, y;
						double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
						double Vyz;
						for (int j = 0; j < number_y; j++)
						{
							north_x = ymin + dy*j;
							//z = ymin+dy*j;//计算沿东向的切片
							for (int k = 0; k < number_x; k++)
							{
								east_y = xmin + k*dx;
								//z = xmin +dx*k;//沿北向的切片
								//将计算点变换到圆柱体坐标系中
								x = north_x; y = east_y;
								xx0 = x0-x; yy0 = y0-y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
								xx02 = xx0*xx0; yy02 = yy0*yy0;
								rxy = (xx02 + yy02);
								rxy2 = rxy*rxy;
								r1 = sqrt(xx02 + yy02 + zz1*zz1);
								r2 = sqrt(xx02 + yy02 + zz2*zz2);
								//--------------------------------------------------------------------
								Vyz = yy0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
								Vyz = Vyz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
								grav[j][k] += Vyz;
							}
						}
					}
					else			//非直立模型
					{
						//1. 圆柱体的计算倾角和偏角
						double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
						double L, L_H;
						L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
						L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
						direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
						direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
						double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
						cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
						//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
						dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
						//printf("%lf\n",dot_horizontal);
						D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
						if (dot_horizontal < 0)
						{
							if (D < 0)	//第三象限
							{
								D = -(D + PI);
							}
							else		//第四象限
							{
								D = PI - D;
							}
						}
						//printf("%lf\n", D*180/PI);
						I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

						//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
						double TransMat[9];
						double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
						TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
						TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
						TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

						//3. 根据变换矩阵计算变换后的点
						double newPos1[3], newPos2[3];
						Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
						Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
						//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
						double point_y[3]; point_y[0] = 0; point_y[1] = 1; point_y[2] = 0;
						double point_z[3]; point_z[0] = 0; point_z[1] = 0; point_z[2] = 1;
						double YDirectCos[3], ZDirectCos[3];
						Mat_Multiply(TransMat, point_y, YDirectCos, 3, 3);//北向在模型坐标系下的方向余弦
						Mat_Multiply(TransMat, point_z, ZDirectCos, 3, 3);
						//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
						//4. 将观测坐标转换到模型坐标系进行计算
						double x0 = newPos2[0];
						double y0 = newPos2[1];
						z1 = newPos1[2], z2 = newPos2[2];
						//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
						double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
						angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
						if (angle_I < 0)				//当磁化方向为向上时
						{
							angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
						}
						double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
						//********************
						//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
						//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
						double Mx = M*cos(angle_I)*sin(angle_D);
						double My = M*cos(angle_I)*cos(angle_D);
						double Mz = M*sin(angle_I);
						double north_x, east_y, x, y;
						double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
						double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
						double SurveyPoint[3], TransPoint[3];
						for (int j = 0; j < number_y; j++)
						{
							north_x = ymin + dy*j;
							for (int k = 0; k < number_x; k++)
							{
								east_y = xmin + k*dx;
								//将计算点变换到圆柱体坐标系中
								SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
								Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
								x = TransPoint[0]; y = TransPoint[1];
								xx0 = x0-x; yy0 = y0-y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
								xx02 = xx0*xx0; yy02 = yy0*yy0;
								rxy = (xx02 + yy02);
								rxy2 = rxy*rxy;
								r1 = sqrt(xx02 + yy02 + zz1*zz1);
								r2 = sqrt(xx02 + yy02 + zz2*zz2);
								//--------------------------------------------------------------------
								Vxx = (-xx02 / pow(zz2 + r2, 2.0) / r2 / r2 + (yy02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
									- (-xx02 / pow(zz1 + r1, 2.0) / r1 / r1 + (yy02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
								Vxy = (yy0*xx0 / pow(zz1 + r1, 2.0) / r1 / r1 + xx0*yy0 / (zz1 + r1) / pow(r1, 3.0))
									- (yy0*xx0 / pow(zz2 + r2, 2.0) / r2 / r2 + xx0*yy0 / (zz2 + r2) / pow(r2, 3.0));//*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
								Vxz = xx0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
								Vyy = (-yy02 / pow(zz2 + r2, 2.0) / r2 / r2 + (xx02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
									- (-yy02 / pow(zz1 + r1, 2.0) / r1 / r1 + (xx02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
								Vyz = yy0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
								Vzz = (zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));
								grav[j][k] += (ZDirectCos[0] * (Vxx*YDirectCos[0] + Vxy*YDirectCos[1] + Vxz*YDirectCos[2]) +
									ZDirectCos[1] * (Vxy*YDirectCos[0] + Vyy*YDirectCos[1] + Vyz*YDirectCos[2]) +
									ZDirectCos[2] * (Vxz*YDirectCos[0] + Vyz*YDirectCos[1] + Vzz*YDirectCos[2]))*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
							}
						}
					}
				}
			}
			break;
			case FORWARD_Vzz:
			{
				for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
				{
					if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
					{
						for (int ii = 0; ii < 3; ii++)
						{
							double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
							forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
							forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
						}
					}
					//0. 与坐标系无关的量
					double r = forwardinfo.model.cylinder_vec[i].Radius;
					double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
					double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
					double S = PI*r*r;
					double EarthMag = forwardinfo.EarthMag;
					double M = S*CHL*EarthMag;											//与系数消去了U0
					//1. 根据线模型的两个端点计算转换矩阵
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
					forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
					temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
					forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
					double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
					double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
					
					if (x1 == x2 && y1 == y2)		//直立模型
					{
						double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
						double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
						double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
						double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
						double Mx = M*cos(angle_I)*sin(angle_D);
						double My = M*cos(angle_I)*cos(angle_D);
						double Mz = M*sin(angle_I);
						double north_x, east_y, x, y;
						double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
						double Vzz;
						for (int j = 0; j < number_y; j++)
						{
							north_x = ymin + dy*j;
							//z = ymin+dy*j;//计算沿东向的切片
							for (int k = 0; k < number_x; k++)
							{
								east_y = xmin + k*dx;
								//z = xmin +dx*k;//沿北向的切片
								//将计算点变换到圆柱体坐标系中
								x = north_x; y = east_y;
								xx0 = x0-x; yy0 = y0-y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
								xx02 = xx0*xx0; yy02 = yy0*yy0;
								rxy = (xx02 + yy02);
								rxy2 = rxy*rxy;
								r1 = sqrt(xx02 + yy02 + zz1*zz1);
								r2 = sqrt(xx02 + yy02 + zz2*zz2);
								//--------------------------------------------------------------------
								Vzz = (zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));
								Vzz = Vzz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
								grav[j][k] += Vzz;
							}
						}
					}
					else			//非直立模型
					{
						//1. 圆柱体的计算倾角和偏角
						double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
						double L, L_H;
						L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
						L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
						direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
						direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
						double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
						cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
						//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
						dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
						//printf("%lf\n",dot_horizontal);
						D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
						if (dot_horizontal < 0)
						{
							if (D < 0)	//第三象限
							{
								D = -(D + PI);
							}
							else		//第四象限
							{
								D = PI - D;
							}
						}
						//printf("%lf\n", D*180/PI);
						I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

						//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
						double TransMat[9];
						double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
						TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
						TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
						TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

						//3. 根据变换矩阵计算变换后的点
						double newPos1[3], newPos2[3];
						Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
						Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
						//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
						double point_z[3]; point_z[0] = 0; point_z[1] = 0; point_z[2] = 1;
						double DirectCos[3];
						Mat_Multiply(TransMat, point_z, DirectCos, 3, 3);//北向在模型坐标系下的方向余弦
						//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
						//4. 将观测坐标转换到模型坐标系进行计算
						double x0 = newPos2[0];
						double y0 = newPos2[1];
						z1 = newPos1[2], z2 = newPos2[2];
						//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
						double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
						angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
						if (angle_I < 0)				//当磁化方向为向上时
						{
							angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
						}
						double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
						//********************
						//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
						//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
						double Mx = M*cos(angle_I)*sin(angle_D);
						double My = M*cos(angle_I)*cos(angle_D);
						double Mz = M*sin(angle_I);
						double north_x, east_y, x, y;
						double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
						double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
						double SurveyPoint[3], TransPoint[3];
						for (int j = 0; j < number_y; j++)
						{
							north_x = ymin + dy*j;
							for (int k = 0; k < number_x; k++)
							{
								east_y = xmin + k*dx;
								//将计算点变换到圆柱体坐标系中
								SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
								Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
								x = TransPoint[0]; y = TransPoint[1];
								xx0 = x0-x; yy0 = y0-y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
								xx02 = xx0*xx0; yy02 = yy0*yy0;
								rxy = (xx02 + yy02);
								rxy2 = rxy*rxy;
								r1 = sqrt(xx02 + yy02 + zz1*zz1);
								r2 = sqrt(xx02 + yy02 + zz2*zz2);
								//--------------------------------------------------------------------
								Vxx = (-xx02 / pow(zz2 + r2, 2.0) / r2 / r2 + (yy02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
									- (-xx02 / pow(zz1 + r1, 2.0) / r1 / r1 + (yy02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
								Vxy = (yy0*xx0 / pow(zz1 + r1, 2.0) / r1 / r1 + xx0*yy0 / (zz1 + r1) / pow(r1, 3.0))
									- (yy0*xx0 / pow(zz2 + r2, 2.0) / r2 / r2 + xx0*yy0 / (zz2 + r2) / pow(r2, 3.0));//*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
								Vxz = xx0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
								Vyy = (-yy02 / pow(zz2 + r2, 2.0) / r2 / r2 + (xx02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
									- (-yy02 / pow(zz1 + r1, 2.0) / r1 / r1 + (xx02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
								Vyz = yy0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
								Vzz = (zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));
								grav[j][k] += (DirectCos[0] * (Vxx*DirectCos[0] + Vxy*DirectCos[1] + Vxz*DirectCos[2]) +
									DirectCos[1] * (Vxy*DirectCos[0] + Vyy*DirectCos[1] + Vyz*DirectCos[2]) +
									DirectCos[2] * (Vxz*DirectCos[0] + Vyz*DirectCos[1] + Vzz*DirectCos[2]))*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
							}
						}
					}
				}
			}
			break;
		default:
			MessageBox(NULL, _T("这是重力位及其导数正演函数，请输入正确的type"), _T("有限长倾斜线模型正演-错误提示"), MB_OK);
			return 0;

	}
	GetGrdMinMax(grav, forwardinfo0.model.grddatainfo);
	return 1;
}

int _3DFiniteLine_mag(double** mag, RegularGeometry3DForward& forwardinfo0, int type)
{
	RegularGeometry3DForward forwardinfo = forwardinfo0;
	//首先将mag赋值为0，后面累加
	int number_x = forwardinfo.model.grddatainfo.AutoGetNumber_x();
	int number_y = forwardinfo.model.grddatainfo.AutoGetNumber_y();
	Assign_Array2(mag, number_y, number_x, 0);
	//坐标范围
	double xmin = forwardinfo.model.grddatainfo.m_AxisBounds[0];
	double xmax = forwardinfo.model.grddatainfo.m_AxisBounds[1];
	double ymin = forwardinfo.model.grddatainfo.m_AxisBounds[2];
	double ymax = forwardinfo.model.grddatainfo.m_AxisBounds[3];
	double dx = forwardinfo.model.grddatainfo.m_Dx;
	double dy = forwardinfo.model.grddatainfo.m_Dy;
	switch (type)
	{
		case FORWARD_Hax:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
				{
					for (int ii = 0; ii < 3; ii++)
					{
						double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
						forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
						forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
					}
				}
				//0. 与坐标系无关的量
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double S = PI*r*r;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				//1. 根据线模型的两个端点计算转换矩阵
				double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
				forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
				temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
				forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
				double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
				
				if (x1 == x2 && y1 == y2)		//直立模型
				{
					double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
					double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vxy,Vxx,Vxz;
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						//z = ymin+dy*j;//计算沿东向的切片
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//z = xmin +dx*k;//沿北向的切片
							//将计算点变换到圆柱体坐标系中
							x = north_x; y = east_y;
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							//Vxx = (yy02 - xx02) / rxy2*(zz1 / r1 - zz2 / r2) - xx02 / rxy*(zz1/pow(r1,3.0)-zz2/pow(r2,3.0));//化简的公式
							Vxx = (-xx02 / pow(zz2 + r2, 2.0) / r2 / r2 + (yy02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
								- (-xx02 / pow(zz1 + r1, 2.0) / r1 / r1 + (yy02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
							//Vxx = Vxx*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

							//Vxy = -2 * xx0*yy0 / rxy2*(zz1 / r1-zz2 / r2) - xx0*yy0 / rxy*(zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));//化简的公式
							Vxy = (yy0*xx0 / pow(zz1 + r1, 2.0) / r1 / r1 + xx0*yy0 / (zz1 + r1) / pow(r1, 3.0))
								- (yy0*xx0 / pow(zz2 + r2, 2.0) / r2 / r2 + xx0*yy0 / (zz2 + r2) / pow(r2, 3.0));
							//Vxy = Vxy*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E					//有问题

							Vxz = xx0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0)); //Vxz = Vxz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

							mag[j][k] += 1.0 / 4.0 / PI*(Mx*Vxx + My*Vxy + Mz*Vxz);
						}
					}
				}
				else			//非直立模型
				{
					//1. 圆柱体的计算倾角和偏角
					double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
					double L, L_H;
					L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
					L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
					direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
					direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
					double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
					cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
					//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
					dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
					//printf("%lf\n",dot_horizontal);
					D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
					if (dot_horizontal < 0)
					{
						if (D < 0)	//第三象限
						{
							D = -(D + PI);
						}
						else		//第四象限
						{
							D = PI - D;
						}
					}
					//printf("%lf\n", D*180/PI);
					I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

					//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
					double TransMat[9];
					double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
					TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
					TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
					TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

					//3. 根据变换矩阵计算变换后的点
					double newPos1[3], newPos2[3];
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
					//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
					double point_x[3]; point_x[0] = 1; point_x[1] = 0; point_x[2] = 0;
					double point_y[3]; point_y[0] = 0; point_y[1] = 1; point_y[2] = 0;
					double point_z[3]; point_z[0] = 0; point_z[1] = 0; point_z[2] = 1;
					double XDirectCos[3], YDirectCos[3], ZDirectCos[3];
					Mat_Multiply(TransMat, point_x, XDirectCos, 3, 3);//北向在模型坐标系下的方向余弦
					Mat_Multiply(TransMat, point_y, YDirectCos, 3, 3);
					Mat_Multiply(TransMat, point_z, ZDirectCos, 3, 3);
					//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
					//4. 将观测坐标转换到模型坐标系进行计算
					double x0 = newPos2[0];
					double y0 = newPos2[1];
					z1 = newPos1[2], z2 = newPos2[2];
					//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
					//********************
					//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
					//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
					double Haxx, Haxy, Haxz;
					double SurveyPoint[3], TransPoint[3];
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//将计算点变换到圆柱体坐标系中
							SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
							Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
							x = TransPoint[0]; y = TransPoint[1];
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							Vxx = (-xx02 / pow(zz2 + r2, 2.0) / r2 / r2 + (yy02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
								- (-xx02 / pow(zz1 + r1, 2.0) / r1 / r1 + (yy02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
							Vxy = (yy0*xx0 / pow(zz1 + r1, 2.0) / r1 / r1 + xx0*yy0 / (zz1 + r1) / pow(r1, 3.0))
								- (yy0*xx0 / pow(zz2 + r2, 2.0) / r2 / r2 + xx0*yy0 / (zz2 + r2) / pow(r2, 3.0));//*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
							Vxz = xx0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
							Vyy = (-yy02 / pow(zz2 + r2, 2.0) / r2 / r2 + (xx02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
								- (-yy02 / pow(zz1 + r1, 2.0) / r1 / r1 + (xx02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
							Vyz = yy0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
							Vzz = (zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));
							Haxx = (XDirectCos[0] * (Vxx*XDirectCos[0] + Vxy*XDirectCos[1] + Vxz*XDirectCos[2]) +
								XDirectCos[1] * (Vxy*XDirectCos[0] + Vyy*XDirectCos[1] + Vyz*XDirectCos[2]) +
								XDirectCos[2] * (Vxz*XDirectCos[0] + Vyz*XDirectCos[1] + Vzz*XDirectCos[2]));//E
							Haxy = (YDirectCos[0] * (Vxx*XDirectCos[0] + Vxy*XDirectCos[1] + Vxz*XDirectCos[2]) +
								YDirectCos[1] * (Vxy*XDirectCos[0] + Vyy*XDirectCos[1] + Vyz*XDirectCos[2]) +
								YDirectCos[2] * (Vxz*XDirectCos[0] + Vyz*XDirectCos[1] + Vzz*XDirectCos[2]));//E
							Haxz = (ZDirectCos[0] * (Vxx*XDirectCos[0] + Vxy*XDirectCos[1] + Vxz*XDirectCos[2]) +
								ZDirectCos[1] * (Vxy*XDirectCos[0] + Vyy*XDirectCos[1] + Vyz*XDirectCos[2]) +
								ZDirectCos[2] * (Vxz*XDirectCos[0] + Vyz*XDirectCos[1] + Vzz*XDirectCos[2]));//E

							mag[j][k] += 1.0 / 4.0 / PI*(Mx*Haxx + My*Haxy + Mz*Haxz);//nT
						}
					}
				}
			}
		}
		break;
		case FORWARD_Hay:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
				{
					for (int ii = 0; ii < 3; ii++)
					{
						double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
						forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
						forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
					}
				}
				//0. 与坐标系无关的量
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double S = PI*r*r;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				//1. 根据线模型的两个端点计算转换矩阵
				double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
				forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
				temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
				forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
				double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
				
				if (x1 == x2 && y1 == y2)		//直立模型
				{
					double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
					double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vxy, Vyy, Vyz;
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						//z = ymin+dy*j;//计算沿东向的切片
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//z = xmin +dx*k;//沿北向的切片
							//将计算点变换到圆柱体坐标系中
							x = north_x; y = east_y;
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							//Vyy = (xx02-yy02) / rxy2*(zz1 / r1 - zz2 / r2) - yy02 / rxy*(zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));//化简的公式
							Vyy = (-yy02 / pow(zz2 + r2, 2.0) / r2 / r2 + (xx02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
								- (-yy02 / pow(zz1 + r1, 2.0) / r1 / r1 + (xx02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));

							//Vxy = -2 * xx0*yy0 / rxy2*(zz1 / r1-zz2 / r2) - xx0*yy0 / rxy*(zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));//化简的公式
							Vxy = (yy0*xx0 / pow(zz1 + r1, 2.0) / r1 / r1 + xx0*yy0 / (zz1 + r1) / pow(r1, 3.0))
								- (yy0*xx0 / pow(zz2 + r2, 2.0) / r2 / r2 + xx0*yy0 / (zz2 + r2) / pow(r2, 3.0));

							Vyz = yy0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0)); //Vyz = Vyz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

							mag[j][k] += 1.0 / 4.0 / PI*(Mx*Vxy + My*Vyy + Mz*Vyz);
						}
					}
				}
				else			//非直立模型
				{
					//1. 圆柱体的计算倾角和偏角
					double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
					double L, L_H;
					L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
					L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
					direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
					direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
					double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
					cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
					//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
					dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
					//printf("%lf\n",dot_horizontal);
					D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
					if (dot_horizontal < 0)
					{
						if (D < 0)	//第三象限
						{
							D = -(D + PI);
						}
						else		//第四象限
						{
							D = PI - D;
						}
					}
					//printf("%lf\n", D*180/PI);
					I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

					//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
					double TransMat[9];
					double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
					TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
					TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
					TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

					//3. 根据变换矩阵计算变换后的点
					double newPos1[3], newPos2[3];
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
					//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
					double point_x[3]; point_x[0] = 1; point_x[1] = 0; point_x[2] = 0;
					double point_y[3]; point_y[0] = 0; point_y[1] = 1; point_y[2] = 0;
					double point_z[3]; point_z[0] = 0; point_z[1] = 0; point_z[2] = 1;
					double XDirectCos[3], YDirectCos[3], ZDirectCos[3];
					Mat_Multiply(TransMat, point_x, XDirectCos, 3, 3);//北向在模型坐标系下的方向余弦
					Mat_Multiply(TransMat, point_y, YDirectCos, 3, 3);
					Mat_Multiply(TransMat, point_z, ZDirectCos, 3, 3);
					//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
					//4. 将观测坐标转换到模型坐标系进行计算
					double x0 = newPos2[0];
					double y0 = newPos2[1];
					z1 = newPos1[2], z2 = newPos2[2];
					//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
					//********************
					//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
					//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
					double Haxy, Hayy, Hazy;
					double SurveyPoint[3], TransPoint[3];
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//将计算点变换到圆柱体坐标系中
							SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
							Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
							x = TransPoint[0]; y = TransPoint[1];
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							Vxx = (-xx02 / pow(zz2 + r2, 2.0) / r2 / r2 + (yy02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
								- (-xx02 / pow(zz1 + r1, 2.0) / r1 / r1 + (yy02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
							Vxy = (yy0*xx0 / pow(zz1 + r1, 2.0) / r1 / r1 + xx0*yy0 / (zz1 + r1) / pow(r1, 3.0))
								- (yy0*xx0 / pow(zz2 + r2, 2.0) / r2 / r2 + xx0*yy0 / (zz2 + r2) / pow(r2, 3.0));//*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
							Vxz = xx0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
							Vyy = (-yy02 / pow(zz2 + r2, 2.0) / r2 / r2 + (xx02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
								- (-yy02 / pow(zz1 + r1, 2.0) / r1 / r1 + (xx02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
							Vyz = yy0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
							Vzz = (zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));
							Haxy = (XDirectCos[0] * (Vxx*YDirectCos[0] + Vxy*YDirectCos[1] + Vxz*YDirectCos[2]) +
								XDirectCos[1] * (Vxy*YDirectCos[0] + Vyy*YDirectCos[1] + Vyz*YDirectCos[2]) +
								XDirectCos[2] * (Vxz*YDirectCos[0] + Vyz*YDirectCos[1] + Vzz*YDirectCos[2]));//E
							Hayy = (YDirectCos[0] * (Vxx*YDirectCos[0] + Vxy*YDirectCos[1] + Vxz*YDirectCos[2]) +
								YDirectCos[1] * (Vxy*YDirectCos[0] + Vyy*YDirectCos[1] + Vyz*YDirectCos[2]) +
								YDirectCos[2] * (Vxz*YDirectCos[0] + Vyz*YDirectCos[1] + Vzz*YDirectCos[2]));//E
							Hazy = (ZDirectCos[0] * (Vxx*YDirectCos[0] + Vxy*YDirectCos[1] + Vxz*YDirectCos[2]) +
								ZDirectCos[1] * (Vxy*YDirectCos[0] + Vyy*YDirectCos[1] + Vyz*YDirectCos[2]) +
								ZDirectCos[2] * (Vxz*YDirectCos[0] + Vyz*YDirectCos[1] + Vzz*YDirectCos[2]));//E

							mag[j][k] += 1.0 / 4.0 / PI*(Mx*Haxy + My*Hayy + Mz*Hazy);//nT
						}
					}
				}
			}
		}
		break;
		case FORWARD_Za:
		{
			for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
			{
				if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
				{
					for (int ii = 0; ii < 3; ii++)
					{
						double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
						forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
						forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
					}
				}
				//0. 与坐标系无关的量
				double r = forwardinfo.model.cylinder_vec[i].Radius;
				double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
				double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
				double S = PI*r*r;
				double EarthMag = forwardinfo.EarthMag;
				double M = S*CHL*EarthMag;											//与系数消去了U0
				//1. 根据线模型的两个端点计算转换矩阵
				double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
				forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
				temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
				forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
				double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
				double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
				
				if (x1 == x2 && y1 == y2)		//直立模型
				{
					double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
					double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vxz, Vyz, Vzz;
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						//z = ymin+dy*j;//计算沿东向的切片
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//z = xmin +dx*k;//沿北向的切片
							//将计算点变换到圆柱体坐标系中
							x = north_x; y = east_y;
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							Vxz = xx0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0)); //Vxz = Vxz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

							Vyz = yy0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));// Vyz = Vyz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

							Vzz = (zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0)); //Vzz = Vzz*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E

							mag[j][k] += 1.0 / 4.0 / PI*(Mx*Vxz + My*Vyz + Mz*Vzz);
						}
					}
				}
				else			//非直立模型
				{
					//1. 圆柱体的计算倾角和偏角
					double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
					double L, L_H;
					L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
					L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
					direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
					direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
					double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
					cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
					//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
					dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
					//printf("%lf\n",dot_horizontal);
					D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
					if (dot_horizontal < 0)
					{
						if (D < 0)	//第三象限
						{
							D = -(D + PI);
						}
						else		//第四象限
						{
							D = PI - D;
						}
					}
					//printf("%lf\n", D*180/PI);
					I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

					//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
					double TransMat[9];
					double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
					TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
					TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
					TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

					//3. 根据变换矩阵计算变换后的点
					double newPos1[3], newPos2[3];
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
					Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
					//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
					double point_x[3]; point_x[0] = 1; point_x[1] = 0; point_x[2] = 0;
					double point_y[3]; point_y[0] = 0; point_y[1] = 1; point_y[2] = 0;
					double point_z[3]; point_z[0] = 0; point_z[1] = 0; point_z[2] = 1;
					double XDirectCos[3], YDirectCos[3], ZDirectCos[3];
					Mat_Multiply(TransMat, point_x, XDirectCos, 3, 3);//北向在模型坐标系下的方向余弦
					Mat_Multiply(TransMat, point_y, YDirectCos, 3, 3);
					Mat_Multiply(TransMat, point_z, ZDirectCos, 3, 3);
					//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
					//4. 将观测坐标转换到模型坐标系进行计算
					double x0 = newPos2[0];
					double y0 = newPos2[1];
					z1 = newPos1[2], z2 = newPos2[2];
					//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
					double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
					double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
					//********************
					//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
					//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
					double Mx = M*cos(angle_I)*sin(angle_D);
					double My = M*cos(angle_I)*cos(angle_D);
					double Mz = M*sin(angle_I);
					double north_x, east_y, x, y;
					double r1, r2, xx0, yy0, zz1, zz2, xx02, yy02, rxy2, rxy;
					double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
					double Haxz, Hayz, Hazz;
					double SurveyPoint[3], TransPoint[3];
					for (int j = 0; j < number_y; j++)
					{
						north_x = ymin + dy*j;
						for (int k = 0; k < number_x; k++)
						{
							east_y = xmin + k*dx;
							//将计算点变换到圆柱体坐标系中
							SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
							Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
							x = TransPoint[0]; y = TransPoint[1];
							xx0 = x0-x; yy0 = y0-y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
							xx02 = xx0*xx0; yy02 = yy0*yy0;
							rxy = (xx02 + yy02);
							rxy2 = rxy*rxy;
							r1 = sqrt(xx02 + yy02 + zz1*zz1);
							r2 = sqrt(xx02 + yy02 + zz2*zz2);
							//--------------------------------------------------------------------
							Vxx = (-xx02 / pow(zz2 + r2, 2.0) / r2 / r2 + (yy02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
								- (-xx02 / pow(zz1 + r1, 2.0) / r1 / r1 + (yy02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
							Vxy = (yy0*xx0 / pow(zz1 + r1, 2.0) / r1 / r1 + xx0*yy0 / (zz1 + r1) / pow(r1, 3.0))
								- (yy0*xx0 / pow(zz2 + r2, 2.0) / r2 / r2 + xx0*yy0 / (zz2 + r2) / pow(r2, 3.0));//*G*forwardinfo.model.cylinder_vec[i].Density*S*1E12;//E
							Vxz = xx0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
							Vyy = (-yy02 / pow(zz2 + r2, 2.0) / r2 / r2 + (xx02 + zz2*zz2) / (zz2 + r2) / pow(r2, 3.0))
								- (-yy02 / pow(zz1 + r1, 2.0) / r1 / r1 + (xx02 + zz1*zz1) / (zz1 + r1) / pow(r1, 3.0));
							Vyz = yy0*(1.0 / pow(r2, 3.0) - 1.0 / pow(r1, 3.0));
							Vzz = (zz1 / pow(r1, 3.0) - zz2 / pow(r2, 3.0));
							Haxz = (XDirectCos[0] * (Vxx*ZDirectCos[0] + Vxy*ZDirectCos[1] + Vxz*ZDirectCos[2]) +
								XDirectCos[1] * (Vxy*ZDirectCos[0] + Vyy*ZDirectCos[1] + Vyz*ZDirectCos[2]) +
								XDirectCos[2] * (Vxz*ZDirectCos[0] + Vyz*ZDirectCos[1] + Vzz*ZDirectCos[2]));//E
							Hayz = (YDirectCos[0] * (Vxx*ZDirectCos[0] + Vxy*ZDirectCos[1] + Vxz*ZDirectCos[2]) +
								YDirectCos[1] * (Vxy*ZDirectCos[0] + Vyy*ZDirectCos[1] + Vyz*ZDirectCos[2]) +
								YDirectCos[2] * (Vxz*ZDirectCos[0] + Vyz*ZDirectCos[1] + Vzz*ZDirectCos[2]));//E
							Hazz = (ZDirectCos[0] * (Vxx*ZDirectCos[0] + Vxy*ZDirectCos[1] + Vxz*ZDirectCos[2]) +
								ZDirectCos[1] * (Vxy*ZDirectCos[0] + Vyy*ZDirectCos[1] + Vyz*ZDirectCos[2]) +
								ZDirectCos[2] * (Vxz*ZDirectCos[0] + Vyz*ZDirectCos[1] + Vzz*ZDirectCos[2]));//E

							mag[j][k] += 1.0 / 4.0 / PI*(Mx*Haxz + My*Hayz + Mz*Hazz);//nT
						}
					}
				}
			}
		}
		break;
	default:
		MessageBox(NULL, _T("这是磁异常三分量正演函数，请输入正确的type"), _T("有限长倾斜线模型正演-错误提示"), MB_OK);
		return 0;

	}
	GetGrdMinMax(mag, forwardinfo0.model.grddatainfo);
	return 1;
}

int _3DFiniteLine_magT(double** mag, RegularGeometry3DForward& forwardinfo0, int type)
{
	RegularGeometry3DForward forwardinfo = forwardinfo0;
	//首先将grav赋值为0，后面累加
	int number_x = forwardinfo.model.grddatainfo.AutoGetNumber_x();
	int number_y = forwardinfo.model.grddatainfo.AutoGetNumber_y();
	double** Hax = CreateArray2(number_y, number_x);
	double** Hay = CreateArray2(number_y, number_x);
	double** Za = CreateArray2(number_y, number_x);
	Assign_Array2(Hax, number_y, number_x, 0);
	Assign_Array2(Hay, number_y, number_x, 0);
	Assign_Array2(Za, number_y, number_x, 0);
	Assign_Array2(mag, number_y, number_x, 0);
	double I = forwardinfo.EarthAngle_I / 180.0*PI;
	double D = forwardinfo.EarthAngle_D / 180.0*PI;

	//计算Hax
	_3DFiniteLine_mag(Hax, forwardinfo, FORWARD_Hax);
	//计算Hay
	_3DFiniteLine_mag(Hay, forwardinfo, FORWARD_Hay);
	//计算Za
	_3DFiniteLine_mag(Za, forwardinfo, FORWARD_Za);
	switch (type)
	{
		case FORWARD_Ta:
		{
			//计算Ta
			for (int i = 0; i < number_y; i++)
			{
				for (int j = 0; j < number_x; j++)
				{
					mag[i][j] = Hax[i][j] * cos(I)*cos(D) + Hay[i][j] * cos(I)*sin(D) + Za[i][j] * sin(I);
				}
			}
		}
		break;
		case FORWARD_Module:
		{
			//计算模量
			for (int i = 0; i < number_y; i++)
			{
				for (int j = 0; j < number_x; j++)
				{
					mag[i][j] = sqrt(Hax[i][j] * Hax[i][j] + Hay[i][j] * Hay[i][j] + Za[i][j] * Za[i][j]);
				}
			}
		}
		break;
		default:
			MessageBox(NULL, _T("这是Ta和模量正演，输入正确的类型参数"), _T("有限长线模型-出错提示"), MB_OK);
		return 0;
	}

	GetGrdMinMax(mag, forwardinfo0.model.grddatainfo);
	//销毁二维数组
	DeleteArray2(Hax, number_y, number_x);
	DeleteArray2(Hay, number_y, number_x);
	DeleteArray2(Za, number_y, number_x);

	return 0;
}

int _3DFiniteVercicalCylinder_grav(double** grav, RegularGeometry3DForward& forwardinfo, int type, int theta_num)
{
	//首先将mag赋值为0，后面累加
	int number_x = forwardinfo.model.grddatainfo.AutoGetNumber_x();
	int number_y = forwardinfo.model.grddatainfo.AutoGetNumber_y();
	Assign_Array2(grav, number_y, number_x, 0);
	//坐标范围
	double xmin = forwardinfo.model.grddatainfo.m_AxisBounds[0];
	double xmax = forwardinfo.model.grddatainfo.m_AxisBounds[1];
	double ymin = forwardinfo.model.grddatainfo.m_AxisBounds[2];
	double ymax = forwardinfo.model.grddatainfo.m_AxisBounds[3];
	double dx = forwardinfo.model.grddatainfo.m_Dx;
	double dy = forwardinfo.model.grddatainfo.m_Dy;
	switch (type)
	{
	case FORWARD_V:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
			double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos2[1] != forwardinfo.model.cylinder_vec[i].Pos1[1] && forwardinfo.model.cylinder_vec[i].Pos2[0] != forwardinfo.model.cylinder_vec[i].Pos1[0])
			{
				CString strinfo; strinfo.Format(_T("第%d个模型不是直立圆柱体"), i);
				MessageBox(NULL, strinfo, _T("有限长直立圆柱体正演-错误提示"), MB_OK);
				return 0;
			}
			double x0 = forwardinfo.model.cylinder_vec[i].Pos2[1];//带入公式中计算的x，y其实并不是surfer绘图中的x，y而是值北向坐标和东向坐标
			double y0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
			//==================================================上面是平面旋转的代码=====================================================================
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double S = PI*r*r;
			double H = (z1 - z2);
			double north_x, east_y, x, y;
			double rr2, xx0, yy0, zz1, zz2, xx02, yy02, temp1, fenzi, fenmu, temp21, temp22;
			double V;
			double dtheta = 2 * PI / theta_num;
			double theta;
			//north_x = 140;//东向切面
			//east_y = 90;//北向切面
			for (int j = 0; j < number_y; j++)
			{
				north_x = ymin + dy*j;
				//z = ymin + dy*j;//计算沿东向的切片
				for (int k = 0; k < number_x; k++)
				{
					east_y = xmin + k*dx;
					//z = xmin +dx*k;//沿北向的切片
					//将计算点变换到圆柱体坐标系中
					x = north_x; y = east_y;
					xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
					rr2 = r*r;
					xx02 = xx0*xx0; yy02 = yy0*yy0;
					V = 0;
					for (int kk = 0; kk < theta_num; kk++)//计算积分
					{
						theta = kk*dtheta;

						temp1 = xx0*cos(theta) + yy0*sin(theta);
						fenzi = rr2 + r*temp1 + 1E-10;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
						fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + 1E-10;
						/*if (fenmu == 0 && fenzi == 0)
						{
						printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

						}*/
						temp21 = sqrt(fenmu + zz1*zz1);
						temp22 = sqrt(fenmu + zz2*zz2);
						//V = V + fenzi/fenmu*(-zz1 / 2.0 * temp21 + fenmu / 2.0 * log(-zz1 + temp21) + zz2 / 2 * temp22 - fenmu / 2 * log(-zz2 + temp22) + zz1 / 2 * abs(zz1) - zz2 / 2 * abs(zz2))*dtheta;
						V = V + fenzi / fenmu/2.0*(temp22*zz2-temp21*zz1+fenmu*log((-zz1+temp21)/(-zz2+temp22))+zz1*fabs(zz1)-zz2*fabs(zz2))*dtheta;
					}
					//--------------------------------------------------------------------
					V = V*G*forwardinfo.model.cylinder_vec[i].Density*1E8;
					grav[j][k] += V;			//
				}
			}
		}
	}
		break;
	case FORWARD_Vz:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
			double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos2[1] != forwardinfo.model.cylinder_vec[i].Pos1[1] && forwardinfo.model.cylinder_vec[i].Pos2[0] != forwardinfo.model.cylinder_vec[i].Pos1[0])
			{
				CString strinfo; strinfo.Format(_T("第%d个模型不是直立圆柱体"),i);
				MessageBox(NULL, strinfo, _T("有限长直立圆柱体正演-错误提示"), MB_OK);
				return 0;
			}
			double x0 = forwardinfo.model.cylinder_vec[i].Pos2[1];//带入公式中计算的x，y其实并不是surfer绘图中的x，y而是值北向坐标和东向坐标
			double y0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
			//==================================================上面是平面旋转的代码=====================================================================
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double S = PI*r*r;
			double H = (z1 - z2);
			double north_x, east_y, x, y;
			double rr2, xx0, yy0, zz1, zz2, xx02, yy02,temp1,fenzi,fenmu,temp21,temp22;
			double Vz;
			double dtheta =2.0 * PI / theta_num;
			double theta;
			//north_x = 140;//东向切面
			//east_y = 90;//北向切面
			for (int j = 0; j < number_y; j++)
			{
				north_x = ymin + dy*j;
				//z = ymin + dy*j;//计算沿东向的切片
				for (int k = 0; k < number_x; k++)
				{
					east_y = xmin + k*dx;
					//z = xmin +dx*k;//沿北向的切片
					//将计算点变换到圆柱体坐标系中
					x = north_x; y = east_y;
					xx0 = x0-x; yy0 = y0-y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
					rr2 = r*r;
					xx02 = xx0*xx0; yy02 = yy0*yy0;
					Vz = 0;
					for (int kk = 0; kk < theta_num; kk++)//计算积分
					{
						theta = kk*dtheta;

						temp1 = xx0*cos(theta) + yy0*sin(theta);
						fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
						fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
						temp21 = sqrt(fenmu + zz1*zz1);
						temp22 = sqrt(fenmu + zz2*zz2);
						Vz += fenzi / fenmu*(temp21 - temp22 - fabs(zz1) + fabs(zz2))*dtheta;
					}
					//--------------------------------------------------------------------
					Vz = Vz*G*forwardinfo.model.cylinder_vec[i].Density*1E8;
					grav[j][k] += Vz;			//
				}
			}
		}
	}
		break;
	case FORWARD_Vzz:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
			double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos2[1] != forwardinfo.model.cylinder_vec[i].Pos1[1] && forwardinfo.model.cylinder_vec[i].Pos2[0] != forwardinfo.model.cylinder_vec[i].Pos1[0])
			{
				CString strinfo; strinfo.Format(_T("第%d个模型不是直立圆柱体"), i);
				MessageBox(NULL, strinfo, _T("有限长直立圆柱体正演-错误提示"), MB_OK);
				return 0;
			}
			double x0 = forwardinfo.model.cylinder_vec[i].Pos2[1];//带入公式中计算的x，y其实并不是surfer绘图中的x，y而是值北向坐标和东向坐标
			double y0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
			//==================================================上面是平面旋转的代码=====================================================================
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double S = PI*r*r;
			double H = (z1 - z2);
			double north_x, east_y, x, y;
			double rr2, xx0, yy0, zz1, zz2, xx02, yy02, temp1, fenzi, fenmu, temp21, temp22;
			double Vzz;
			double dtheta = 2 * PI / theta_num;
			double theta;
			//north_x = 140;//东向切面
			//east_y = 90;//北向切面
			for (int j = 0; j < number_y; j++)
			{
				north_x = ymin + dy*j;
				//z = ymin + dy*j;//计算沿东向的切片
				for (int k = 0; k < number_x; k++)
				{
					east_y = xmin + k*dx;
					//z = xmin +dx*k;//沿北向的切片
					//将计算点变换到圆柱体坐标系中
					x = north_x; y = east_y;
					xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
					rr2 = r*r;
					xx02 = xx0*xx0; yy02 = yy0*yy0;
					Vzz = 0;
					for (int kk = 0; kk < theta_num; kk++)//计算积分
					{
						theta = kk*dtheta;

						temp1 = xx0*cos(theta) + yy0*sin(theta);
						fenzi = rr2 + r*temp1 + 1E-10;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
						fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + 1E-10;
						/*if (fenmu == 0 && fenzi == 0)
						{
						printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

						}*/
						temp21 = sqrt(fenmu + zz1*zz1);
						temp22 = sqrt(fenmu + zz2*zz2);
						Vzz += fenzi / fenmu*(zz2/temp22-zz1/temp21)*dtheta;
						//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
					}
					//--------------------------------------------------------------------
					Vzz = Vzz*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
					grav[j][k] += Vzz;			//
				}
			}
		}
	}
		break;
	case FORWARD_Vxz:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
			double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos2[1] != forwardinfo.model.cylinder_vec[i].Pos1[1] && forwardinfo.model.cylinder_vec[i].Pos2[0] != forwardinfo.model.cylinder_vec[i].Pos1[0])
			{
				CString strinfo; strinfo.Format(_T("第%d个模型不是直立圆柱体"), i);
				MessageBox(NULL, strinfo, _T("有限长直立圆柱体正演-错误提示"), MB_OK);
				return 0;
			}
			double x0 = forwardinfo.model.cylinder_vec[i].Pos2[1];//带入公式中计算的x，y其实并不是surfer绘图中的x，y而是值北向坐标和东向坐标
			double y0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
			//==================================================上面是平面旋转的代码=====================================================================
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double S = PI*r*r;
			double H = (z1 - z2);
			double north_x, east_y, x, y;
			double rr2, xx0, yy0, zz1, zz2, xx02, yy02, L1,L2,vpt,temp1,fenmu;
			double Vxz;
			double dtheta = 2 * PI / theta_num;
			double theta;
			//north_x = 140;//东向切面
			//east_y = 90;//北向切面
			for (int j = 0; j < number_y; j++)
			{
				north_x = ymin + dy*j;
				//z = ymin + dy*j;//计算沿东向的切片
				for (int k = 0; k < number_x; k++)
				{
					east_y = xmin + k*dx;
					//z = xmin +dx*k;//沿北向的切片
					//将计算点变换到圆柱体坐标系中
					x = north_x; y = east_y;
					xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
					rr2 = r*r;
					xx02 = xx0*xx0; yy02 = yy0*yy0;
					Vxz = 0;
					for (int kk = 0; kk < theta_num; kk++)//计算积分
					{
						theta = kk*dtheta;

						temp1 = xx0*cos(theta) + yy0*sin(theta);
						fenmu = rr2 + xx02 + yy02 + 2 * r*temp1;
						L1 = sqrt(fenmu + zz1*zz1);
						L2 = sqrt(fenmu + zz2*zz2);
						vpt = r*cos(theta);
						Vxz += (1.0/L2-1.0/L1)*vpt*dtheta;
						//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
					}
					//--------------------------------------------------------------------
					Vxz = Vxz*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
					grav[j][k] += Vxz;			//
				}
			}
		}
	}
		break;
	case FORWARD_Vyz:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
			double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos2[1] != forwardinfo.model.cylinder_vec[i].Pos1[1] && forwardinfo.model.cylinder_vec[i].Pos2[0] != forwardinfo.model.cylinder_vec[i].Pos1[0])
			{
				CString strinfo; strinfo.Format(_T("第%d个模型不是直立圆柱体"), i);
				MessageBox(NULL, strinfo, _T("有限长直立圆柱体正演-错误提示"), MB_OK);
				return 0;
			}
			double x0 = forwardinfo.model.cylinder_vec[i].Pos2[1];//带入公式中计算的x，y其实并不是surfer绘图中的x，y而是值北向坐标和东向坐标
			double y0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
			//==================================================上面是平面旋转的代码=====================================================================
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double S = PI*r*r;
			double H = (z1 - z2);
			double north_x, east_y, x, y;
			double rr2, xx0, yy0, zz1, zz2, xx02, yy02, L1, L2, upt, temp1, fenmu;
			double Vyz;
			double dtheta = 2 * PI / theta_num;
			double theta;
			//north_x = 140;//东向切面
			//east_y = 90;//北向切面
			for (int j = 0; j < number_y; j++)
			{
				north_x = ymin + dy*j;
				//z = ymin + dy*j;//计算沿东向的切片
				for (int k = 0; k < number_x; k++)
				{
					east_y = xmin + k*dx;
					//z = xmin +dx*k;//沿北向的切片
					//将计算点变换到圆柱体坐标系中
					x = north_x; y = east_y;
					xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
					rr2 = r*r;
					xx02 = xx0*xx0; yy02 = yy0*yy0;
					Vyz = 0;
					for (int kk = 0; kk < theta_num; kk++)//计算积分
					{
						theta = kk*dtheta;

						temp1 = xx0*cos(theta) + yy0*sin(theta);
						fenmu = rr2 + xx02 + yy02 + 2 * r*temp1;
						L1 = sqrt(fenmu + zz1*zz1);
						L2 = sqrt(fenmu + zz2*zz2);
						upt = -r*sin(theta);
						Vyz += (1.0 / L1 - 1.0 / L2)*upt*dtheta;
						//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
					}
					//--------------------------------------------------------------------
					Vyz = Vyz*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
					grav[j][k] += Vyz;			//
				}
			}
		}
	}
		break;
	case FORWARD_Vx:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
			double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos2[1] != forwardinfo.model.cylinder_vec[i].Pos1[1] && forwardinfo.model.cylinder_vec[i].Pos2[0] != forwardinfo.model.cylinder_vec[i].Pos1[0])
			{
				CString strinfo; strinfo.Format(_T("第%d个模型不是直立圆柱体"), i);
				MessageBox(NULL, strinfo, _T("有限长直立圆柱体正演-错误提示"), MB_OK);
				return 0;
			}
			double x0 = forwardinfo.model.cylinder_vec[i].Pos2[1];//带入公式中计算的x，y其实并不是surfer绘图中的x，y而是值北向坐标和东向坐标
			double y0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
			//==================================================上面是平面旋转的代码=====================================================================
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double north_x, east_y, x, y;
			double rr2, xx0, yy0, zz1, zz2, xx02, yy02, L1, L2, vpt, temp1, fenmu;
			double Vx;
			double dtheta = 2 * PI / theta_num;
			double theta;
			//north_x = 140;//东向切面
			//east_y = 90;//北向切面
			for (int j = 0; j < number_y; j++)
			{
				north_x = ymin + dy*j;
				//z = ymin + dy*j;//计算沿东向的切片
				for (int k = 0; k < number_x; k++)
				{
					east_y = xmin + k*dx;
					//z = xmin +dx*k;//沿北向的切片
					//将计算点变换到圆柱体坐标系中
					x = north_x; y = east_y;
					xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
					rr2 = r*r;
					xx02 = xx0*xx0; yy02 = yy0*yy0;
					Vx = 0;
					for (int kk = 0; kk < theta_num; kk++)//计算积分
					{
						theta = kk*dtheta;

						temp1 = xx0*cos(theta) + yy0*sin(theta);
						fenmu = rr2 + xx02 + yy02 + 2 * r*temp1+1E-10;
						L1 = sqrt(fenmu + zz1*zz1);
						L2 = sqrt(fenmu + zz2*zz2);
						vpt = r*cos(theta);
						Vx += 0.5*(Sign(zz2)*log((L2 - fabs(zz2)) / (L2 + fabs(zz2))) - Sign(zz1)*log((L1 - fabs(zz1)) / (L1 + fabs(zz1))))*vpt*dtheta;
						//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
					}
					//--------------------------------------------------------------------
					Vx = Vx*G*forwardinfo.model.cylinder_vec[i].Density*1E8;//E
					grav[j][k] += Vx;			//
				}
			}
		}
	}
		break;
	case FORWARD_Vxx:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
			double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos2[1] != forwardinfo.model.cylinder_vec[i].Pos1[1] && forwardinfo.model.cylinder_vec[i].Pos2[0] != forwardinfo.model.cylinder_vec[i].Pos1[0])
			{
				CString strinfo; strinfo.Format(_T("第%d个模型不是直立圆柱体"), i);
				MessageBox(NULL, strinfo, _T("有限长直立圆柱体正演-错误提示"), MB_OK);
				return 0;
			}
			double x0 = forwardinfo.model.cylinder_vec[i].Pos2[1];//带入公式中计算的x，y其实并不是surfer绘图中的x，y而是值北向坐标和东向坐标
			double y0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
			//==================================================上面是平面旋转的代码=====================================================================
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double S = PI*r*r;
			double H = (z1 - z2);
			double north_x, east_y, x, y;
			double rr2, xx0, yy0, zz1, zz2, xx02, yy02, temp1, fenzi, fenmu, temp21, temp22, utvpt;
			double Vxx;
			double dtheta = 2 * PI / theta_num;
			double theta;
			//north_x = 140;//东向切面
			//east_y = 90;//北向切面
			for (int j = 0; j < number_y; j++)
			{
				north_x = ymin + dy*j;
				//z = ymin + dy*j;//计算沿东向的切片
				for (int k = 0; k < number_x; k++)
				{
					east_y = xmin + k*dx;
					//z = xmin +dx*k;//沿北向的切片
					//将计算点变换到圆柱体坐标系中
					x = north_x; y = east_y;
					xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
					rr2 = r*r;
					xx02 = xx0*xx0; yy02 = yy0*yy0;
					Vxx = 0;
					for (int kk = 0; kk < theta_num; kk++)//计算积分
					{
						theta = kk*dtheta;

						temp1 = xx0*cos(theta) + yy0*sin(theta);
						fenzi = rr2 + r*temp1 + 1E-10;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
						fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + 1E-10;
						/*if (fenmu == 0 && fenzi == 0)
						{
						printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

						}*/
						utvpt = (r*cos(theta) + xx0)*r*cos(theta);
						temp21 = sqrt(fenmu + zz1*zz1);
						temp22 = sqrt(fenmu + zz2*zz2);
						Vxx += (utvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
						//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
					}
					//--------------------------------------------------------------------
					Vxx = Vxx*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
					grav[j][k] += Vxx;			//
				}
			}
		}
	}
		break;
	case FORWARD_Vxy:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
			double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos2[1] != forwardinfo.model.cylinder_vec[i].Pos1[1] && forwardinfo.model.cylinder_vec[i].Pos2[0] != forwardinfo.model.cylinder_vec[i].Pos1[0])
			{
				CString strinfo; strinfo.Format(_T("第%d个模型不是直立圆柱体"), i);
				MessageBox(NULL, strinfo, _T("有限长直立圆柱体正演-错误提示"), MB_OK);
				return 0;
			}
			double x0 = forwardinfo.model.cylinder_vec[i].Pos2[1];//带入公式中计算的x，y其实并不是surfer绘图中的x，y而是值北向坐标和东向坐标
			double y0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
			//==================================================上面是平面旋转的代码=====================================================================
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double S = PI*r*r;
			double H = (z1 - z2);
			double north_x, east_y, x, y;
			double rr2, xx0, yy0, zz1, zz2, xx02, yy02, temp1, fenzi, fenmu, temp21, temp22, vtvpt;
			double Vxy;
			double dtheta = 2 * PI / theta_num;
			double theta;
			//north_x = 140;//东向切面
			//east_y = 90;//北向切面
			for (int j = 0; j < number_y; j++)
			{
				north_x = ymin + dy*j;
				//z = ymin + dy*j;//计算沿东向的切片
				for (int k = 0; k < number_x; k++)
				{
					east_y = xmin + k*dx;
					//z = xmin +dx*k;//沿北向的切片
					//将计算点变换到圆柱体坐标系中
					x = north_x; y = east_y;
					xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
					rr2 = r*r;
					xx02 = xx0*xx0; yy02 = yy0*yy0;
					Vxy = 0;
					for (int kk = 0; kk < theta_num; kk++)//计算积分
					{
						theta = kk*dtheta;

						temp1 = xx0*cos(theta) + yy0*sin(theta);
						fenzi = rr2 + r*temp1 + 1E-10;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
						fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + 1E-10;
						/*if (fenmu == 0 && fenzi == 0)
						{
						printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

						}*/
						vtvpt = (r*sin(theta) + yy0)*r*cos(theta);
						temp21 = sqrt(fenmu + zz1*zz1);
						temp22 = sqrt(fenmu + zz2*zz2);
						Vxy += (vtvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
						//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
					}
					//--------------------------------------------------------------------
					Vxy = Vxy*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
					grav[j][k] += Vxy;			//
				}
			}
		}
	}
		break;
	case FORWARD_Vyy:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
			double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos2[1] != forwardinfo.model.cylinder_vec[i].Pos1[1] && forwardinfo.model.cylinder_vec[i].Pos2[0] != forwardinfo.model.cylinder_vec[i].Pos1[0])
			{
				CString strinfo; strinfo.Format(_T("第%d个模型不是直立圆柱体"), i);
				MessageBox(NULL, strinfo, _T("有限长直立圆柱体正演-错误提示"), MB_OK);
				return 0;
			}
			double x0 = forwardinfo.model.cylinder_vec[i].Pos2[1];//带入公式中计算的x，y其实并不是surfer绘图中的x，y而是值北向坐标和东向坐标
			double y0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
			//==================================================上面是平面旋转的代码=====================================================================
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double S = PI*r*r;
			double H = (z1 - z2);
			double north_x, east_y, x, y;
			double rr2, xx0, yy0, zz1, zz2, xx02, yy02, temp1, fenzi, fenmu, temp21, temp22, vtupt;
			double Vyy;
			double dtheta = 2 * PI / theta_num;
			double theta;
			//north_x = 140;//东向切面
			//east_y = 90;//北向切面
			for (int j = 0; j < number_y; j++)
			{
				north_x = ymin + dy*j;
				//z = ymin + dy*j;//计算沿东向的切片
				for (int k = 0; k < number_x; k++)
				{
					east_y = xmin + k*dx;
					//z = xmin +dx*k;//沿北向的切片
					//将计算点变换到圆柱体坐标系中
					x = north_x; y = east_y;
					xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
					rr2 = r*r;
					xx02 = xx0*xx0; yy02 = yy0*yy0;
					Vyy = 0;
					for (int kk = 0; kk < theta_num; kk++)//计算积分
					{
						theta = kk*dtheta;

						temp1 = xx0*cos(theta) + yy0*sin(theta);
						fenzi = rr2 + r*temp1 + 1E-10;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
						fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + 1E-10;
						/*if (fenmu == 0 && fenzi == 0)
						{
						printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

						}*/
						vtupt = -(r*sin(theta) + yy0)*r*sin(theta);
						temp21 = sqrt(fenmu + zz1*zz1);
						temp22 = sqrt(fenmu + zz2*zz2);
						Vyy += (vtupt / fenmu*(zz2 / temp22-zz1 / temp21 ))*dtheta;
						//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
					}
					//--------------------------------------------------------------------
					Vyy = Vyy*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
					grav[j][k] += Vyy;			//
				}
			}
		}
	}
		break;
	case FORWARD_Vy:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
			double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos2[1] != forwardinfo.model.cylinder_vec[i].Pos1[1] && forwardinfo.model.cylinder_vec[i].Pos2[0] != forwardinfo.model.cylinder_vec[i].Pos1[0])
			{
				CString strinfo; strinfo.Format(_T("第%d个模型不是直立圆柱体"), i);
				MessageBox(NULL, strinfo, _T("有限长直立圆柱体正演-错误提示"), MB_OK);
				return 0;
			}
			double x0 = forwardinfo.model.cylinder_vec[i].Pos2[1];//带入公式中计算的x，y其实并不是surfer绘图中的x，y而是值北向坐标和东向坐标
			double y0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
			//==================================================上面是平面旋转的代码=====================================================================
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double S = PI*r*r;
			double H = (z1 - z2);
			double north_x, east_y, x, y;
			double rr2, xx0, yy0, zz1, zz2, xx02, yy02, L1, L2, upt, temp1, fenmu;
			double Vy;
			double dtheta = 2 * PI / theta_num;
			double theta;
			//north_x = 140;//东向切面
			//east_y = 90;//北向切面
			for (int j = 0; j < number_y; j++)
			{
				north_x = ymin + dy*j;
				//z = ymin + dy*j;//计算沿东向的切片
				for (int k = 0; k < number_x; k++)
				{
					east_y = xmin + k*dx;
					//z = xmin +dx*k;//沿北向的切片
					//将计算点变换到圆柱体坐标系中
					x = north_x; y = east_y;
					xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
					rr2 = r*r;
					xx02 = xx0*xx0; yy02 = yy0*yy0;
					Vy = 0;
					for (int kk = 0; kk < theta_num; kk++)//计算积分
					{
						theta = kk*dtheta;

						temp1 = xx0*cos(theta) + yy0*sin(theta);
						fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + 1E-10;
						L1 = sqrt(fenmu + zz1*zz1);
						L2 = sqrt(fenmu + zz2*zz2);
						upt = -r*sin(theta);
						Vy += 0.5*(Sign(zz1)*log((L1 - fabs(zz1)) / (L1 + fabs(zz1))) - Sign(zz2)*log((L2 - fabs(zz2)) / (L2 + fabs(zz2))))*upt*dtheta;
						//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
					}
					//--------------------------------------------------------------------
					Vy = Vy*G*forwardinfo.model.cylinder_vec[i].Density*1E8;//E
					grav[j][k] += Vy;			//
				}
			}
		}
	}
		break;
	default:
		MessageBox(NULL, _T("这是重力位及其导数正演函数，请输入正确的type"), _T("有限长直立圆柱体正演-错误提示"), MB_OK);
		return 0;

	}
	GetGrdMinMax(grav, forwardinfo.model.grddatainfo);
	return 1;
}

int _3DFiniteVercicalCylinder_mag(double** mag, RegularGeometry3DForward& forwardinfo, int type, int theta_num)
{
	//首先将mag赋值为0，后面累加
	int number_x = forwardinfo.model.grddatainfo.AutoGetNumber_x();
	int number_y = forwardinfo.model.grddatainfo.AutoGetNumber_y();
	Assign_Array2(mag, number_y, number_x, 0);
	//坐标范围
	double xmin = forwardinfo.model.grddatainfo.m_AxisBounds[0];
	double xmax = forwardinfo.model.grddatainfo.m_AxisBounds[1];
	double ymin = forwardinfo.model.grddatainfo.m_AxisBounds[2];
	double ymax = forwardinfo.model.grddatainfo.m_AxisBounds[3];
	double dx = forwardinfo.model.grddatainfo.m_Dx;
	double dy = forwardinfo.model.grddatainfo.m_Dy;
	switch (type)
	{
	case FORWARD_Za:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
			double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos2[1] != forwardinfo.model.cylinder_vec[i].Pos1[1] && forwardinfo.model.cylinder_vec[i].Pos2[0] != forwardinfo.model.cylinder_vec[i].Pos1[0])
			{
				CString strinfo; strinfo.Format(_T("第%d个模型不是直立圆柱体"), i);
				MessageBox(NULL, strinfo, _T("有限长直立圆柱体正演-错误提示"), MB_OK);
				return 0;
			}
			double x0 = forwardinfo.model.cylinder_vec[i].Pos2[1];//带入公式中计算的x，y其实并不是surfer绘图中的x，y而是值北向坐标和东向坐标
			double y0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
			//==================================================上面是平面旋转的代码=====================================================================
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
			double EarthMag = forwardinfo.EarthMag;
			double M = CHL*EarthMag;											//与系数消去了U0
			double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
			double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
			double Mx = M*cos(angle_I)*sin(angle_D);
			double My = M*cos(angle_I)*cos(angle_D);
			double Mz = M*sin(angle_I);
			double H = (z1 - z2);
			double north_x, east_y, x, y;
			double rr2, xx0, yy0, zz1, zz2, xx02, yy02, temp1, fenzi, fenmu, temp21, temp22, vpt, upt;
			double Vzz,Vxz,Vyz;
			double dtheta = 2 * PI / theta_num;
			double theta;
			//north_x = 140;//东向切面
			//east_y = 90;//北向切面
			for (int j = 0; j < number_y; j++)
			{
				north_x = ymin + dy*j;
				//z = ymin + dy*j;//计算沿东向的切片
				for (int k = 0; k < number_x; k++)
				{
					east_y = xmin + k*dx;
					//z = xmin +dx*k;//沿北向的切片
					//将计算点变换到圆柱体坐标系中
					x = north_x; y = east_y;
					xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
					rr2 = r*r;
					xx02 = xx0*xx0; yy02 = yy0*yy0;
					Vzz = 0; Vxz = 0; Vyz = 0;
					for (int kk = 0; kk < theta_num; kk++)//计算积分
					{
						theta = kk*dtheta;

						temp1 = xx0*cos(theta) + yy0*sin(theta);
						fenzi = rr2 + r*temp1 + 1E-10;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
						fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + 1E-10;
						/*if (fenmu == 0 && fenzi == 0)
						{
						printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

						}*/
						temp21 = sqrt(fenmu + zz1*zz1);
						temp22 = sqrt(fenmu + zz2*zz2);
						Vzz += fenzi / fenmu*(zz2 / temp22 - zz1 / temp21)*dtheta;
						vpt = r*cos(theta);
						Vxz += (1.0 / temp22 - 1.0 / temp21)*vpt*dtheta;
						upt = -r*sin(theta);
						Vyz += (1.0 / temp21 - 1.0 / temp22)*upt*dtheta;
					}
					
					//--------------------------------------------------------------------
					mag[j][k] += 1.0 / 4.0 / PI*(Mx*Vxz + My*Vyz + Mz*Vzz);
				}
			}
		}
	}
	break;
	case FORWARD_Hax:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
			double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos2[1] != forwardinfo.model.cylinder_vec[i].Pos1[1] && forwardinfo.model.cylinder_vec[i].Pos2[0] != forwardinfo.model.cylinder_vec[i].Pos1[0])
			{
				CString strinfo; strinfo.Format(_T("第%d个模型不是直立圆柱体"), i);
				MessageBox(NULL, strinfo, _T("有限长直立圆柱体正演-错误提示"), MB_OK);
				return 0;
			}
			double x0 = forwardinfo.model.cylinder_vec[i].Pos2[1];//带入公式中计算的x，y其实并不是surfer绘图中的x，y而是值北向坐标和东向坐标
			double y0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
			//==================================================上面是平面旋转的代码=====================================================================
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
			double EarthMag = forwardinfo.EarthMag;
			double M = CHL*EarthMag;											//与系数消去了U0
			double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
			double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
			double Mx = M*cos(angle_I)*sin(angle_D);
			double My = M*cos(angle_I)*cos(angle_D);
			double Mz = M*sin(angle_I);
			double H = (z1 - z2);
			double north_x, east_y, x, y;
			double rr2, xx0, yy0, zz1, zz2, xx02, yy02, temp1, fenzi, fenmu, temp21, temp22, vpt, utvpt, vtvpt;
			double Vxx, Vxz, Vxy;
			double dtheta = 2 * PI / theta_num;
			double theta;
			//north_x = 140;//东向切面
			//east_y = 90;//北向切面
			for (int j = 0; j < number_y; j++)
			{
				north_x = ymin + dy*j;
				//z = ymin + dy*j;//计算沿东向的切片
				for (int k = 0; k < number_x; k++)
				{
					east_y = xmin + k*dx;
					//z = xmin +dx*k;//沿北向的切片
					//将计算点变换到圆柱体坐标系中
					x = north_x; y = east_y;
					xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
					rr2 = r*r;
					xx02 = xx0*xx0; yy02 = yy0*yy0;
					Vxx = 0; Vxy = 0; Vxz = 0;
					for (int kk = 0; kk < theta_num; kk++)//计算积分
					{
						theta = kk*dtheta;

						temp1 = xx0*cos(theta) + yy0*sin(theta);
						fenzi = rr2 + r*temp1 + 1E-10;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
						fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + 1E-10;
						/*if (fenmu == 0 && fenzi == 0)
						{
						printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

						}*/
						temp21 = sqrt(fenmu + zz1*zz1);
						temp22 = sqrt(fenmu + zz2*zz2);
						utvpt = (r*cos(theta) + xx0)*r*cos(theta);
						Vxx += (utvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
						vpt = r*cos(theta);
						Vxz += (1.0 / temp22 - 1.0 / temp21)*vpt*dtheta;
						vtvpt = (r*sin(theta) + yy0)*r*cos(theta);
						Vxy += (vtvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
					}
					//--------------------------------------------------------------------
					mag[j][k] += 1.0 / 4.0 / PI*(Mx*Vxx + My*Vxy + Mz*Vxz);
				}
			}
		}
	}
	break;
	case FORWARD_Hay:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			double z1 = forwardinfo.model.cylinder_vec[i].Pos1[2];
			double z2 = forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos2[1] != forwardinfo.model.cylinder_vec[i].Pos1[1] && forwardinfo.model.cylinder_vec[i].Pos2[0] != forwardinfo.model.cylinder_vec[i].Pos1[0])
			{
				CString strinfo; strinfo.Format(_T("第%d个模型不是直立圆柱体"), i);
				MessageBox(NULL, strinfo, _T("有限长直立圆柱体正演-错误提示"), MB_OK);
				return 0;
			}
			double x0 = forwardinfo.model.cylinder_vec[i].Pos2[1];//带入公式中计算的x，y其实并不是surfer绘图中的x，y而是值北向坐标和东向坐标
			double y0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
			//==================================================上面是平面旋转的代码=====================================================================
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
			double EarthMag = forwardinfo.EarthMag;
			double M = CHL*EarthMag;											//与系数消去了U0
			double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
			double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
			double Mx = M*cos(angle_I)*sin(angle_D);
			double My = M*cos(angle_I)*cos(angle_D);
			double Mz = M*sin(angle_I);
			double H = (z1 - z2);
			double north_x, east_y, x, y;
			double rr2, xx0, yy0, zz1, zz2, xx02, yy02, temp1, fenzi, fenmu, temp21, temp22, utvpt, vtvpt, vtupt, upt;
			double Vyy, Vxy, Vyz;
			double dtheta = 2 * PI / theta_num;
			double theta;
			//north_x = 140;//东向切面
			//east_y = 90;//北向切面
			for (int j = 0; j < number_y; j++)
			{
				north_x = ymin + dy*j;
				//z = ymin + dy*j;//计算沿东向的切片
				for (int k = 0; k < number_x; k++)
				{
					east_y = xmin + k*dx;
					//z = xmin +dx*k;//沿北向的切片
					//将计算点变换到圆柱体坐标系中
					x = north_x; y = east_y;
					xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
					rr2 = r*r;
					xx02 = xx0*xx0; yy02 = yy0*yy0;
					Vyy = 0; Vxy = 0; Vyz = 0;
					for (int kk = 0; kk < theta_num; kk++)//计算积分
					{
						theta = kk*dtheta;

						temp1 = xx0*cos(theta) + yy0*sin(theta);
						fenzi = rr2 + r*temp1 + 1E-10;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
						fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + 1E-10;
						/*if (fenmu == 0 && fenzi == 0)
						{
						printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

						}*/
						temp21 = sqrt(fenmu + zz1*zz1);
						temp22 = sqrt(fenmu + zz2*zz2);
						utvpt = (r*cos(theta) + xx0)*r*cos(theta);
						vtupt = -(r*sin(theta) + yy0)*r*sin(theta);
						Vyy += (vtupt / fenmu*(zz2 / temp22 - zz1 / temp21))*dtheta;
						vtvpt = (r*sin(theta) + yy0)*r*cos(theta);
						Vxy += (vtvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
						upt = -r*sin(theta);
						Vyz += (1.0 / temp21 - 1.0 / temp22)*upt*dtheta;
					}
					//--------------------------------------------------------------------
					mag[j][k] += 1.0 / 4.0 / PI*(Mx*Vyy + My*Vxy + Mz*Vyz);
				}
			}
		}
	}
		break;
	default:
		MessageBox(NULL, _T("这是磁异常分量正演函数，请输入正确的type"), _T("有限长直立圆柱体正演-错误提示"), MB_OK);
		return 0;

	}
	GetGrdMinMax(mag, forwardinfo.model.grddatainfo);
	return 1;
}

int _3DFiniteVercicalCylinder_magT(double** mag, RegularGeometry3DForward& forwardinfo, int type, int theta_num)
{
	//首先将grav赋值为0，后面累加
	int number_x = forwardinfo.model.grddatainfo.GetNumber_x();
	int number_y = forwardinfo.model.grddatainfo.GetNumber_y();
	double** Hax = CreateArray2(number_y, number_x);
	double** Hay = CreateArray2(number_y, number_x);
	double** Za = CreateArray2(number_y, number_x);
	Assign_Array2(Hax, number_y, number_x, 0);
	Assign_Array2(Hay, number_y, number_x, 0);
	Assign_Array2(Za, number_y, number_x, 0);
	Assign_Array2(mag, number_y, number_x, 0);
	double I = forwardinfo.EarthAngle_I / 180.0*PI;
	double D = forwardinfo.EarthAngle_D / 180.0*PI;

	//计算Hax
	_3DFiniteVercicalCylinder_mag(Hax, forwardinfo, FORWARD_Hax,theta_num);
	//计算Hay
	_3DFiniteVercicalCylinder_mag(Hay, forwardinfo, FORWARD_Hay, theta_num);
	//计算Za
	_3DFiniteVercicalCylinder_mag(Za, forwardinfo, FORWARD_Za, theta_num);
	switch (type)
	{
	case FORWARD_Ta:
	{
		//计算Ta
		for (int i = 0; i < number_y; i++)
		{
			for (int j = 0; j < number_x; j++)
			{
				mag[i][j] = Hax[i][j] * cos(I)*cos(D) + Hay[i][j] * cos(I)*sin(D) + Za[i][j] * sin(I);
			}
		}
	}
		break;
	case FORWARD_Module:
	{
		//计算模量
		for (int i = 0; i < number_y; i++)
		{
			for (int j = 0; j < number_x; j++)
			{
				mag[i][j] = sqrt(Hax[i][j] * Hax[i][j] + Hay[i][j] * Hay[i][j] + Za[i][j] * Za[i][j]);
			}
		}
	}
		break;
	default:
		MessageBox(NULL, _T("这是Ta和模量正演，输入正确的类型参数"), _T("出错提示"), MB_OK);
		return 0;
	}

	GetGrdMinMax(mag, forwardinfo.model.grddatainfo);
	//销毁二维数组
	DeleteArray2(Hax, number_y, number_x);
	DeleteArray2(Hay, number_y, number_x);
	DeleteArray2(Za, number_y, number_x);

	return 0;
}

int _3DFiniteCylinder_Grav(double** grav, RegularGeometry3DForward& forwardinfo0, int type, int theta_num)
{
	RegularGeometry3DForward forwardinfo = forwardinfo0;
	//首先将mag赋值为0，后面累加
	int number_x = forwardinfo.model.grddatainfo.AutoGetNumber_x();
	int number_y = forwardinfo.model.grddatainfo.AutoGetNumber_y();
	Assign_Array2(grav, number_y, number_x, 0);
	//坐标范围
	double xmin = forwardinfo.model.grddatainfo.m_AxisBounds[0];
	double xmax = forwardinfo.model.grddatainfo.m_AxisBounds[1];
	double ymin = forwardinfo.model.grddatainfo.m_AxisBounds[2];
	double ymax = forwardinfo.model.grddatainfo.m_AxisBounds[3];
	double dx = forwardinfo.model.grddatainfo.m_Dx;
	double dy = forwardinfo.model.grddatainfo.m_Dy;
	switch (type)
	{
	case FORWARD_V:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			forwardinfo.model.cylinder_vec[i].Pos1[2] = -forwardinfo.model.cylinder_vec[i].Pos1[2];
			forwardinfo.model.cylinder_vec[i].Pos2[2] - forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
			{
				for (int ii = 0; ii < 3; ii++)
				{
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
					forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
					forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
				}
			}
			//0. 与坐标系无关的量
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double z = -forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
			double S = PI*r*r;
			double EarthMag = forwardinfo.EarthMag;
			double M = S*CHL*EarthMag;											//与系数消去了U0
			//1. 根据线模型的两个端点计算转换矩阵
			double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
			forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
			temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
			forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
			double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = (forwardinfo.model.cylinder_vec[i].Pos1[2]);
			double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = (forwardinfo.model.cylinder_vec[i].Pos2[2]);

			if (x1 == x2 && y1 == y2)		//直立模型
			{
				printf("直立圆柱体\n");
				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI; 
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double rr2, temp1, fenzi, fenmu, temp21, temp22;
				double V;
				double dtheta = 2 * PI / theta_num;
				double theta;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin + dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						rr2 = r*r;
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						V = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + 1E-10;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + 1E-10;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							//V = V + fenzi/fenmu*(-zz1 / 2.0 * temp21 + fenmu / 2.0 * log(-zz1 + temp21) + zz2 / 2 * temp22 - fenmu / 2 * log(-zz2 + temp22) + zz1 / 2 * abs(zz1) - zz2 / 2 * abs(zz2))*dtheta;
							V = V + fenzi / fenmu / 2.0*(temp22*zz2 - temp21*zz1 + fenmu*log((-zz1 + temp21) / (-zz2 + temp22)) + zz1*fabs(zz1) - zz2*fabs(zz2))*dtheta;
						}
						//--------------------------------------------------------------------
						V = V*G*forwardinfo.model.cylinder_vec[i].Density*1E8;
						grav[j][k] += V;			//
					}
				}
			}
			else			//非直立模型
			{
				printf("倾斜圆柱体\n");
				//1. 圆柱体的计算倾角和偏角
				double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
				double L, L_H;
				L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
				L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
				direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
				direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
				double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
				cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
				//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
				dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
				//printf("%lf\n",dot_horizontal);
				D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
				if (dot_horizontal < 0)
				{
					if (D < 0)	//第三象限
					{
						D = -(D + PI);
					}
					else		//第四象限
					{
						D = PI - D;
					}
				}
				//printf("%lf\n", D*180/PI);
				I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

				//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
				double TransMat[9];
				double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
				TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
				TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
				TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

				//3. 根据变换矩阵计算变换后的点
				double newPos1[3], newPos2[3];
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);

				//4. 将观测坐标转换到模型坐标系进行计算
				double x0 = newPos2[0];
				double y0 = newPos2[1];
				z1 = newPos1[2], z2 = newPos2[2];
				//printf("%lf %lf \n%lf \n%lf \n",x0,y0,z1,z2);
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
				angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
				if (angle_I < 0)				//当磁化方向为向上时
				{
					angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
				}
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
				//********************
				//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
				//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double rr2, temp1, fenzi, fenmu, temp21, temp22;
				double V;
				double dtheta = 2 * PI / theta_num;
				double theta;
				double SurveyPoint[3], TransPoint[3];
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//将计算点变换到圆柱体坐标系中
						SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
						Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
						x = TransPoint[0]; y = TransPoint[1];
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						rr2 = r*r;
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						V = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + 1E-10;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + 1E-10;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							//V = V + fenzi/fenmu*(-zz1 / 2.0 * temp21 + fenmu / 2.0 * log(-zz1 + temp21) + zz2 / 2 * temp22 - fenmu / 2 * log(-zz2 + temp22) + zz1 / 2 * abs(zz1) - zz2 / 2 * abs(zz2))*dtheta;
							V = V + fenzi / fenmu / 2.0*(temp22*zz2 - temp21*zz1 + fenmu*log((-zz1 + temp21) / (-zz2 + temp22)) + zz1*fabs(zz1) - zz2*fabs(zz2))*dtheta;
						}
						//--------------------------------------------------------------------
						V = V*G*forwardinfo.model.cylinder_vec[i].Density*1E8;
						grav[j][k] += V;			//
					}
				}
			}
		}
	}
	break;
	case FORWARD_Vx:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			forwardinfo.model.cylinder_vec[i].Pos1[2] = -forwardinfo.model.cylinder_vec[i].Pos1[2];
			forwardinfo.model.cylinder_vec[i].Pos2[2] - forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
			{
				for (int ii = 0; ii < 3; ii++)
				{
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
					forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
					forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
				}
			}
			//0. 与坐标系无关的量
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double rr2 = r*r;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
			double S = PI*r*r;
			double EarthMag = forwardinfo.EarthMag;
			double M = S*CHL*EarthMag;											//与系数消去了U0
			//1. 根据线模型的两个端点计算转换矩阵
			double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
			forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
			temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
			forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
			double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = (forwardinfo.model.cylinder_vec[i].Pos1[2]);
			double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = (forwardinfo.model.cylinder_vec[i].Pos2[2]);

			if (x1 == x2 && y1 == y2)		//直立模型
			{
				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vx;
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, L1, L2, fenmu, vpt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vx = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							L1 = sqrt(fenmu + zz1*zz1);
							L2 = sqrt(fenmu + zz2*zz2);
							vpt = r*cos(theta);
							Vx += 0.5*(Sign(zz2)*log((L2 - fabs(zz2)) / (L2 + fabs(zz2))) - Sign(zz1)*log((L1 - fabs(zz1)) / (L1 + fabs(zz1))))*vpt*dtheta;
						}
						//--------------------------------------------------------------------
						Vx = Vx*G*forwardinfo.model.cylinder_vec[i].Density*1E8;//E
						grav[j][k] += Vx;			//
					}
				}
			}
			else			//非直立模型
			{
				//1. 圆柱体的计算倾角和偏角
				double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
				double L, L_H;
				L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
				L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
				direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
				direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
				double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
				cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
				//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
				dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
				//printf("%lf\n",dot_horizontal);
				D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
				if (dot_horizontal < 0)
				{
					if (D < 0)	//第三象限
					{
						D = -(D + PI);
					}
					else		//第四象限
					{
						D = PI - D;
					}
				}
				//printf("%lf\n", D*180/PI);
				I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

				//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
				double TransMat[9];
				double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
				TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
				TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
				TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

				//3. 根据变换矩阵计算变换后的点
				double newPos1[3], newPos2[3];
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
				//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
				double point_x[3]; point_x[0] = 1; point_x[1] = 0; point_x[2] = 0;
				double DirectCos[3];
				Mat_Multiply(TransMat, point_x, DirectCos, 3, 3);//北向在模型坐标系下的方向余弦
				//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
				//4. 将观测坐标转换到模型坐标系进行计算
				double x0 = newPos2[0];
				double y0 = newPos2[1];
				z1 = newPos1[2], z2 = newPos2[2];
				//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
				angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
				if (angle_I < 0)				//当磁化方向为向上时
				{
					angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
				}
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
				//********************
				//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
				//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vx, Vy, Vz;
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, L1, L2, fenmu, fenzi, vpt, upt;
				double SurveyPoint[3], TransPoint[3];
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//将计算点变换到圆柱体坐标系中
						SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
						Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
						x = TransPoint[0]; y = TransPoint[1];
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vx = 0; Vy = 0; Vz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							L1 = sqrt(fenmu + zz1*zz1);
							L2 = sqrt(fenmu + zz2*zz2);
							vpt = r*cos(theta); upt = -r*sin(theta);
							Vx += 0.5*(Sign(zz2)*log((L2 - fabs(zz2)) / (L2 + fabs(zz2))) - Sign(zz1)*log((L1 - fabs(zz1)) / (L1 + fabs(zz1))))*vpt*dtheta;
							Vy += 0.5*(Sign(zz1)*log((L1 - fabs(zz1)) / (L1 + fabs(zz1))) - Sign(zz2)*log((L2 - fabs(zz2)) / (L2 + fabs(zz2))))*upt*dtheta;
							Vz += fenzi / fenmu*(L1 - L2 - fabs(zz1) + fabs(zz2))*dtheta;
						}
						grav[j][k] += (Vx*DirectCos[0] + Vy*DirectCos[1] + Vz*DirectCos[2])*G*forwardinfo.model.cylinder_vec[i].Density*1E8;//mGal/m
					}
				}
			}
		}
	}
	break;
	case FORWARD_Vy:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			forwardinfo.model.cylinder_vec[i].Pos1[2] = -forwardinfo.model.cylinder_vec[i].Pos1[2];
			forwardinfo.model.cylinder_vec[i].Pos2[2] - forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
			{
				for (int ii = 0; ii < 3; ii++)
				{
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
					forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
					forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
				}
			}
			//0. 与坐标系无关的量
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double rr2 = r*r;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
			double S = PI*r*r;
			double EarthMag = forwardinfo.EarthMag;
			double M = S*CHL*EarthMag;											//与系数消去了U0
			//1. 根据线模型的两个端点计算转换矩阵
			double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
			forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
			temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
			forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
			double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = (forwardinfo.model.cylinder_vec[i].Pos1[2]);
			double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = (forwardinfo.model.cylinder_vec[i].Pos2[2]);

			if (x1 == x2 && y1 == y2)		//直立模型
			{
				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vy;
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, L1, L2, fenmu, upt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vy = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							L1 = sqrt(fenmu + zz1*zz1);
							L2 = sqrt(fenmu + zz2*zz2);
							upt = -r*sin(theta);
							Vy += 0.5*(Sign(zz1)*log((L1 - fabs(zz1)) / (L1 + fabs(zz1))) - Sign(zz2)*log((L2 - fabs(zz2)) / (L2 + fabs(zz2))))*upt*dtheta;
						}
						//--------------------------------------------------------------------
						Vy = Vy*G*forwardinfo.model.cylinder_vec[i].Density*1E8;//E
						grav[j][k] += Vy;			//
					}
				}
			}
			else			//非直立模型
			{
				//1. 圆柱体的计算倾角和偏角
				double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
				double L, L_H;
				L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
				L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
				direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
				direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
				double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
				cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
				//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
				dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
				//printf("%lf\n",dot_horizontal);
				D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
				if (dot_horizontal < 0)
				{
					if (D < 0)	//第三象限
					{
						D = -(D + PI);
					}
					else		//第四象限
					{
						D = PI - D;
					}
				}
				//printf("%lf\n", D*180/PI);
				I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

				//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
				double TransMat[9];
				double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
				TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
				TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
				TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

				//3. 根据变换矩阵计算变换后的点
				double newPos1[3], newPos2[3];
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
				//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
				double point_y[3]; point_y[0] = 0; point_y[1] = 1; point_y[2] = 0;
				double DirectCos[3];
				Mat_Multiply(TransMat, point_y, DirectCos, 3, 3);//北向在模型坐标系下的方向余弦
				//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
				//4. 将观测坐标转换到模型坐标系进行计算
				double x0 = newPos2[0];
				double y0 = newPos2[1];
				z1 = newPos1[2], z2 = newPos2[2];
				//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
				angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
				if (angle_I < 0)				//当磁化方向为向上时
				{
					angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
				}
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
				//********************
				//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
				//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vx, Vy, Vz;
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, L1, L2, fenmu, fenzi, vpt, upt;
				double SurveyPoint[3], TransPoint[3];
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//将计算点变换到圆柱体坐标系中
						SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
						Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
						x = TransPoint[0]; y = TransPoint[1];
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vx = 0; Vy = 0; Vz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							L1 = sqrt(fenmu + zz1*zz1);
							L2 = sqrt(fenmu + zz2*zz2);
							vpt = r*cos(theta); upt = -r*sin(theta);
							Vx += 0.5*(Sign(zz2)*log((L2 - fabs(zz2)) / (L2 + fabs(zz2))) - Sign(zz1)*log((L1 - fabs(zz1)) / (L1 + fabs(zz1))))*vpt*dtheta;
							Vy += 0.5*(Sign(zz1)*log((L1 - fabs(zz1)) / (L1 + fabs(zz1))) - Sign(zz2)*log((L2 - fabs(zz2)) / (L2 + fabs(zz2))))*upt*dtheta;
							Vz += fenzi / fenmu*(L1 - L2 - fabs(zz1) + fabs(zz2))*dtheta;
						}
						grav[j][k] += (Vx*DirectCos[0] + Vy*DirectCos[1] + Vz*DirectCos[2])*G*forwardinfo.model.cylinder_vec[i].Density*1E8;//mGal/m
					}
				}
			}
		}
	}
	break;
	case FORWARD_Vz:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			forwardinfo.model.cylinder_vec[i].Pos1[2] = -forwardinfo.model.cylinder_vec[i].Pos1[2];
			forwardinfo.model.cylinder_vec[i].Pos2[2] - forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
			{
				for (int ii = 0; ii < 3; ii++)
				{
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
					forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
					forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
				}
			}
			//0. 与坐标系无关的量
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double rr2 = r*r;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
			double S = PI*r*r;
			double EarthMag = forwardinfo.EarthMag;
			double M = S*CHL*EarthMag;											//与系数消去了U0
			//1. 根据线模型的两个端点计算转换矩阵
			double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
			forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
			temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
			forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
			double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = (forwardinfo.model.cylinder_vec[i].Pos1[2]);
			double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = (forwardinfo.model.cylinder_vec[i].Pos2[2]);

			if (x1 == x2 && y1 == y2)		//直立模型
			{
				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vz;
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, fenmu, fenzi,temp21,temp22;;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vz += fenzi / fenmu*(temp21 - temp22 - fabs(zz1) + fabs(zz2))*dtheta;
						}
						//--------------------------------------------------------------------
						Vz = Vz*G*forwardinfo.model.cylinder_vec[i].Density*1E8;//E
						grav[j][k] += Vz;			//
					}
				}
			}
			else			//非直立模型
			{
				//1. 圆柱体的计算倾角和偏角
				double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
				double L, L_H;
				L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
				L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
				direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
				direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
				double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
				cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
				//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
				dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
				//printf("%lf\n",dot_horizontal);
				D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
				if (dot_horizontal < 0)
				{
					if (D < 0)	//第三象限
					{
						D = -(D + PI);
					}
					else		//第四象限
					{
						D = PI - D;
					}
				}
				//printf("%lf\n", D*180/PI);
				I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

				//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
				double TransMat[9];
				double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
				TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
				TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
				TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

				//3. 根据变换矩阵计算变换后的点
				double newPos1[3], newPos2[3];
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
				//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
				double point_z[3]; point_z[0] = 0; point_z[1] = 0; point_z[2] = 1;
				double DirectCos[3];
				Mat_Multiply(TransMat, point_z, DirectCos, 3, 3);//北向在模型坐标系下的方向余弦
				//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
				//4. 将观测坐标转换到模型坐标系进行计算
				double x0 = newPos2[0];
				double y0 = newPos2[1];
				z1 = newPos1[2], z2 = newPos2[2];
				//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
				angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
				if (angle_I < 0)				//当磁化方向为向上时
				{
					angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
				}
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
				//********************
				//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
				//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vx, Vy, Vz;
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, L1, L2, fenmu, fenzi, vpt, upt;
				double SurveyPoint[3], TransPoint[3];
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//将计算点变换到圆柱体坐标系中
						SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
						Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
						x = TransPoint[0]; y = TransPoint[1];
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vx = 0; Vy = 0; Vz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							L1 = sqrt(fenmu + zz1*zz1);
							L2 = sqrt(fenmu + zz2*zz2);
							vpt = r*cos(theta); upt = -r*sin(theta);
							Vx += 0.5*(Sign(zz2)*log((L2 - fabs(zz2)) / (L2 + fabs(zz2))) - Sign(zz1)*log((L1 - fabs(zz1)) / (L1 + fabs(zz1))))*vpt*dtheta;
							Vy += 0.5*(Sign(zz1)*log((L1 - fabs(zz1)) / (L1 + fabs(zz1))) - Sign(zz2)*log((L2 - fabs(zz2)) / (L2 + fabs(zz2))))*upt*dtheta;
							Vz += fenzi / fenmu*(L1 - L2 - fabs(zz1) + fabs(zz2))*dtheta;
						}
						grav[j][k] += (Vx*DirectCos[0] + Vy*DirectCos[1] + Vz*DirectCos[2])*G*forwardinfo.model.cylinder_vec[i].Density*1E8;//mGal/m
					}
				}
			}
		}
	}
	break;
	case FORWARD_Vxx:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			forwardinfo.model.cylinder_vec[i].Pos1[2] = -forwardinfo.model.cylinder_vec[i].Pos1[2];
			forwardinfo.model.cylinder_vec[i].Pos2[2] - forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
			{
				for (int ii = 0; ii < 3; ii++)
				{
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
					forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
					forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
				}
			}
			//0. 与坐标系无关的量
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double rr2 = r*r;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			//1. 根据线模型的两个端点计算转换矩阵
			double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
			forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
			temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
			forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
			double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = (forwardinfo.model.cylinder_vec[i].Pos1[2]);
			double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = (forwardinfo.model.cylinder_vec[i].Pos2[2]);

			if (x1 == x2 && y1 == y2)		//直立模型
			{
				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vxx;
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi,utvpt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vxx = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							utvpt = (r*cos(theta) + xx0)*r*cos(theta);
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vxx += (utvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}
						//--------------------------------------------------------------------
						Vxx = Vxx*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
						grav[j][k] += Vxx;			//
					}
				}
			}
			else			//非直立模型
			{
				//1. 圆柱体的计算倾角和偏角
				double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
				double L, L_H;
				L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
				L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
				direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
				direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
				double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
				cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
				//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
				dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
				//printf("%lf\n",dot_horizontal);
				D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
				if (dot_horizontal < 0)
				{
					if (D < 0)	//第三象限
					{
						D = -(D + PI);
					}
					else		//第四象限
					{
						D = PI - D;
					}
				}
				//printf("%lf\n", D*180/PI);
				I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

				//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
				double TransMat[9];
				double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
				TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
				TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
				TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

				//3. 根据变换矩阵计算变换后的点
				double newPos1[3], newPos2[3];
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
				//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
				double point_x[3]; point_x[0] = 1; point_x[1] = 0; point_x[2] = 0;
				double DirectCos[3];
				Mat_Multiply(TransMat, point_x, DirectCos, 3, 3);//北向在模型坐标系下的方向余弦
				//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
				//4. 将观测坐标转换到模型坐标系进行计算
				double x0 = newPos2[0];
				double y0 = newPos2[1];
				z1 = newPos1[2], z2 = newPos2[2];
				//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
				angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
				if (angle_I < 0)				//当磁化方向为向上时
				{
					angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
				}
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
				//********************
				//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
				//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
				double SurveyPoint[3], TransPoint[3];
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi, utvpt, vtvpt, vpt, upt, vtupt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//将计算点变换到圆柱体坐标系中
						SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
						Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
						x = TransPoint[0]; y = TransPoint[1];
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vxx = 0; Vxy = 0; Vxz = 0; Vyy = 0; Vyz = 0; Vzz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							vpt = r*cos(theta); upt = -r*sin(theta);
							utvpt = (r*cos(theta) + xx0)*vpt;
							vtvpt = (r*sin(theta) + yy0)*vpt;
							vtupt = (r*sin(theta) + yy0)*upt;
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vxx += (utvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxy += (vtvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxz += (1.0 / temp22 - 1.0 / temp21)*vpt*dtheta;
							Vyy += (vtupt / fenmu*(zz2 / temp22 - zz1 / temp21))*dtheta;
							Vyz += (1.0 / temp21 - 1.0 / temp22)*upt*dtheta;
							Vzz += fenzi / fenmu*(zz2 / temp22 - zz1 / temp21)*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}
						//--------------------------------------------------------------------
						grav[j][k] += (DirectCos[0] * (Vxx*DirectCos[0] + Vxy*DirectCos[1] + Vxz*DirectCos[2]) +
									   DirectCos[1] * (Vxy*DirectCos[0] + Vyy*DirectCos[1] + Vyz*DirectCos[2]) +
									   DirectCos[2] * (Vxz*DirectCos[0] + Vyz*DirectCos[1] + Vzz*DirectCos[2]))*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
					}
				}
			}
		}
	}
	break;
	case FORWARD_Vxy:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			forwardinfo.model.cylinder_vec[i].Pos1[2] = -forwardinfo.model.cylinder_vec[i].Pos1[2];
			forwardinfo.model.cylinder_vec[i].Pos2[2] - forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
			{
				for (int ii = 0; ii < 3; ii++)
				{
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
					forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
					forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
				}
			}
			//0. 与坐标系无关的量
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double rr2 = r*r;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			//1. 根据线模型的两个端点计算转换矩阵
			double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
			forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
			temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
			forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
			double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = (forwardinfo.model.cylinder_vec[i].Pos1[2]);
			double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = (forwardinfo.model.cylinder_vec[i].Pos2[2]);

			if (x1 == x2 && y1 == y2)		//直立模型
			{
				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vxy;
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi, vtvpt,vpt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vxy = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							vpt = r*cos(theta); 
							vtvpt = (r*sin(theta) + yy0)*vpt;
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vxy += (vtvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}
						//--------------------------------------------------------------------
						Vxy = Vxy*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
						grav[j][k] += Vxy;			//
					}
				}
			}
			else			//非直立模型
			{
				//1. 圆柱体的计算倾角和偏角
				double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
				double L, L_H;
				L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
				L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
				direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
				direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
				double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
				cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
				//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
				dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
				//printf("%lf\n",dot_horizontal);
				D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
				if (dot_horizontal < 0)
				{
					if (D < 0)	//第三象限
					{
						D = -(D + PI);
					}
					else		//第四象限
					{
						D = PI - D;
					}
				}
				//printf("%lf\n", D*180/PI);
				I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

				//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
				double TransMat[9];
				double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
				TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
				TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
				TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

				//3. 根据变换矩阵计算变换后的点
				double newPos1[3], newPos2[3];
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
				//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
				double point_x[3]; point_x[0] = 1; point_x[1] = 0; point_x[2] = 0;
				double point_y[3]; point_y[0] = 0; point_y[1] = 1; point_y[2] = 0;
				double XDirectCos[3], YDirectCos[3];
				Mat_Multiply(TransMat, point_x, XDirectCos, 3, 3);//北向在模型坐标系下的方向余弦
				Mat_Multiply(TransMat, point_y, YDirectCos, 3, 3);
				//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
				//4. 将观测坐标转换到模型坐标系进行计算
				double x0 = newPos2[0];
				double y0 = newPos2[1];
				z1 = newPos1[2], z2 = newPos2[2];
				//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
				angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
				if (angle_I < 0)				//当磁化方向为向上时
				{
					angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
				}
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
				//********************
				//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
				//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
				double SurveyPoint[3], TransPoint[3];
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi, utvpt, vtvpt, vpt, upt, vtupt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//将计算点变换到圆柱体坐标系中
						SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
						Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
						x = TransPoint[0]; y = TransPoint[1];
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vxx = 0; Vxy = 0; Vxz = 0; Vyy = 0; Vyz = 0; Vzz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							vpt = r*cos(theta); upt = -r*sin(theta);
							utvpt = (r*cos(theta) + xx0)*vpt;
							vtvpt = (r*sin(theta) + yy0)*vpt;
							vtupt = (r*sin(theta) + yy0)*upt;
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vxx += (utvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxy += (vtvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxz += (1.0 / temp22 - 1.0 / temp21)*vpt*dtheta;
							Vyy += (vtupt / fenmu*(zz2 / temp22 - zz1 / temp21))*dtheta;
							Vyz += (1.0 / temp21 - 1.0 / temp22)*upt*dtheta;
							Vzz += fenzi / fenmu*(zz2 / temp22 - zz1 / temp21)*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}
						//--------------------------------------------------------------------
						grav[j][k] += (YDirectCos[0] * (Vxx*XDirectCos[0] + Vxy*XDirectCos[1] + Vxz*XDirectCos[2]) +
							YDirectCos[1] * (Vxy*XDirectCos[0] + Vyy*XDirectCos[1] + Vyz*XDirectCos[2]) +
							YDirectCos[2] * (Vxz*XDirectCos[0] + Vyz*XDirectCos[1] + Vzz*XDirectCos[2]))*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
					}
				}
			}
		}
	}
	break;
	case FORWARD_Vxz:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			forwardinfo.model.cylinder_vec[i].Pos1[2] = -forwardinfo.model.cylinder_vec[i].Pos1[2];
			forwardinfo.model.cylinder_vec[i].Pos2[2] - forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
			{
				for (int ii = 0; ii < 3; ii++)
				{
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
					forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
					forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
				}
			}
			//0. 与坐标系无关的量
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double rr2 = r*r;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			//1. 根据线模型的两个端点计算转换矩阵
			double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
			forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
			temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
			forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
			double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = (forwardinfo.model.cylinder_vec[i].Pos1[2]);
			double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = (forwardinfo.model.cylinder_vec[i].Pos2[2]);

			if (x1 == x2 && y1 == y2)		//直立模型
			{
				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vxz;
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi, vtvpt, vpt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vxz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							vpt = r*cos(theta);
							vtvpt = (r*sin(theta) + yy0)*vpt;
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vxz += (1.0 / temp22 - 1.0 / temp21)*vpt*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}
						//--------------------------------------------------------------------
						Vxz = Vxz*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
						grav[j][k] += Vxz;			//
					}
				}
			}
			else			//非直立模型
			{
				//1. 圆柱体的计算倾角和偏角
				double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
				double L, L_H;
				L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
				L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
				direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
				direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
				double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
				cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
				//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
				dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
				//printf("%lf\n",dot_horizontal);
				D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
				if (dot_horizontal < 0)
				{
					if (D < 0)	//第三象限
					{
						D = -(D + PI);
					}
					else		//第四象限
					{
						D = PI - D;
					}
				}
				//printf("%lf\n", D*180/PI);
				I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

				//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
				double TransMat[9];
				double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
				TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
				TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
				TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

				//3. 根据变换矩阵计算变换后的点
				double newPos1[3], newPos2[3];
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
				//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
				double point_x[3]; point_x[0] = 1; point_x[1] = 0; point_x[2] = 0;
				double point_z[3]; point_z[0] = 0; point_z[1] = 0; point_z[2] = 1;
				double XDirectCos[3], ZDirectCos[3];
				Mat_Multiply(TransMat, point_x, XDirectCos, 3, 3);//北向在模型坐标系下的方向余弦
				Mat_Multiply(TransMat, point_z, ZDirectCos, 3, 3);
				//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
				//4. 将观测坐标转换到模型坐标系进行计算
				double x0 = newPos2[0];
				double y0 = newPos2[1];
				z1 = newPos1[2], z2 = newPos2[2];
				//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
				angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
				if (angle_I < 0)				//当磁化方向为向上时
				{
					angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
				}
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
				//********************
				//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
				//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
				double SurveyPoint[3], TransPoint[3];
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi, utvpt, vtvpt, vpt, upt, vtupt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//将计算点变换到圆柱体坐标系中
						SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
						Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
						x = TransPoint[0]; y = TransPoint[1];
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vxx = 0; Vxy = 0; Vxz = 0; Vyy = 0; Vyz = 0; Vzz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							vpt = r*cos(theta); upt = -r*sin(theta);
							utvpt = (r*cos(theta) + xx0)*vpt;
							vtvpt = (r*sin(theta) + yy0)*vpt;
							vtupt = (r*sin(theta) + yy0)*upt;
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vxx += (utvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxy += (vtvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxz += (1.0 / temp22 - 1.0 / temp21)*vpt*dtheta;
							Vyy += (vtupt / fenmu*(zz2 / temp22 - zz1 / temp21))*dtheta;
							Vyz += (1.0 / temp21 - 1.0 / temp22)*upt*dtheta;
							Vzz += fenzi / fenmu*(zz2 / temp22 - zz1 / temp21)*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}
						//--------------------------------------------------------------------
						grav[j][k] += (ZDirectCos[0] * (Vxx*XDirectCos[0] + Vxy*XDirectCos[1] + Vxz*XDirectCos[2]) +
							ZDirectCos[1] * (Vxy*XDirectCos[0] + Vyy*XDirectCos[1] + Vyz*XDirectCos[2]) +
							ZDirectCos[2] * (Vxz*XDirectCos[0] + Vyz*XDirectCos[1] + Vzz*XDirectCos[2]))*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
					}
				}
			}
		}
	}
	break;
	case FORWARD_Vyy:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			forwardinfo.model.cylinder_vec[i].Pos1[2] = -forwardinfo.model.cylinder_vec[i].Pos1[2];
			forwardinfo.model.cylinder_vec[i].Pos2[2] - forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
			{
				for (int ii = 0; ii < 3; ii++)
				{
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
					forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
					forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
				}
			}
			//0. 与坐标系无关的量
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double rr2 = r*r;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			//1. 根据线模型的两个端点计算转换矩阵
			double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
			forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
			temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
			forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
			double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = (forwardinfo.model.cylinder_vec[i].Pos1[2]);
			double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = (forwardinfo.model.cylinder_vec[i].Pos2[2]);

			if (x1 == x2 && y1 == y2)		//直立模型
			{
				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vyy;
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi, vtupt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vyy = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							vtupt = -(r*sin(theta) + yy0)*r*sin(theta);
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vyy += (vtupt / fenmu*(zz2 / temp22 - zz1 / temp21))*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}
						//--------------------------------------------------------------------
						Vyy = Vyy*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
						grav[j][k] += Vyy;			//
					}
				}
			}
			else			//非直立模型
			{
				//1. 圆柱体的计算倾角和偏角
				double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
				double L, L_H;
				L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
				L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
				direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
				direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
				double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
				cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
				//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
				dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
				//printf("%lf\n",dot_horizontal);
				D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
				if (dot_horizontal < 0)
				{
					if (D < 0)	//第三象限
					{
						D = -(D + PI);
					}
					else		//第四象限
					{
						D = PI - D;
					}
				}
				//printf("%lf\n", D*180/PI);
				I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

				//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
				double TransMat[9];
				double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
				TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
				TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
				TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

				//3. 根据变换矩阵计算变换后的点
				double newPos1[3], newPos2[3];
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
				//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
				double point_y[3]; point_y[0] = 0; point_y[1] = 1; point_y[2] = 0;
				double DirectCos[3];
				Mat_Multiply(TransMat, point_y, DirectCos, 3, 3);//北向在模型坐标系下的方向余弦
				//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
				//4. 将观测坐标转换到模型坐标系进行计算
				double x0 = newPos2[0];
				double y0 = newPos2[1];
				z1 = newPos1[2], z2 = newPos2[2];
				//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
				angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
				if (angle_I < 0)				//当磁化方向为向上时
				{
					angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
				}
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
				//********************
				//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
				//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
				double SurveyPoint[3], TransPoint[3];
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi, utvpt, vtvpt, vpt, upt, vtupt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//将计算点变换到圆柱体坐标系中
						SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
						Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
						x = TransPoint[0]; y = TransPoint[1];
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vxx = 0; Vxy = 0; Vxz = 0; Vyy = 0; Vyz = 0; Vzz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							vpt = r*cos(theta); upt = -r*sin(theta);
							utvpt = (r*cos(theta) + xx0)*vpt;
							vtvpt = (r*sin(theta) + yy0)*vpt;
							vtupt = (r*sin(theta) + yy0)*upt;
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vxx += (utvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxy += (vtvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxz += (1.0 / temp22 - 1.0 / temp21)*vpt*dtheta;
							Vyy += (vtupt / fenmu*(zz2 / temp22 - zz1 / temp21))*dtheta;
							Vyz += (1.0 / temp21 - 1.0 / temp22)*upt*dtheta;
							Vzz += fenzi / fenmu*(zz2 / temp22 - zz1 / temp21)*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}
						//--------------------------------------------------------------------
						grav[j][k] += (DirectCos[0] * (Vxx*DirectCos[0] + Vxy*DirectCos[1] + Vxz*DirectCos[2]) +
							DirectCos[1] * (Vxy*DirectCos[0] + Vyy*DirectCos[1] + Vyz*DirectCos[2]) +
							DirectCos[2] * (Vxz*DirectCos[0] + Vyz*DirectCos[1] + Vzz*DirectCos[2]))*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
					}
				}
			}
		}
	}
	break;
	case FORWARD_Vzz:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			forwardinfo.model.cylinder_vec[i].Pos1[2] = -forwardinfo.model.cylinder_vec[i].Pos1[2];
			forwardinfo.model.cylinder_vec[i].Pos2[2] - forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
			{
				for (int ii = 0; ii < 3; ii++)
				{
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
					forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
					forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
				}
			}
			//0. 与坐标系无关的量
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double rr2 = r*r;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			//1. 根据线模型的两个端点计算转换矩阵
			double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
			forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
			temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
			forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
			double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = (forwardinfo.model.cylinder_vec[i].Pos1[2]);
			double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = (forwardinfo.model.cylinder_vec[i].Pos2[2]);

			if (x1 == x2 && y1 == y2)		//直立模型
			{
				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vzz;
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vzz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vzz += fenzi / fenmu*(zz2 / temp22 - zz1 / temp21)*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}
						//--------------------------------------------------------------------
						Vzz = Vzz*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
						grav[j][k] += Vzz;			//
					}
				}
			}
			else			//非直立模型
			{
				//1. 圆柱体的计算倾角和偏角
				double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
				double L, L_H;
				L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
				L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
				direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
				direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
				double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
				cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
				//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
				dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
				//printf("%lf\n",dot_horizontal);
				D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
				if (dot_horizontal < 0)
				{
					if (D < 0)	//第三象限
					{
						D = -(D + PI);
					}
					else		//第四象限
					{
						D = PI - D;
					}
				}
				//printf("%lf\n", D*180/PI);
				I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

				//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
				double TransMat[9];
				double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
				TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
				TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
				TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

				//3. 根据变换矩阵计算变换后的点
				double newPos1[3], newPos2[3];
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
				//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
				double point_z[3]; point_z[0] = 0; point_z[1] = 0; point_z[2] = 1;
				double DirectCos[3];
				Mat_Multiply(TransMat, point_z, DirectCos, 3, 3);//北向在模型坐标系下的方向余弦
				//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
				//4. 将观测坐标转换到模型坐标系进行计算
				double x0 = newPos2[0];
				double y0 = newPos2[1];
				z1 = newPos1[2], z2 = newPos2[2];
				//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
				angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
				if (angle_I < 0)				//当磁化方向为向上时
				{
					angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
				}
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
				//********************
				//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
				//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
				double SurveyPoint[3], TransPoint[3];
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi, utvpt, vtvpt, vpt, upt, vtupt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//将计算点变换到圆柱体坐标系中
						SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
						Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
						x = TransPoint[0]; y = TransPoint[1];
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vxx = 0; Vxy = 0; Vxz = 0; Vyy = 0; Vyz = 0; Vzz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							vpt = r*cos(theta); upt = -r*sin(theta);
							utvpt = (r*cos(theta) + xx0)*vpt;
							vtvpt = (r*sin(theta) + yy0)*vpt;
							vtupt = (r*sin(theta) + yy0)*upt;
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vxx += (utvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxy += (vtvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxz += (1.0 / temp22 - 1.0 / temp21)*vpt*dtheta;
							Vyy += (vtupt / fenmu*(zz2 / temp22 - zz1 / temp21))*dtheta;
							Vyz += (1.0 / temp21 - 1.0 / temp22)*upt*dtheta;
							Vzz += fenzi / fenmu*(zz2 / temp22 - zz1 / temp21)*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}
						//--------------------------------------------------------------------
						grav[j][k] += (DirectCos[0] * (Vxx*DirectCos[0] + Vxy*DirectCos[1] + Vxz*DirectCos[2]) +
							DirectCos[1] * (Vxy*DirectCos[0] + Vyy*DirectCos[1] + Vyz*DirectCos[2]) +
							DirectCos[2] * (Vxz*DirectCos[0] + Vyz*DirectCos[1] + Vzz*DirectCos[2]))*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
					}
				}
			}
		}
	}
	break;
	case FORWARD_Vyz:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			forwardinfo.model.cylinder_vec[i].Pos1[2] = -forwardinfo.model.cylinder_vec[i].Pos1[2];
			forwardinfo.model.cylinder_vec[i].Pos2[2] - forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
			{
				for (int ii = 0; ii < 3; ii++)
				{
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
					forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
					forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
				}
			}
			//0. 与坐标系无关的量
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double rr2 = r*r;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			//1. 根据线模型的两个端点计算转换矩阵
			double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
			forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
			temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
			forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
			double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = (forwardinfo.model.cylinder_vec[i].Pos1[2]);
			double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = (forwardinfo.model.cylinder_vec[i].Pos2[2]);

			if (x1 == x2 && y1 == y2)		//直立模型
			{
				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vyz;
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi, upt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vyz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							upt = -r*sin(theta);
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vyz += (1.0 / temp21 - 1.0 / temp22)*upt*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}
						//--------------------------------------------------------------------
						Vyz = Vyz*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
						grav[j][k] += Vyz;			//
					}
				}
			}
			else			//非直立模型
			{
				//1. 圆柱体的计算倾角和偏角
				double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
				double L, L_H;
				L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
				L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
				direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
				direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
				double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
				cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
				//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
				dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
				//printf("%lf\n",dot_horizontal);
				D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
				if (dot_horizontal < 0)
				{
					if (D < 0)	//第三象限
					{
						D = -(D + PI);
					}
					else		//第四象限
					{
						D = PI - D;
					}
				}
				//printf("%lf\n", D*180/PI);
				I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

				//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
				double TransMat[9];
				double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
				TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
				TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
				TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

				//3. 根据变换矩阵计算变换后的点
				double newPos1[3], newPos2[3];
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
				//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
				double point_y[3]; point_y[0] = 0; point_y[1] = 1; point_y[2] = 0;
				double point_z[3]; point_z[0] = 0; point_z[1] = 0; point_z[2] = 1;
				double YDirectCos[3], ZDirectCos[3];
				Mat_Multiply(TransMat, point_y, YDirectCos, 3, 3);//北向在模型坐标系下的方向余弦
				Mat_Multiply(TransMat, point_z, ZDirectCos, 3, 3);
				//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
				//4. 将观测坐标转换到模型坐标系进行计算
				double x0 = newPos2[0];
				double y0 = newPos2[1];
				z1 = newPos1[2], z2 = newPos2[2];
				//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I;//化为弧度
				angle_I = PI / 2.0 - fabs((PI / 2.0 - angle_I) - I);		//************可能有点问题，到时候数值实验中遇到不对的地方再检查这里吧
				if (angle_I < 0)				//当磁化方向为向上时
				{
					angle_I = PI / 2.0 - (fabs(angle_I) + PI / 2.0 - I);
				}
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D - D;//模型坐标系中的磁化偏角为观测坐标系磁化偏角减去模型偏角
				//********************
				//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
				//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
				double SurveyPoint[3], TransPoint[3];
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi, utvpt, vtvpt, vpt, upt, vtupt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//将计算点变换到圆柱体坐标系中
						SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
						Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
						x = TransPoint[0]; y = TransPoint[1];
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vxx = 0; Vxy = 0; Vxz = 0; Vyy = 0; Vyz = 0; Vzz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							vpt = r*cos(theta); upt = -r*sin(theta);
							utvpt = (r*cos(theta) + xx0)*vpt;
							vtvpt = (r*sin(theta) + yy0)*vpt;
							vtupt = (r*sin(theta) + yy0)*upt;
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vxx += (utvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxy += (vtvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxz += (1.0 / temp22 - 1.0 / temp21)*vpt*dtheta;
							Vyy += (vtupt / fenmu*(zz2 / temp22 - zz1 / temp21))*dtheta;
							Vyz += (1.0 / temp21 - 1.0 / temp22)*upt*dtheta;
							Vzz += fenzi / fenmu*(zz2 / temp22 - zz1 / temp21)*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}
						//--------------------------------------------------------------------
						grav[j][k] += (ZDirectCos[0] * (Vxx*YDirectCos[0] + Vxy*YDirectCos[1] + Vxz*YDirectCos[2]) +
							ZDirectCos[1] * (Vxy*YDirectCos[0] + Vyy*YDirectCos[1] + Vyz*YDirectCos[2]) +
							ZDirectCos[2] * (Vxz*YDirectCos[0] + Vyz*YDirectCos[1] + Vzz*YDirectCos[2]))*G*forwardinfo.model.cylinder_vec[i].Density*1E12;//E
					}
				}
			}
		}
	}
		break;
	default:
		MessageBox(NULL, _T("这是重力位及其导数正演函数，请输入正确的type"), _T("有限长倾斜圆柱体型正演-错误提示"), MB_OK);
		return 0;

	}
	GetGrdMinMax(grav, forwardinfo0.model.grddatainfo);
	return 1;
}

int _3DFiniteCylinder_mag(double** mag, RegularGeometry3DForward& forwardinfo0, int type, int theta_num)
{
	RegularGeometry3DForward forwardinfo = forwardinfo0;
	//首先将mag赋值为0，后面累加
	int number_x = forwardinfo.model.grddatainfo.AutoGetNumber_x();
	int number_y = forwardinfo.model.grddatainfo.AutoGetNumber_y();
	Assign_Array2(mag, number_y, number_x, 0);
	//坐标范围
	double xmin = forwardinfo.model.grddatainfo.m_AxisBounds[0];
	double xmax = forwardinfo.model.grddatainfo.m_AxisBounds[1];
	double ymin = forwardinfo.model.grddatainfo.m_AxisBounds[2];
	double ymax = forwardinfo.model.grddatainfo.m_AxisBounds[3];
	double dx = forwardinfo.model.grddatainfo.m_Dx;
	double dy = forwardinfo.model.grddatainfo.m_Dy;
	switch (type)
	{
	case FORWARD_Hax:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			forwardinfo.model.cylinder_vec[i].Pos1[2] = -forwardinfo.model.cylinder_vec[i].Pos1[2];
			forwardinfo.model.cylinder_vec[i].Pos2[2] - forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
			{
				for (int ii = 0; ii < 3; ii++)
				{
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
					forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
					forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
				}
			}
			//0. 与坐标系无关的量
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double rr2 = r*r;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
			double EarthMag = forwardinfo.EarthMag;
			double M = CHL*EarthMag;											//与系数消去了U0
			//1. 根据线模型的两个端点计算转换矩阵
			double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
			forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
			temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
			forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
			double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = (forwardinfo.model.cylinder_vec[i].Pos1[2]);
			double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = (forwardinfo.model.cylinder_vec[i].Pos2[2]);

			if (x1 == x2 && y1 == y2)		//直立模型
			{
				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vxy, Vxx, Vxz;
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi, utvpt, vtvpt, vpt, upt, vtupt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vxx = 0; Vxy = 0; Vxz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							vpt = r*cos(theta); upt = -r*sin(theta);
							utvpt = (r*cos(theta) + xx0)*vpt;
							vtvpt = (r*sin(theta) + yy0)*vpt;
							vtupt = (r*sin(theta) + yy0)*upt;
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vxx += (utvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxy += (vtvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxz += (1.0 / temp22 - 1.0 / temp21)*vpt*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}

						mag[j][k] += 1.0 / 4.0 / PI*(Mx*Vxx + My*Vxy + Mz*Vxz);
					}
				}
			}
			else			//非直立模型
			{
				//1. 圆柱体的计算倾角和偏角
				double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
				double L, L_H;
				L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
				L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
				direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
				direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
				double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
				cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
				//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
				dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
				//printf("%lf\n",dot_horizontal);
				D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
				if (dot_horizontal < 0)
				{
					if (D < 0)	//第三象限
					{
						D = -(D + PI);
					}
					else		//第四象限
					{
						D = PI - D;
					}
				}
				//printf("%lf\n", D*180/PI);
				I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

				//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
				double TransMat[9];
				double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
				TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
				TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
				TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

				//3. 根据变换矩阵计算变换后的点
				double newPos1[3], newPos2[3];
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
				//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
				double point_x[3]; point_x[0] = 1; point_x[1] = 0; point_x[2] = 0;
				double point_y[3]; point_y[0] = 0; point_y[1] = 1; point_y[2] = 0;
				double point_z[3]; point_z[0] = 0; point_z[1] = 0; point_z[2] = 1;
				double XDirectCos[3], YDirectCos[3], ZDirectCos[3];
				Mat_Multiply(TransMat, point_x, XDirectCos, 3, 3);//北向在模型坐标系下的方向余弦
				Mat_Multiply(TransMat, point_y, YDirectCos, 3, 3);
				Mat_Multiply(TransMat, point_z, ZDirectCos, 3, 3);
				//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
				//4. 将观测坐标转换到模型坐标系进行计算
				double x0 = newPos2[0];
				double y0 = newPos2[1];
				z1 = newPos1[2], z2 = newPos2[2];
				//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				//********************
				//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
				//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
				double Haxx, Haxy, Haxz;
				double SurveyPoint[3], TransPoint[3];
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi, utvpt, vtvpt, vpt, upt, vtupt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//将计算点变换到圆柱体坐标系中
						SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
						Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
						x = TransPoint[0]; y = TransPoint[1];
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vxx = 0; Vxy = 0; Vxz = 0; Vyy = 0; Vyz = 0; Vzz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							vpt = r*cos(theta); upt = -r*sin(theta);
							utvpt = (r*cos(theta) + xx0)*vpt;
							vtvpt = (r*sin(theta) + yy0)*vpt;
							vtupt = (r*sin(theta) + yy0)*upt;
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vxx += (utvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxy += (vtvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxz += (1.0 / temp22 - 1.0 / temp21)*vpt*dtheta;
							Vyy += (vtupt / fenmu*(zz2 / temp22 - zz1 / temp21))*dtheta;
							Vyz += (1.0 / temp21 - 1.0 / temp22)*upt*dtheta;
							Vzz += fenzi / fenmu*(zz2 / temp22 - zz1 / temp21)*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}

						Haxx = (XDirectCos[0] * (Vxx*XDirectCos[0] + Vxy*XDirectCos[1] + Vxz*XDirectCos[2]) +
								XDirectCos[1] * (Vxy*XDirectCos[0] + Vyy*XDirectCos[1] + Vyz*XDirectCos[2]) +
								XDirectCos[2] * (Vxz*XDirectCos[0] + Vyz*XDirectCos[1] + Vzz*XDirectCos[2]));//E
						Haxy = (YDirectCos[0] * (Vxx*XDirectCos[0] + Vxy*XDirectCos[1] + Vxz*XDirectCos[2]) +
								YDirectCos[1] * (Vxy*XDirectCos[0] + Vyy*XDirectCos[1] + Vyz*XDirectCos[2]) +
								YDirectCos[2] * (Vxz*XDirectCos[0] + Vyz*XDirectCos[1] + Vzz*XDirectCos[2]));//E
						Haxz = (ZDirectCos[0] * (Vxx*XDirectCos[0] + Vxy*XDirectCos[1] + Vxz*XDirectCos[2]) +
								ZDirectCos[1] * (Vxy*XDirectCos[0] + Vyy*XDirectCos[1] + Vyz*XDirectCos[2]) +
								ZDirectCos[2] * (Vxz*XDirectCos[0] + Vyz*XDirectCos[1] + Vzz*XDirectCos[2]));//E

						mag[j][k] += 1.0 / 4.0 / PI*(Mx*Haxx + My*Haxy + Mz*Haxz);//nT
					}
				}
			}
		}
	}
	break;
	case FORWARD_Hay:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			forwardinfo.model.cylinder_vec[i].Pos1[2] = -forwardinfo.model.cylinder_vec[i].Pos1[2];
			forwardinfo.model.cylinder_vec[i].Pos2[2] - forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
			{
				for (int ii = 0; ii < 3; ii++)
				{
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
					forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
					forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
				}
			}
			//0. 与坐标系无关的量
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double rr2 = r*r;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
			double EarthMag = forwardinfo.EarthMag;
			double M = CHL*EarthMag;											//与系数消去了U0
			//1. 根据线模型的两个端点计算转换矩阵
			double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
			forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
			temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
			forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
			double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = (forwardinfo.model.cylinder_vec[i].Pos1[2]);
			double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = (forwardinfo.model.cylinder_vec[i].Pos2[2]);

			if (x1 == x2 && y1 == y2)		//直立模型
			{
				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vxy, Vyy, Vyz;
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi, utvpt, vtvpt, vpt, upt, vtupt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vyy = 0; Vxy = 0; Vyz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							vpt = r*cos(theta); upt = -r*sin(theta);
							utvpt = (r*cos(theta) + xx0)*vpt;
							vtvpt = (r*sin(theta) + yy0)*vpt;
							vtupt = (r*sin(theta) + yy0)*upt;
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vxy += (vtvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vyy += (vtupt / fenmu*(zz2 / temp22 - zz1 / temp21))*dtheta;
							Vyz += (1.0 / temp21 - 1.0 / temp22)*upt*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}

						mag[j][k] += 1.0 / 4.0 / PI*(Mx*Vxy + My*Vyy + Mz*Vyz);
					}
				}
			}
			else			//非直立模型
			{
				//1. 圆柱体的计算倾角和偏角
				double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
				double L, L_H;
				L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
				L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
				direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
				direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
				double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
				cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
				//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
				dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
				//printf("%lf\n",dot_horizontal);
				D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
				if (dot_horizontal < 0)
				{
					if (D < 0)	//第三象限
					{
						D = -(D + PI);
					}
					else		//第四象限
					{
						D = PI - D;
					}
				}
				//printf("%lf\n", D*180/PI);
				I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

				//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
				double TransMat[9];
				double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
				TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
				TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
				TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

				//3. 根据变换矩阵计算变换后的点
				double newPos1[3], newPos2[3];
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
				//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
				double point_x[3]; point_x[0] = 1; point_x[1] = 0; point_x[2] = 0;
				double point_y[3]; point_y[0] = 0; point_y[1] = 1; point_y[2] = 0;
				double point_z[3]; point_z[0] = 0; point_z[1] = 0; point_z[2] = 1;
				double XDirectCos[3], YDirectCos[3], ZDirectCos[3];
				Mat_Multiply(TransMat, point_x, XDirectCos, 3, 3);//北向在模型坐标系下的方向余弦
				Mat_Multiply(TransMat, point_y, YDirectCos, 3, 3);
				Mat_Multiply(TransMat, point_z, ZDirectCos, 3, 3);
				//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
				//4. 将观测坐标转换到模型坐标系进行计算
				double x0 = newPos2[0];
				double y0 = newPos2[1];
				z1 = newPos1[2], z2 = newPos2[2];
				//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				//********************
				//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
				//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
				double Hayy, Haxy, Hayz;
				double SurveyPoint[3], TransPoint[3];
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi, utvpt, vtvpt, vpt, upt, vtupt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//将计算点变换到圆柱体坐标系中
						SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
						Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
						x = TransPoint[0]; y = TransPoint[1];
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vxx = 0; Vxy = 0; Vxz = 0; Vyy = 0; Vyz = 0; Vzz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							vpt = r*cos(theta); upt = -r*sin(theta);
							utvpt = (r*cos(theta) + xx0)*vpt;
							vtvpt = (r*sin(theta) + yy0)*vpt;
							vtupt = (r*sin(theta) + yy0)*upt;
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vxx += (utvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxy += (vtvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxz += (1.0 / temp22 - 1.0 / temp21)*vpt*dtheta;
							Vyy += (vtupt / fenmu*(zz2 / temp22 - zz1 / temp21))*dtheta;
							Vyz += (1.0 / temp21 - 1.0 / temp22)*upt*dtheta;
							Vzz += fenzi / fenmu*(zz2 / temp22 - zz1 / temp21)*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}

						Hayy = (YDirectCos[0] * (Vxx*YDirectCos[0] + Vxy*YDirectCos[1] + Vxz*YDirectCos[2]) +
								YDirectCos[1] * (Vxy*YDirectCos[0] + Vyy*YDirectCos[1] + Vyz*YDirectCos[2]) +
								YDirectCos[2] * (Vxz*YDirectCos[0] + Vyz*YDirectCos[1] + Vzz*YDirectCos[2]));//E
						Haxy = (XDirectCos[0] * (Vxx*YDirectCos[0] + Vxy*YDirectCos[1] + Vxz*YDirectCos[2]) +
								XDirectCos[1] * (Vxy*YDirectCos[0] + Vyy*YDirectCos[1] + Vyz*YDirectCos[2]) +
								XDirectCos[2] * (Vxz*YDirectCos[0] + Vyz*YDirectCos[1] + Vzz*YDirectCos[2]));//E
						Hayz = (ZDirectCos[0] * (Vxx*YDirectCos[0] + Vxy*YDirectCos[1] + Vxz*YDirectCos[2]) +
								ZDirectCos[1] * (Vxy*YDirectCos[0] + Vyy*YDirectCos[1] + Vyz*YDirectCos[2]) +
								ZDirectCos[2] * (Vxz*YDirectCos[0] + Vyz*YDirectCos[1] + Vzz*YDirectCos[2]));//E

						mag[j][k] += 1.0 / 4.0 / PI*(Mx*Haxy + My*Hayy + Mz*Hayz);//nT
					}
				}
			}
		}
	}
	break;
	case FORWARD_Za:
	{
		for (int i = 0; i < (int)forwardinfo.model.cylinder_vec.size(); i++)
		{
			forwardinfo.model.cylinder_vec[i].Pos1[2] = -forwardinfo.model.cylinder_vec[i].Pos1[2];
			forwardinfo.model.cylinder_vec[i].Pos2[2] - forwardinfo.model.cylinder_vec[i].Pos2[2];
			if (forwardinfo.model.cylinder_vec[i].Pos1[2]>forwardinfo.model.cylinder_vec[i].Pos2[2])		//判断，确保第一个点是z值较小的点
			{
				for (int ii = 0; ii < 3; ii++)
				{
					double temp = forwardinfo.model.cylinder_vec[i].Pos1[ii];
					forwardinfo.model.cylinder_vec[i].Pos1[ii] = forwardinfo.model.cylinder_vec[i].Pos2[ii];
					forwardinfo.model.cylinder_vec[i].Pos2[ii] = temp;
				}
			}
			//0. 与坐标系无关的量
			double r = forwardinfo.model.cylinder_vec[i].Radius;
			double rr2 = r*r;
			double z = forwardinfo.model.grddatainfo.m_Height_data;//加上正演高度
			double CHL = forwardinfo.model.cylinder_vec[i].CiHuaLv;
			double EarthMag = forwardinfo.EarthMag;
			double M = CHL*EarthMag;											//与系数消去了U0
			//1. 根据线模型的两个端点计算转换矩阵
			double temp = forwardinfo.model.cylinder_vec[i].Pos1[0]; forwardinfo.model.cylinder_vec[i].Pos1[0] = forwardinfo.model.cylinder_vec[i].Pos1[1];
			forwardinfo.model.cylinder_vec[i].Pos1[1] = temp;
			temp = forwardinfo.model.cylinder_vec[i].Pos2[0]; forwardinfo.model.cylinder_vec[i].Pos2[0] = forwardinfo.model.cylinder_vec[i].Pos2[1];
			forwardinfo.model.cylinder_vec[i].Pos2[1] = temp;
			double x1 = forwardinfo.model.cylinder_vec[i].Pos1[0], y1 = forwardinfo.model.cylinder_vec[i].Pos1[1], z1 = (forwardinfo.model.cylinder_vec[i].Pos1[2]);
			double x2 = forwardinfo.model.cylinder_vec[i].Pos2[0], y2 = forwardinfo.model.cylinder_vec[i].Pos2[1], z2 = (forwardinfo.model.cylinder_vec[i].Pos2[2]);

			if (x1 == x2 && y1 == y2)		//直立模型
			{
				double x0 = forwardinfo.model.cylinder_vec[i].Pos2[0];
				double y0 = forwardinfo.model.cylinder_vec[i].Pos2[1];
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vxz, Vyz, Vzz;
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi,vpt, upt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					//z = ymin+dy*j;//计算沿东向的切片
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//z = xmin +dx*k;//沿北向的切片
						//将计算点变换到圆柱体坐标系中
						x = north_x; y = east_y;
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - z; zz2 = z2 - z;//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vxz = 0; Vyz = 0; Vzz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							vpt = r*cos(theta); upt = -r*sin(theta);
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vxz += (1.0 / temp22 - 1.0 / temp21)*vpt*dtheta;
							Vyz += (1.0 / temp21 - 1.0 / temp22)*upt*dtheta;
							Vzz += fenzi / fenmu*(zz2 / temp22 - zz1 / temp21)*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}

						mag[j][k] += 1.0 / 4.0 / PI*(Mx*Vxz + My*Vyz + Mz*Vzz);
					}
				}
			}
			else			//非直立模型
			{
				//1. 圆柱体的计算倾角和偏角
				double direct_vector_horizontal[3], direct_vector_x[3];//圆柱体在水平面投影的单位向量（小z指向大z）
				double L, L_H;
				L_H = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));//水平面投影长度
				L = sqrt(L_H*L_H + (z1 - z2)*(z1 - z2));	//圆柱体长度
				direct_vector_horizontal[0] = (x2 - x1) / L_H; direct_vector_horizontal[1] = (y2 - y1) / L_H; direct_vector_horizontal[2] = 0;
				direct_vector_x[0] = 1; direct_vector_x[1] = 0; direct_vector_x[2] = 0;
				double *cross_horizontal, dot_horizontal, D, I;//x与水平投影方向的叉乘及点乘，偏角,倾角
				cross_horizontal = Cross(direct_vector_x, direct_vector_horizontal);
				//printf("%lf  %lf  %lf ", cross_horizontal[0], cross_horizontal[1], cross_horizontal[2]);
				dot_horizontal = VectorDot(direct_vector_x, direct_vector_horizontal);
				//printf("%lf\n",dot_horizontal);
				D = asin(cross_horizontal[2]);	//计算圆柱体的偏角
				if (dot_horizontal < 0)
				{
					if (D < 0)	//第三象限
					{
						D = -(D + PI);
					}
					else		//第四象限
					{
						D = PI - D;
					}
				}
				//printf("%lf\n", D*180/PI);
				I = fabs(asin(L_H / L));// printf("%lf\n", I* 180 / PI);//圆柱体与z轴的夹角

				//2. 计算圆柱体坐标系和观测坐标系之间的变换矩阵
				double TransMat[9];
				double cosD = cos(D), sinD = sin(D), cosI = cos(I), sinI = sin(I);
				TransMat[0] = cosI*cosD; TransMat[1] = cosI*sinD; TransMat[2] = -sinI;
				TransMat[3] = -sinD; TransMat[4] = cosD; TransMat[5] = 0;
				TransMat[6] = sinI*cosD; TransMat[7] = sinI*sinD; TransMat[8] = cosI;

				//3. 根据变换矩阵计算变换后的点
				double newPos1[3], newPos2[3];
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos1, newPos1, 3, 3);
				Mat_Multiply(TransMat, forwardinfo.model.cylinder_vec[i].Pos2, newPos2, 3, 3);
				//计算观测坐标系中x轴（也就是北向）在模型坐标系中的方向向量，方便计算方向导数
				double point_x[3]; point_x[0] = 1; point_x[1] = 0; point_x[2] = 0;
				double point_y[3]; point_y[0] = 0; point_y[1] = 1; point_y[2] = 0;
				double point_z[3]; point_z[0] = 0; point_z[1] = 0; point_z[2] = 1;
				double XDirectCos[3], YDirectCos[3], ZDirectCos[3];
				Mat_Multiply(TransMat, point_x, XDirectCos, 3, 3);//北向在模型坐标系下的方向余弦
				Mat_Multiply(TransMat, point_y, YDirectCos, 3, 3);
				Mat_Multiply(TransMat, point_z, ZDirectCos, 3, 3);
				//printf("%lf %lf %lf\n%lf\n", NewPoint_x[0], NewPoint_x[1], NewPoint_x[2], sqrt(NewPoint_x[0] * NewPoint_x[0] + NewPoint_x[1] * NewPoint_x[1] + NewPoint_x[2] * NewPoint_x[2]));
				//4. 将观测坐标转换到模型坐标系进行计算
				double x0 = newPos2[0];
				double y0 = newPos2[1];
				z1 = newPos1[2], z2 = newPos2[2];
				//printf("%lf %lf \n%lf \n%lf \n", x0, y0, z1, z2);
				double angle_I = forwardinfo.model.cylinder_vec[i].Angle_I / 180.0*PI;//化为弧度
				double angle_D = forwardinfo.model.cylinder_vec[i].Angle_D / 180.0*PI;
				//********************
				//********************最好这里测试一下转换到模型坐标系后的磁化倾角和偏角
				//printf("变换后的磁化偏角: %lf \n变换后的磁化倾角: %lf\n", angle_D/PI*180, angle_I/PI*180);
				double Mx = M*cos(angle_I)*sin(angle_D);
				double My = M*cos(angle_I)*cos(angle_D);
				double Mz = M*sin(angle_I);
				double north_x, east_y, x, y;
				double xx0, yy0, zz1, zz2, xx02, yy02;
				double Vxx, Vxy, Vxz, Vyy, Vyz, Vzz;
				double Haxz, Hayz, Hazz;
				double SurveyPoint[3], TransPoint[3];
				double dtheta = 2 * PI / theta_num;
				double theta, temp1, temp21, temp22, fenmu, fenzi, utvpt, vtvpt, vpt, upt, vtupt;
				for (int j = 0; j < number_y; j++)
				{
					north_x = ymin + dy*j;
					for (int k = 0; k < number_x; k++)
					{
						east_y = xmin + k*dx;
						//将计算点变换到圆柱体坐标系中
						SurveyPoint[0] = north_x; SurveyPoint[1] = east_y; SurveyPoint[2] = z;
						Mat_Multiply(TransMat, SurveyPoint, TransPoint, 3, 3);
						x = TransPoint[0]; y = TransPoint[1];
						xx0 = x0 - x; yy0 = y0 - y; zz1 = z1 - TransPoint[2]; zz2 = z2 - TransPoint[2];//这里的xx0和yy0都是考虑的了南北方向与xy的关系，下面的x都是代表北向坐标；而y0表示传入的按照surfer图的坐标也表示北向的
						xx02 = xx0*xx0; yy02 = yy0*yy0;
						Vxx = 0; Vxy = 0; Vxz = 0; Vyy = 0; Vyz = 0; Vzz = 0;
						for (int kk = 0; kk < theta_num; kk++)//计算积分
						{
							theta = kk*dtheta;

							temp1 = xx0*cos(theta) + yy0*sin(theta);
							fenzi = rr2 + r*temp1 + EPS;//当观测点在x，y轴与圆的交点处时，会出现分子分母同时为0的情况，经推导这是个同届无穷小，同加一个小量就行
							fenmu = rr2 + xx02 + yy02 + 2 * r*temp1 + EPS;
							/*if (fenmu == 0 && fenzi == 0)
							{
							printf("奇异点: %lf   %lf   %lf\n",x,y,theta);

							}*/
							vpt = r*cos(theta); upt = -r*sin(theta);
							utvpt = (r*cos(theta) + xx0)*vpt;
							vtvpt = (r*sin(theta) + yy0)*vpt;
							vtupt = (r*sin(theta) + yy0)*upt;
							temp21 = sqrt(fenmu + zz1*zz1);
							temp22 = sqrt(fenmu + zz2*zz2);
							Vxx += (utvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxy += (vtvpt / fenmu*(zz1 / temp21 - zz2 / temp22))*dtheta;
							Vxz += (1.0 / temp22 - 1.0 / temp21)*vpt*dtheta;
							Vyy += (vtupt / fenmu*(zz2 / temp22 - zz1 / temp21))*dtheta;
							Vyz += (1.0 / temp21 - 1.0 / temp22)*upt*dtheta;
							Vzz += fenzi / fenmu*(zz2 / temp22 - zz1 / temp21)*dtheta;
							//printf("theta: %lf\ntemp1: %lf\nfenzi: %lf\nfenmu: %lf\ntemp21: %lf\ntemp22: %lf\nV: %lf\n", theta,temp1,fenzi,fenmu,temp21,temp22,V); return 0;
						}

						Haxz = (XDirectCos[0] * (Vxx*ZDirectCos[0] + Vxy*ZDirectCos[1] + Vxz*ZDirectCos[2]) +
							XDirectCos[1] * (Vxy*ZDirectCos[0] + Vyy*ZDirectCos[1] + Vyz*ZDirectCos[2]) +
							XDirectCos[2] * (Vxz*ZDirectCos[0] + Vyz*ZDirectCos[1] + Vzz*ZDirectCos[2]));//E
						Hayz = (YDirectCos[0] * (Vxx*ZDirectCos[0] + Vxy*ZDirectCos[1] + Vxz*ZDirectCos[2]) +
							YDirectCos[1] * (Vxy*ZDirectCos[0] + Vyy*ZDirectCos[1] + Vyz*ZDirectCos[2]) +
							YDirectCos[2] * (Vxz*ZDirectCos[0] + Vyz*ZDirectCos[1] + Vzz*ZDirectCos[2]));//E
						Haxz = (ZDirectCos[0] * (Vxx*ZDirectCos[0] + Vxy*ZDirectCos[1] + Vxz*ZDirectCos[2]) +
							ZDirectCos[1] * (Vxy*ZDirectCos[0] + Vyy*ZDirectCos[1] + Vyz*ZDirectCos[2]) +
							ZDirectCos[2] * (Vxz*ZDirectCos[0] + Vyz*ZDirectCos[1] + Vzz*ZDirectCos[2]));//E

						mag[j][k] += 1.0 / 4.0 / PI*(Mx*Haxz + My*Hayz + Mz*Hazz);//nT
					}
				}
			}
		}
	}
		break;
	default:
		MessageBox(NULL, _T("这是磁异常三分量正演函数，请输入正确的type"), _T("有限长倾斜圆柱体模型正演-错误提示"), MB_OK);
		return 0;

	}
	GetGrdMinMax(mag, forwardinfo0.model.grddatainfo);
	return 1;
}

int _3DFiniteCylinder_magT(double** mag, RegularGeometry3DForward& forwardinfo, int type, int theta_num)
{
	//首先将grav赋值为0，后面累加
	int number_x = forwardinfo.model.grddatainfo.GetNumber_x();
	int number_y = forwardinfo.model.grddatainfo.GetNumber_y();
	double** Hax = CreateArray2(number_y, number_x);
	double** Hay = CreateArray2(number_y, number_x);
	double** Za = CreateArray2(number_y, number_x);
	Assign_Array2(Hax, number_y, number_x, 0);
	Assign_Array2(Hay, number_y, number_x, 0);
	Assign_Array2(Za, number_y, number_x, 0);
	Assign_Array2(mag, number_y, number_x, 0);
	double I = forwardinfo.EarthAngle_I / 180.0*PI;
	double D = forwardinfo.EarthAngle_D / 180.0*PI;

	//计算Hax
	_3DFiniteCylinder_mag(Hax, forwardinfo, FORWARD_Hax, theta_num);
	//计算Hay
	_3DFiniteCylinder_mag(Hay, forwardinfo, FORWARD_Hay, theta_num);
	//计算Za
	_3DFiniteCylinder_mag(Za, forwardinfo, FORWARD_Za, theta_num);
	switch (type)
	{
	case FORWARD_Ta:
	{
		//计算Ta
		for (int i = 0; i < number_y; i++)
		{
			for (int j = 0; j < number_x; j++)
			{
				mag[i][j] = Hax[i][j] * cos(I)*cos(D) + Hay[i][j] * cos(I)*sin(D) + Za[i][j] * sin(I);
			}
		}
	}
		break;
	case FORWARD_Module:
	{
		//计算模量
		for (int i = 0; i < number_y; i++)
		{
			for (int j = 0; j < number_x; j++)
			{
				mag[i][j] = sqrt(Hax[i][j] * Hax[i][j] + Hay[i][j] * Hay[i][j] + Za[i][j] * Za[i][j]);
			}
		}
	}
		break;
	default:
		MessageBox(NULL, _T("这是Ta和模量正演，输入正确的类型参数"), _T("有限长倾斜圆柱体模型正演-出错提示"), MB_OK);
		return 0;
	}

	GetGrdMinMax(mag, forwardinfo.model.grddatainfo);
	//销毁二维数组
	DeleteArray2(Hax, number_y, number_x);
	DeleteArray2(Hay, number_y, number_x);
	DeleteArray2(Za, number_y, number_x);

	return 0;
}