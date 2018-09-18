#include<iostream.h>
#include<stdio.h>
#include<math.h>
#include<string>
#include<fstream>
#include<sstream>
using namespace std;

const double I=45; // 正常场磁倾角 单位：度
const double A=0;  // 正常场磁偏角 单位：度

string filename; // 原始磁异常
string filename_hax,filename_hay,filename_za;  // 用于保存计算的三分量Hax、Hay和Za
string filename_Ta,filename_R, filename_E, filename_Q, filename_L; // 用于保存五种磁异常模量
string Filename_Ta,Filename_R, Filename_E, Filename_Q, Filename_L; // 用于保存五种磁异常模量边界识别结果


int XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,Z_MIN,Z_MAX;

double open_grd(string Filename1);
double save_grd(string Filename2,,xn,yn,x_min,x_max,y_min,y_max,data[]);
int dew;
int main()
{
	// 读取文件
	[XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,Z_MIN,Z_MAX,data]=open_grd("组合直立长方体Delta_T_不带剩磁.grd");
    
	// 计算原始异常的三分量：Hax,Hay,Za
	double Hax[YN][XN];
	double Hay[YN][XN];
	double Za[YN][XN];

    Hax=save_grd(filename_hax,XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,Data[]);
	Hay=save_grd(filename_hay,XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,Data[]);
	Za=save_grd(filename_za,XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,Data[]);

    // 利用余弦扩编进行求导运算，求导函数分别为：
	// Derivative_direction_X、Derivative_direction_Y、Derivative_direction_Z
    
    
    // 三分量Hax、Hay和Za的各方向的导数分别为:Hax_x、Hax_y和Hax_z
	//                                        Hay_x、Hay_y和Hay_z
	//                                        Za_x、Za_y和Za_z


    // 模量Ta
	int i,j;
	double Ta[YN][XN];
	for (i=1;i<=YN;i++)
	{
		for (j=1;j<=XN;j++)
		{
			Ta[i][j]=sqrt(Hax[i][j]*Hax[i][j]+Hay[i][j]*Hay[i][j]+Za[i][j]*Za[i][j]);
		}
	}


    // 模量R
    double Tax[YN][XN];	
	double Tay[YN][XN];
	double Taz[YN][XN];
	double R[YN][XN];
    for (i=1;i<=YN;i++)
	{
		for (j=1;j<=XN;j++)
		{
			Tax[i][j]=(Hax[i][j]*Hax_x[i][j]+Hay[i][j]*Hay_x[i][j]+Za[i][j]*Za_x[i][j])/Ta[i][j];
		    Tay[i][j]=(Hax[i][j]*Hax_y[i][j]+Hay[i][j]*Hay_y[i][j]+Za[i][j]*Za_y[i][j])/Ta[i][j];
			Taz[i][j]=(Hax[i][j]*Hax_z[i][j]+Hay[i][j]*Hay_z[i][j]+Za[i][j]*Za_z[i][j])/Ta[i][j];
			R[i][j]=sqrt(Tax[i][j]*Tax[i][j]+Tay[i][j]*Tay[i][j]+Taz[i][j]*Taz[i][j]);
		}
	}


   // 模量E
	double Theta_Hx[YN][XN];
    double Theta_Hy[YN][XN];
	double Theta_Za[YN][XN];
	double E[YN][XN];
    for (i=1;i<=YN;i++)
	{
		for (j=1;j<=XN;j++)
		{
			Theta_Hx[i][j]=Hax_x[i][j]*Hax_x[i][j]+Hax_y[i][j]*Hax_y[i][j]+Hax_z[i][j]*Hax_z[i][j];
			Theta_Hy[i][j]=Hay_x[i][j]*Hay_x[i][j]+Hay_y[i][j]*Hay_y[i][j]+Hay_z[i][j]*Hay_z[i][j];
			Theta_Za[i][j]=Za_x[i][j]*Za_x[i][j]+Za_y[i][j]*Za_y[i][j]+Za_z[i][j]*Za_z[i][j];
            E[i][j]=sqrt((Theta_Hx[i][j]+Theta_Hy[i][j]+Theta_Za[i][j])/2);
		}
	}

 
    // 模量Q
    double Q[YN][XN];
	for (i=1;i<=YN;i++)
	{
		for (j=1;i<=XN;j++)
		{
			Q[i][j]=sqrt(fabs(2*E[i][j]*E[i][j]-R[i][j]*R[i][j]));
		}
	}

   // 模量L
	double L[YN][XN];
	for (i=1;i<=YN;i++)
	{
		for (j=1;j<=XN;j++)
		{
			L[i][j]=Q[i][j]*Q[i][j]/Ta[i][j];
		}
	}

    // 各模量的方向导数分别为：Ta_x、Ta_y、Ta_z; R_x、R_y、R_z
    //                         E_x、E_y、E_z; Q_x、Q_y、Q_z; L_x、L_y、L_z
    
	Ta_x=
    Ta_y=
    Ta_z=

    R_x=
    R_y=
    R_z=

    E_x=
    E_y=
    E_z=

    Q_x=
	Q_y=
	Q_z=

	L_x=
	L_y=
	L_z=

    double Ta_Alplitude[YN][XN];
	double R_Alplitude[YN][XN];
	double E_Alplitude[YN][XN];
	double Q_Alplitude[YN][XN];
	double L_Alplitude[YN][XN];

	double Ta_Tdr_improved[YN][XN]; // 改进型的Tilt导数
	double R_Tdr_improved[YN][XN];
	double E_Tdr_improved[YN][XN];
	double Q_Tdr_improved[YN][XN];
	double L_Tdr_improved[YN][XN];

    for (i=1;i<=YN;i++)
	{
		for (j=1;j<=XN;j++)
		{
			// 各模量的总梯度模
			Ta_Alplitude[i][j]=sqrt(Ta_x[i][j]*Ta_x[i][j]+Ta_y[i][j]*Ta_y[i][j]+Ta_z[i][j]*Ta_z[i][j]);
			R_Alplitude[i][j]=sqrt(R_x[i][j]*R_x[i][j]+R_y[i][j]*R_y[i][j]+R_z[i][j]*R_z[i][j]);
			E_Alplitude[i][j]=sqrt(E_x[i][j]*E_x[i][j]+E_y[i][j]*E_y[i][j]+E_z[i][j]*E_z[i][j]);
			Q_Alplitude[i][j]=sqrt(Q_x[i][j]*Q_x[i][j]+Q_y[i][j]*Q_y[i][j]+Q_z[i][j]*Q_z[i][j]);
			L_Alplitude[i][j]=sqrt(L_x[i][j]*L_x[i][j]+L_y[i][j]*L_y[i][j]+L_z[i][j]*L_z[i][j]);
            
			// 改进型的Tilt导数，即边界识别结果
			Ta_Tdr_improved[i][j]=atan2(Ta_z[i][j],Ta_Alplitude[i][j]);
            R_Tdr_improved[i][j]=atan2(R_z[i][j],R_Alplitude[i][j]);
			E_Tdr_improved[i][j]=atan2(E_z[i][j],E_Alplitude[i][j]);
			Q_Tdr_improved[i][j]=atan2(Q_z[i][j],Q_Alplitude[i][j]);
			L_Tdr_improved[i][j]=atan2(L_z[i][j],L_Alplitude[i][j]);

            // 进行判断
			if Ta_Tdr_improved(i,j)<0
				Ta_Tdr_improved(i,j)=0;
			else if R_Tdr_improved(i,j)<0
				R_Tdr_improved(i,j)=0;
			else if E_Tdr_improved(i,j)<0
				E_Tdr_improved(i,j)=0;
			else if Q_Tdr_improved(i,j)<0
				Q_Tdr_improved(i,j)=0;
			else if L_Tdr_improved(i,j)<0
				L_Tdr_improved(i,j)=0;

		// 保存磁异常三分量Hax、Hay、Za 网格数据
			save_grd;

		// 保存五种模量网格数据：Ta、R、E、Q、L
			save_grd;
		
		// 保存五种模量边界识别网格数据：Ta、R、E、Q、L
            save_grd;


		}
	}


	delete [] obdata;
    return 0;
}





// 用于读取grd文件
double open_grd(filename1)
{
	ifstream infile;
	FILE *fp;
	if ((fp=fopen(filename1,"rb"))==NULL)
	{
		cout<<"打不开文件！！"<<endl;
		return 0;
	}
    
	string filehead; // 存储头文件DSAA
	getline(infile,filehead);  // 读取grd数据文件头DSAA
	infile>>XN>>YN>>X_MIN>>X_MAX>>Y_MIN>>Y_MAX>>Z_MIN>>Z_MAX;
    double *obdata=new double [XN*YN];
	for (i=1;i<XN*YN;i++)
	{
		infile>>obdata[i];
	}
	infile.close();

	return 0;
}


// 用于保存grd文件
double save_grd(filename2,XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,Data[])
{
	ifstream outGrid;
	FILE *fp;
	if ((fp=fopen(filename2,"w"))==NULL)
	{
		cout<<"打不开文件！！"<<endl;
		return 0;
	}
    outGrid<<"DSAA"<<endl;
	outGrid<<XN<<"  "<<YN<<endl;
	outGrid<<X_MIN<<"  "<<X_MAX<<endl;
    outGrid<<Y_MIN<<"  "<<Y_MAX<<endl;
    
	// 查找数据中的极大值和极小值
	double Z_MIN,Z_MAX;
	Z_MIN=Z_MAX=data[0];

    for(i=1;i<XN*YN;i++)
    {
       if(obdata[i]>MAX)
           MAX=obdata[i];
	   else if(obdata[i]<MIN)
		   MIN=obdata[i];
	}

    outGrid<<MIN<<"  "<<MAX<<endl; // 输出极小值和极大值
	for(i=0;i<datanum;i++)
	{
         outGrid<<obdata[i]<<"  "; // 输出观测数据
	}
    outGrid.close();

	return 0;
}
