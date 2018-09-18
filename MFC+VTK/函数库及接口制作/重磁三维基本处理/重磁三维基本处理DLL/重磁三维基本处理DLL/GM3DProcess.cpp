

#include "GM3DProcess.h"
#include "FFTN.h"
namespace GM3DProcess
{
	int DerivativeY(CGM_GRDDATA indata, CGM_GRDDATA* processdata)
	{
		double tp = 2 * PI;
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;

		double* originaldata = indata.GetData();
		fftw_complex* out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
		FFT2d(number_y, number_x, originaldata, out_c);			//傅里叶变换
		double *bounds = indata.GetBounds();
		double dx = fabs(bounds[1] - bounds[0]) / (number_x - 1);
		double em = dx*(number_x);
		double u;
		double real, imag;
		int number_x_half = (int)(number_x / 2);
		int index_ij;
		for (int i = 0; i < number_y; i++)
		{
			for (int j = 0; j <= number_x_half; j++)
			{
				u = j / em;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0] * tp*u;
				imag = out_c[index_ij][1] * tp*u;
				out_c[index_ij][0] = -imag;
				out_c[index_ij][1] = real;
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				u = -(number_x - j) / em;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0] * tp*u;
				imag = out_c[index_ij][1] * tp*u;
				out_c[index_ij][0] = -imag;
				out_c[index_ij][1] = real;
			}
		}
		IFFT2d(number_y, number_x, out_c, processdata->GetData());
		return 1;
	}

	int DerivativeX(CGM_GRDDATA indata, CGM_GRDDATA* processdata)
	{
		double tp = 2 * PI;
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double* originaldata = indata.GetData();
		fftw_complex* out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
		FFT2d(number_y, number_x, originaldata, out_c);			//傅里叶变换
		double *bounds = indata.GetBounds();
		double dy = fabs(bounds[3] - bounds[2]) / (number_y - 1);
		double en = dy*(number_y);
		double v;
		double real, imag;
		int number_y_half = (int)(number_y / 2);
		int index_ij;
		for (int i = 0; i <= number_y_half; i++)
		{
			v = i / en;
			for (int j = 0; j < number_x; j++)
			{
				index_ij = i*number_x + j;
				real = out_c[index_ij][0] * tp*v;
				imag = out_c[index_ij][1] * tp*v;
				out_c[index_ij][0] = -imag;
				out_c[index_ij][1] = real;
			}
		}
		for (int i = number_y_half + 1; i <number_y; i++)
		{
			v = (i - number_y) / en;
			for (int j = 0; j < number_x; j++)
			{
				index_ij = i*number_x + j;
				real = out_c[index_ij][0] * tp*v;
				imag = out_c[index_ij][1] * tp*v;
				out_c[index_ij][0] = -imag;
				out_c[index_ij][1] = real;
			}
		}
		IFFT2d(number_y, number_x, out_c, processdata->GetData());

		return 1;
	}
	int DerivativeHorizontal(CGM_GRDDATA indata, double* Orientation, CGM_GRDDATA* processdata)
	{
		double cosx = Orientation[0], cosy = Orientation[1];
		double tp = 2 * PI;
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double* originaldata = indata.GetData();
		fftw_complex* out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
		FFT2d(number_y, number_x, originaldata, out_c);			//傅里叶变换
		double *bounds = indata.GetBounds();
		double dy = fabs(bounds[3] - bounds[2]) / (number_y - 1), dx = fabs(bounds[1] - bounds[0]) / (number_x - 1);
		double en = dy*(number_y), em = dx*(number_x);
		double u, v;
		double real, imag;
		int number_x_half = (int)(number_x / 2);
		int number_y_half = (int)(number_y / 2);
		int index_ij;
		//乘以频率域因子改变傅里叶变换系数
		//1. v为正
		for (int i = 0; i <= number_y_half; i++)
		{
			u = i / en;
			for (int j = 0; j <= number_x_half; j++)
			{
				v = j / em;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0] * tp*(u*cosx + v*cosy);
				imag = out_c[index_ij][1] * tp*(u*cosx + v*cosy);
				out_c[index_ij][0] = -imag;
				out_c[index_ij][1] = real;
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				v = (j - number_x) / em;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0] * tp*(u*cosx + v*cosy);
				imag = out_c[index_ij][1] * tp*(u*cosx + v*cosy);
				out_c[index_ij][0] = -imag;
				out_c[index_ij][1] = real;
			}
		}
		for (int i = number_y_half + 1; i <number_y; i++)
		{
			u = (i - number_y) / en;
			for (int j = 0; j <= number_x_half; j++)
			{
				v = j / em;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0] * tp*(u*cosx + v*cosy);
				imag = out_c[index_ij][1] * tp*(u*cosx + v*cosy);
				out_c[index_ij][0] = -imag;
				out_c[index_ij][1] = real;
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				v = (j - number_x) / em;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0] * tp*(u*cosx + v*cosy);
				imag = out_c[index_ij][1] * tp*(u*cosx + v*cosy);
				out_c[index_ij][0] = -imag;
				out_c[index_ij][1] = real;
			}
		}
		//逆变换
		IFFT2d(number_y, number_x, out_c, processdata->GetData());

		return 1;
	}
	int DerivativeHorizontal2(CGM_GRDDATA indata, double aerfa_x, CGM_GRDDATA* processdata)
	{
		double orientation[2];
		orientation[0] = cos(aerfa_x / 180.0*PI);
		orientation[1] = sin(aerfa_x / 180.0*PI);
		DerivativeHorizontal(indata, orientation, processdata);
		return 1;
	}
	int AmplitudeSpectrum(CGM_GRDDATA indata, CGM_GRDDATA* processdata)
	{
		double tp = 2 * PI;
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double* originaldata = indata.GetData();
		fftw_complex* out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
		FFT2d(number_y, number_x, originaldata, out_c);			//傅里叶变换
		double *bounds = indata.GetBounds();
		double dy = fabs(bounds[3] - bounds[2]) / (number_y - 1), dx = fabs(bounds[1] - bounds[0]) / (number_x - 1);
		double en = dy*(number_y), em = dx*(number_x);
		double u, v;
		double real, imag;
		int number_x_half = (int)(number_x / 2);
		int number_y_half = (int)(number_y / 2);
		int index_ij;
		//改变变换后的数据的xy的范围，为频率
		//double* bounds_uv = processdata->GetBounds();
		//bounds_uv[0] = -number_x_half / em; bounds_uv[1] = number_x_half / em;
		//bounds_uv[2] = -number_y_half / en; bounds_uv[3] = number_y_half / en;
		//bounds_uv[4] = -number_y_half / en; bounds_uv[3] = 0;
		//乘以频率域因子改变傅里叶变换系数
		//给出的振幅谱的图像是分为四个象限-FS~0，0~FS
		//1. v为正
		for (int i = 0; i <= number_y_half; i++)
		{
			v = i / en;
			for (int j = 0; j <= number_x_half; j++)
			{
				u = j / em;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				processdata->GetData()[(i + (number_y - number_y_half - 1))*number_x + j + (number_x - number_x_half - 1)] = sqrt(real*real + imag*imag);
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				u = (j - number_x) / em;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				processdata->GetData()[(i + (number_y - number_y_half - 1))*number_x + j - (number_x_half + 1)] = sqrt(real*real + imag*imag);;
			}
		}
		for (int i = number_y_half + 1; i <number_y; i++)
		{
			v = (i - number_y) / en;
			for (int j = 0; j <= number_x_half; j++)
			{
				u = j / em;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				processdata->GetData()[(i - (number_y_half + 1))*number_x + j + (number_x - number_x_half - 1)] = sqrt(real*real + imag*imag);;
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				u = (j - number_x) / em;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				processdata->GetData()[(i - (number_y_half + 1))*number_x + j - (number_x_half + 1)] = sqrt(real*real + imag*imag);;
			}
		}

		return 1;
	}
	int UpwardContinuation(CGM_GRDDATA indata, double height, CGM_GRDDATA* processdata)
	{
		double tp = 2 * PI;
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double* originaldata = indata.GetData();
		fftw_complex* out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
		FFT2d(number_y, number_x, originaldata, out_c);			//傅里叶变换
		double *bounds = indata.GetBounds();
		double dy = fabs(bounds[3] - bounds[2]) / (number_y - 1), dx = fabs(bounds[1] - bounds[0]) / (number_x - 1);
		double en = dy*(number_y), em = dx*(number_x);
		double u, v;
		double factor = -tp*height;
		double factor_conti;
		int number_x_half = (int)(number_x / 2);
		int number_y_half = (int)(number_y / 2);
		int index_ij;
		//乘以频率域因子改变傅里叶变换系数
		//1. v为正
		for (int i = 0; i <= number_y_half; i++)
		{
			v = i / en;
			for (int j = 0; j <= number_x_half; j++)
			{
				u = j / em;
				index_ij = i*number_x + j;
				factor_conti = exp(factor*sqrt(u*u + v*v));
				out_c[index_ij][0] = out_c[index_ij][0] * factor_conti;
				out_c[index_ij][1] = out_c[index_ij][1] * factor_conti;
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				u = (j - number_x) / em;
				index_ij = i*number_x + j;
				factor_conti = exp(factor*sqrt(u*u + v*v));
				out_c[index_ij][0] = out_c[index_ij][0] * factor_conti;
				out_c[index_ij][1] = out_c[index_ij][1] * factor_conti;
			}
		}
		for (int i = number_y_half + 1; i <number_y; i++)
		{
			v = (i - number_y) / en;
			for (int j = 0; j <= number_x_half; j++)
			{
				u = j / em;
				index_ij = i*number_x + j;
				factor_conti = exp(factor*sqrt(u*u + v*v));
				out_c[index_ij][0] = out_c[index_ij][0] * factor_conti;
				out_c[index_ij][1] = out_c[index_ij][1] * factor_conti;
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				u = (j - number_x) / em;
				index_ij = i*number_x + j;
				factor_conti = exp(factor*sqrt(u*u + v*v));
				out_c[index_ij][0] = out_c[index_ij][0] * factor_conti;
				out_c[index_ij][1] = out_c[index_ij][1] * factor_conti;
			}
		}
		//逆变换
		IFFT2d(number_y, number_x, out_c, processdata->GetData());

		return 1;
	}

	int DerivativeVertical(CGM_GRDDATA indata, CGM_GRDDATA* processdata)
	{
		double tp = 2 * PI;
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double* originaldata = indata.GetData();
		fftw_complex* out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
		FFT2d(number_y, number_x, originaldata, out_c);			//傅里叶变换
		double *bounds = indata.GetBounds();
		double dy = fabs(bounds[3] - bounds[2]) / (number_y - 1), dx = fabs(bounds[1] - bounds[0]) / (number_x - 1);
		double en = dy*(number_y), em = dx*(number_x);
		double u, v;
		double factor_vertical;
		int number_x_half = (int)(number_x / 2);
		int number_y_half = (int)(number_y / 2);
		int index_ij;
		//乘以频率域因子改变傅里叶变换系数
		//1. v为正
		for (int i = 0; i <= number_y_half; i++)
		{
			v = i / en;
			for (int j = 0; j <= number_x_half; j++)
			{
				u = j / em;
				index_ij = i*number_x + j;
				factor_vertical = sqrt(u*u + v*v)*tp;
				out_c[index_ij][0] = out_c[index_ij][0] * factor_vertical;
				out_c[index_ij][1] = out_c[index_ij][1] * factor_vertical;
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				u = (j - number_x) / em;
				index_ij = i*number_x + j;
				factor_vertical = sqrt(u*u + v*v)*tp;
				out_c[index_ij][0] = out_c[index_ij][0] * factor_vertical;
				out_c[index_ij][1] = out_c[index_ij][1] * factor_vertical;
			}
		}
		for (int i = number_y_half + 1; i <number_y; i++)
		{
			v = (i - number_y) / en;
			for (int j = 0; j <= number_x_half; j++)
			{
				u = j / em;
				factor_vertical = sqrt(u*u + v*v)*tp;
				index_ij = i*number_x + j;
				out_c[index_ij][0] = out_c[index_ij][0] * factor_vertical;
				out_c[index_ij][1] = out_c[index_ij][1] * factor_vertical;
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				u = (j - number_x) / em;
				index_ij = i*number_x + j;
				factor_vertical = sqrt(u*u + v*v)*tp;
				out_c[index_ij][0] = out_c[index_ij][0] * factor_vertical;
				out_c[index_ij][1] = out_c[index_ij][1] * factor_vertical;
			}
		}
		//逆变换
		IFFT2d(number_y, number_x, out_c, processdata->GetData());

		return 1;
	}

	int Derivative(CGM_GRDDATA indata, double* Orientation, CGM_GRDDATA* processdata)
	{
		double cosx = Orientation[0], cosy = Orientation[1], cosz = Orientation[2];
		double tp = 2.0 * PI;
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double* originaldata = indata.GetData();
		fftw_complex* out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
		FFT2d(number_y, number_x, originaldata, out_c);			//傅里叶变换
		double *bounds = indata.GetBounds();
		double dy = fabs(bounds[3] - bounds[2]) / (number_y - 1), dx = fabs(bounds[1] - bounds[0]) / (number_x - 1);
		double en = dy*(number_y), em = dx*(number_x);
		double u, v;
		double real, imag;
		double factor_vertical;
		int number_x_half = (int)(number_x / 2);
		int number_y_half = (int)(number_y / 2);
		int index_ij;
		//乘以频率域因子改变傅里叶变换系数
		//1. v为正
		for (int i = 0; i <= number_y_half; i++)
		{
			v = i / en;
			for (int j = 0; j <= number_x_half; j++)
			{
				u = j / em;
				index_ij = i*number_x + j;
				factor_vertical = sqrt(u*u + v*v);
				real = out_c[index_ij][0] * tp*(u*cosx + v*cosy + factor_vertical*cosz);
				imag = out_c[index_ij][1] * tp*(u*cosx + v*cosy + factor_vertical*cosz);
				out_c[index_ij][0] = -imag;
				out_c[index_ij][1] = real;
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				u = (j - number_x) / em;
				index_ij = i*number_x + j;
				factor_vertical = sqrt(u*u + v*v);
				real = out_c[index_ij][0] * tp*(u*cosx + v*cosy + factor_vertical*cosz);
				imag = out_c[index_ij][1] * tp*(u*cosx + v*cosy + factor_vertical*cosz);
				out_c[index_ij][0] = -imag;
				out_c[index_ij][1] = real;
			}
		}
		for (int i = number_y_half + 1; i <number_y; i++)
		{
			v = (i - number_y) / en;
			for (int j = 0; j <= number_x_half; j++)
			{
				u = j / em;
				index_ij = i*number_x + j;
				factor_vertical = sqrt(u*u + v*v);
				real = out_c[index_ij][0] * tp*(u*cosx + v*cosy + factor_vertical*cosz);
				imag = out_c[index_ij][1] * tp*(u*cosx + v*cosy + factor_vertical*cosz);
				out_c[index_ij][0] = -imag;
				out_c[index_ij][1] = real;
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				u = (j - number_x) / em;
				index_ij = i*number_x + j;
				factor_vertical = sqrt(u*u + v*v);
				real = out_c[index_ij][0] * tp*(u*cosx + v*cosy + factor_vertical*cosz);
				imag = out_c[index_ij][1] * tp*(u*cosx + v*cosy + factor_vertical*cosz);
				out_c[index_ij][0] = -imag;
				out_c[index_ij][1] = real;
			}
		}
		//逆变换
		IFFT2d(number_y, number_x, out_c, processdata->GetData());

		return 1;
	}

	int Derivative2(CGM_GRDDATA indata, double* Angle, CGM_GRDDATA* processdata)
	{
		double angle_A = Angle[0], angle_I = Angle[1];
		double orientation[3];
		angle_A = angle_A / 180.0*PI;
		angle_I = angle_I / 180.0*PI;
		orientation[0] = cos(angle_I)*cos(angle_A);
		orientation[1] = cos(angle_I)*sin(angle_A);
		orientation[2] = sin(angle_I);
		Derivative(indata, orientation, processdata);

		return 1;
	}

	double* MeanRadialPowerSpectrum(CGM_GRDDATA indata, int& frequenceNum)
	{
		//--------------------
		double tp = 2 * PI;
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double* originaldata = indata.GetData();
		fftw_complex* out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
		FFT2d(number_y, number_x, originaldata, out_c);			//傅里叶变换
		double *bounds = indata.GetBounds();
		double dy = fabs(bounds[3] - bounds[2]) / (number_y - 1), dx = fabs(bounds[1] - bounds[0]) / (number_x - 1);
		double en = dy*(number_y), em = dx*(number_x);
		double u, v;
		int number_x_half = (int)(number_x / 2);
		int number_y_half = (int)(number_y / 2);
		double df = sqrt(1.0 / (em*em) + 1.0 / (en*en));
		double fmax = sqrt(dx*dx + dy*dy) / dx / dy / 2.0;
		frequenceNum = (int)(fmax / df + 0.5);
		double spectrum = 0, sum = 0, kk = 0, r = 0;
		double*spectrumTab = new double[frequenceNum * 2];
		for (int k = 1; k < frequenceNum; k++)
		{
			spectrum = 0; sum = 0; kk = 0;
			for (int i = 0; i < number_y_half; i++)
			{
				v = i / en;
				for (int j = 0; j < number_x_half; j++)
				{
					u = j / em;
					r = sqrt(u*u + v*v);
					if ((r>df*(k - 1)) && (r <= df*k))
					{
						sum += out_c[j + i*number_x][0] * out_c[j + i*number_x][0] + out_c[j + i*number_x][1] * out_c[j + i*number_x][1];
						kk++;
					}
				}
			}

			if (kk != 0)
			{
				spectrum = sum / kk;
				spectrumTab[k] = df*k;
				spectrumTab[frequenceNum + k] = spectrum;
			}
		}
		spectrumTab[0 + frequenceNum] = 0;
		spectrumTab[0] = out_c[0][0] * out_c[0][0] + out_c[0][1] * out_c[0][1];

		return spectrumTab;
	}

	//释放此DLL内部申请的内存空间（1唯情况
	int DeleteArray1D(double* array1d)
	{
		if (array1d)
		{
			delete[] array1d;

		}
		return 1;
	}
	int MagComponentTransT2X_ext(CGM_GRDDATA indata, MagneticComponentTransStruct* magdatainfo, CGM_GRDDATA* processdata)
	{
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		double *bounds = indata.GetBounds();
		//扩边
		CGM_GRDDATA indata_ext;
		int extenNum_x = number_x * 2, extenNum_y = number_y * 2;
		indata_ext.Initialize(extenNum_x, extenNum_y, bounds);
		ExtenBoundary(indata, &indata_ext);
		//处理后的数据存储
		CGM_GRDDATA processdata_ext;
		double *bounds_ext = indata_ext.GetBounds();
		processdata_ext.Initialize(extenNum_x, extenNum_y, bounds_ext);
		//调用函数
		MagComponentTransT2X(indata, magdatainfo, processdata);

		////缩边
		//double* outdata = processdata->GetData();
		//int num_ext_left = (extenNum_x-number_x) / 2;
		//int num_ext_down = (extenNum_y-number_y) / 2;
		//for (int i = 0; i < number_y; i++)
		//{
		//	for (int j = 0; j < number_x; j++)
		//	{
		//		outdata[j + i*number_x] = processdata_ext.GetData()[(j + num_ext_left) + (i + num_ext_down)*extenNum_x];
		//	}
		//}

		////销毁指针
		//indata_ext.Delete();
		//processdata_ext.Delete();
		return 1;
	}
	int MagComponentTransT2X(CGM_GRDDATA indata, MagneticComponentTransStruct* magdatainfo, CGM_GRDDATA* processdata)
	{
		double angle_I = magdatainfo->MagAngle_I_Earth / 180.0*PI;
		double angle_D = magdatainfo->MagAngle_D_Earth / 180.0*PI;
		double L0 = cos(angle_I)*cos(angle_D);
		double M0 = cos(angle_I)*sin(angle_D);
		double N0 = sin(angle_I);

		double tp = 2 * PI;
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double *bounds = indata.GetBounds();
		
		//利用扩边后的数据进行计算
		double* originaldata = indata.GetData();
		fftw_complex* out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
		FFT2d(number_y, number_x, originaldata, out_c);			//傅里叶变换
		
		double dy = fabs(bounds[3] - bounds[2]) / (number_y - 1), dx = fabs(bounds[1] - bounds[0]) / (number_x - 1);
		double en = dy*(number_y), em = dx*(number_x);
		double u, v;
		double real, imag;
		int number_x_half = (int)(number_x / 2);
		int number_y_half = (int)(number_y / 2);
		double factor_real, factor_imag;
		double a, b, c;
		int index_ij;
		//乘以频率域因子改变傅里叶变换系数
		//程序设计中列方向为y方向，而列方向又对应的是北向，在书本中的北向其实才是x方向，因此为了不使结果错误，再循环中让u为列方向的频率，v为行方向的频率，那么
		//书本中的公式就可以直接使用而不会出现错误
		//1. T→X
		for (int i = 0; i <= number_y_half; i++)
		{
			u = i / en + 1E-10;
			for (int j = 0; j <= number_x_half; j++)
			{
				v = j / em + 1E-10;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				b = (L0*u + M0*v);
				a = N0*sqrt(u*u+v*v);
				c = a*a + b*b;
				factor_real = b / c*u;
				factor_imag = a / c*u;

				out_c[index_ij][0] = (factor_real*real-factor_imag*imag);
				out_c[index_ij][1] = (factor_real*imag+factor_imag*real);
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				v = (j - number_x) / em + 1E-10;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				b = (L0*u + M0*v);
				a = N0*sqrt(u*u + v*v);
				c = a*a + b*b;
				factor_real = b / c*u;
				factor_imag = a / c*u;

				out_c[index_ij][0] = (factor_real*real - factor_imag*imag);
				out_c[index_ij][1] = (factor_real*imag + factor_imag*real);
			}
		}
		for (int i = number_y_half + 1; i <number_y; i++)
		{
			u = (i - number_y) / en + 1E-10;
			for (int j = 0; j <= number_x_half; j++)
			{
				v = j / em + 1E-10;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				b = (L0*u + M0*v);
				a = N0*sqrt(u*u + v*v);
				c = a*a + b*b;
				factor_real = b / c*u;
				factor_imag = a / c*u;

				out_c[index_ij][0] = (factor_real*real - factor_imag*imag);
				out_c[index_ij][1] = (factor_real*imag + factor_imag*real);
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				v = (j - number_x) / em + 1E-10;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				b = (L0*u + M0*v);
				a = N0*sqrt(u*u + v*v);
				c = a*a + b*b;
				factor_real = b / c*u;
				factor_imag = a / c*u;

				out_c[index_ij][0] = (factor_real*real - factor_imag*imag);
				out_c[index_ij][1] = (factor_real*imag + factor_imag*real);
			}
		}
		//逆变换
		IFFT2d(number_y, number_x, out_c, processdata->GetData());
		return 1;
	}

	int MagComponentTransT2Y(CGM_GRDDATA indata, MagneticComponentTransStruct* magdatainfo, CGM_GRDDATA* processdata)
	{
		double angle_I = magdatainfo->MagAngle_I_Earth / 180.0*PI;
		double angle_D = magdatainfo->MagAngle_D_Earth / 180.0*PI;
		double L0 = cos(angle_I)*cos(angle_D);
		double M0 = cos(angle_I)*sin(angle_D);
		double N0 = sin(angle_I);

		double tp = 2 * PI;
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double* originaldata = indata.GetData();
		fftw_complex* out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
		FFT2d(number_y, number_x, originaldata, out_c);			//傅里叶变换
		double *bounds = indata.GetBounds();
		double dy = fabs(bounds[3] - bounds[2]) / (number_y - 1), dx = fabs(bounds[1] - bounds[0]) / (number_x - 1);
		double en = dy*(number_y), em = dx*(number_x);
		double u, v;
		double real, imag;
		int number_x_half = (int)(number_x / 2);
		int number_y_half = (int)(number_y / 2);
		double factor_real, factor_imag;
		double a, b, c;
		int index_ij;
		//乘以频率域因子改变傅里叶变换系数
		//程序设计中列方向为y方向，而列方向又对应的是北向，在书本中的北向其实才是x方向，因此为了不使结果错误，再循环中让u为列方向的频率，v为行方向的频率，那么
		//书本中的公式就可以直接使用而不会出现错误
		//1. T→X
		for (int i = 0; i <= number_y_half; i++)
		{
			u = i / en + 1E-10;
			for (int j = 0; j <= number_x_half; j++)
			{
				v = j / em + 1E-10;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				b = (L0*u + M0*v);
				a = N0*sqrt(u*u + v*v);
				c = a*a + b*b;
				factor_real = b / c*v;
				factor_imag = a / c*v;

				out_c[index_ij][0] = (factor_real*real - factor_imag*imag);
				out_c[index_ij][1] = (factor_real*imag + factor_imag*real);
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				v = (j - number_x) / em + 1E-10;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				b = (L0*u + M0*v);
				a = N0*sqrt(u*u + v*v);
				c = a*a + b*b;
				factor_real = b / c*v;
				factor_imag = a / c*v;

				out_c[index_ij][0] = (factor_real*real - factor_imag*imag);
				out_c[index_ij][1] = (factor_real*imag + factor_imag*real);
			}
		}
		for (int i = number_y_half + 1; i <number_y; i++)
		{
			u = (i - number_y) / en + 1E-10;
			for (int j = 0; j <= number_x_half; j++)
			{
				v = j / em + 1E-10;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				b = (L0*u + M0*v);
				a = N0*sqrt(u*u + v*v);
				c = a*a + b*b;
				factor_real = b / c*v;
				factor_imag = a / c*v;

				out_c[index_ij][0] = (factor_real*real - factor_imag*imag);
				out_c[index_ij][1] = (factor_real*imag + factor_imag*real);
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				v = (j - number_x) / em + 1E-10;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				b = (L0*u + M0*v);
				a = N0*sqrt(u*u + v*v);
				c = a*a + b*b;
				factor_real = b / c*v;
				factor_imag = a / c*v;

				out_c[index_ij][0] = (factor_real*real - factor_imag*imag);
				out_c[index_ij][1] = (factor_real*imag + factor_imag*real);
			}
		}
		//逆变换
		IFFT2d(number_y, number_x, out_c, processdata->GetData());

		return 1;
	}
	int MagComponentTransT2Z(CGM_GRDDATA indata, MagneticComponentTransStruct* magdatainfo, CGM_GRDDATA* processdata)
	{
		double angle_I = magdatainfo->MagAngle_I_Earth / 180.0*PI;
		double angle_D = magdatainfo->MagAngle_D_Earth / 180.0*PI;
		double L0 = cos(angle_I)*cos(angle_D);
		double M0 = cos(angle_I)*sin(angle_D);
		double N0 = sin(angle_I);

		double tp = 2 * PI;
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double* originaldata = indata.GetData();
		fftw_complex* out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
		FFT2d(number_y, number_x, originaldata, out_c);			//傅里叶变换
		double *bounds = indata.GetBounds();
		double dy = fabs(bounds[3] - bounds[2]) / (number_y - 1), dx = fabs(bounds[1] - bounds[0]) / (number_x - 1);
		double en = dy*(number_y), em = dx*(number_x);
		double u, v;
		double real, imag;
		int number_x_half = (int)(number_x / 2);
		int number_y_half = (int)(number_y / 2);
		double factor_real, factor_imag;
		double a, b, c,a0;
		int index_ij;
		//乘以频率域因子改变傅里叶变换系数
		//程序设计中列方向为y方向，而列方向又对应的是北向，在书本中的北向其实才是x方向，因此为了不使结果错误，再循环中让u为列方向的频率，v为行方向的频率，那么
		//书本中的公式就可以直接使用而不会出现错误
		//1. T→X
		for (int i = 0; i <= number_y_half; i++)
		{
			u = i / en + 1E-10;
			for (int j = 0; j <= number_x_half; j++)
			{
				v = j / em + 1E-10;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				b = (L0*u + M0*v);
				a0 = sqrt(u*u + v*v);
				a = N0*a0;
				c = a*a + b*b;
				factor_real = a / c*a0;
				factor_imag = -b / c*a0;

				out_c[index_ij][0] = (factor_real*real - factor_imag*imag);
				out_c[index_ij][1] = (factor_real*imag + factor_imag*real);
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				v = (j - number_x) / em + 1E-10;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				b = (L0*u + M0*v);
				a0 = sqrt(u*u + v*v);
				a = N0*a0;
				c = a*a + b*b;
				factor_real = a / c*a0;
				factor_imag = -b / c*a0;

				out_c[index_ij][0] = (factor_real*real - factor_imag*imag);
				out_c[index_ij][1] = (factor_real*imag + factor_imag*real);
			}
		}
		for (int i = number_y_half + 1; i <number_y; i++)
		{
			u = (i - number_y) / en + 1E-10;
			for (int j = 0; j <= number_x_half; j++)
			{
				v = j / em + 1E-10;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				b = (L0*u + M0*v);
				a0 = sqrt(u*u + v*v);
				a = N0*a0;
				c = a*a + b*b;
				factor_real = a / c*a0;
				factor_imag = -b / c*a0;

				out_c[index_ij][0] = (factor_real*real - factor_imag*imag);
				out_c[index_ij][1] = (factor_real*imag + factor_imag*real);
			}

			for (int j = number_x_half + 1; j < number_x; j++)
			{
				v = (j - number_x) / em + 1E-10;
				index_ij = i*number_x + j;
				real = out_c[index_ij][0];
				imag = out_c[index_ij][1];
				b = (L0*u + M0*v);
				a0 = sqrt(u*u + v*v);
				a = N0*a0;
				c = a*a + b*b;
				factor_real = a / c*a0;
				factor_imag = -b / c*a0;

				out_c[index_ij][0] = (factor_real*real - factor_imag*imag);
				out_c[index_ij][1] = (factor_real*imag + factor_imag*real);
			}
		}
		//逆变换
		IFFT2d(number_y, number_x, out_c, processdata->GetData());

		return 1;
	}

	int MagComponentTransT2Ta(CGM_GRDDATA indata, MagneticComponentTransStruct* magdatainfo, CGM_GRDDATA* processdata)
	{

		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double* originaldata = indata.GetData();
		double *bounds = indata.GetBounds();
		//首先计算三分量
		//1.1 计算hax
		CGM_GRDDATA processout_Hax;
		processout_Hax.Initialize(number_x, number_y, bounds);
		MagComponentTransT2X(indata, magdatainfo, &processout_Hax);
		//1.2 计算hay
		CGM_GRDDATA processout_Hay;
		processout_Hay.Initialize(number_x, number_y, bounds);
		MagComponentTransT2Y(indata, magdatainfo, &processout_Hay);
		//1.3 计算Za
		CGM_GRDDATA processout_Za;
		processout_Za.Initialize(number_x, number_y, bounds);
		MagComponentTransT2Z(indata, magdatainfo, &processout_Za);

		//2.0计算模量
		double* outdata = processdata->GetData();
		double* Hax = processout_Hax.GetData();
		double* Hay = processout_Hay.GetData();
		double* Za = processout_Za.GetData();
		for (int i = 0; i < number_y; i++)
		{
			for (int j = 0; j < number_x; j++)
			{
				outdata[j + i*number_x] = sqrt(Hax[j + i*number_x] * Hax[j + i*number_x] + Hay[j + i*number_x] * Hay[j + i*number_x] + Za[j + i*number_x] * Za[j + i*number_x]);
			}
		}

		////释放内存
		processout_Hax.Delete();
		processout_Hay.Delete();
		processout_Za.Delete();
		return 1;
	}

	int MagComponentTransT2R(CGM_GRDDATA indata, MagneticComponentTransStruct* magdatainfo, CGM_GRDDATA* processdata)
	{
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double* originaldata = indata.GetData();
		double *bounds = indata.GetBounds();
		//首先计算三分量
		//1.1 计算hax
		CGM_GRDDATA processout_Hax;
		processout_Hax.Initialize(number_x, number_y, bounds);
		MagComponentTransT2X(indata, magdatainfo, &processout_Hax);
		//1.2 计算hay
		CGM_GRDDATA processout_Hay;
		processout_Hay.Initialize(number_x, number_y, bounds);
		MagComponentTransT2Y(indata, magdatainfo, &processout_Hay);
		//1.3 计算Za
		CGM_GRDDATA processout_Za;
		processout_Za.Initialize(number_x, number_y, bounds);
		MagComponentTransT2Z(indata, magdatainfo, &processout_Za);
		//1.4 计算三分量的三个方向的导数
		CGM_GRDDATA processout_Hax_x;//Hax的导数
		processout_Hax_x.Initialize(number_x, number_y, bounds);
		DerivativeX(processout_Hax, &processout_Hax_x);
		CGM_GRDDATA processout_Hax_y;
		processout_Hax_y.Initialize(number_x, number_y, bounds);
		DerivativeY(processout_Hax, &processout_Hax_y);
		CGM_GRDDATA processout_Hax_z;
		processout_Hax_z.Initialize(number_x, number_y, bounds);
		DerivativeVertical(processout_Hax, &processout_Hax_z);
		CGM_GRDDATA processout_Hay_x;	//Hay的导数
		processout_Hay_x.Initialize(number_x, number_y, bounds);
		DerivativeX(processout_Hay, &processout_Hay_x);
		CGM_GRDDATA processout_Hay_y;
		processout_Hay_y.Initialize(number_x, number_y, bounds);
		DerivativeY(processout_Hay, &processout_Hay_y);
		CGM_GRDDATA processout_Hay_z;
		processout_Hay_z.Initialize(number_x, number_y, bounds);
		DerivativeVertical(processout_Hay, &processout_Hay_z);
		CGM_GRDDATA processout_Za_x;	//Za的导数
		processout_Za_x.Initialize(number_x, number_y, bounds);
		DerivativeX(processout_Za, &processout_Za_x);
		CGM_GRDDATA processout_Za_y;	//Za的导数
		processout_Za_y.Initialize(number_x, number_y, bounds);
		DerivativeY(processout_Za, &processout_Za_y);
		CGM_GRDDATA processout_Za_z;	//Za的导数
		processout_Za_z.Initialize(number_x, number_y, bounds);
		DerivativeVertical(processout_Za, &processout_Za_z);

		//2.0计算模量
		double* outdata = processdata->GetData();
		double* Hax = processout_Hax.GetData();
		double* Hay = processout_Hay.GetData();
		double* Za = processout_Za.GetData();
		double* Za_z=processout_Za_z.GetData();
		double* Za_x=processout_Za_x.GetData();
		double* Za_y = processout_Za_y.GetData();
		double* Hax_x = processout_Hax_x.GetData();
		double* Hax_y = processout_Hax_y.GetData();
		double* Hax_z = processout_Hax_z.GetData();
		double* Hay_x = processout_Hay_x.GetData();
		double* Hay_y = processout_Hay_y.GetData();
		double* Hay_z = processout_Hay_z.GetData();
		double Ta,Tax, Tay, Taz;
		for (int i = 0; i < number_y; i++)
		{
			for (int j = 0; j < number_x; j++)
			{
				Ta = sqrt(Hax[j + i*number_x] * Hax[j + i*number_x] + Hay[j + i*number_x] * Hay[j + i*number_x] + Za[j + i*number_x] * Za[j + i*number_x]);
				Tax= (Hax[j + i*number_x] * Hax_x[j + i*number_x] + Hay[j + i*number_x] * Hay_x[j + i*number_x] + Za[j + i*number_x] * Za_x[j + i*number_x]) / Ta;
				Tay = (Hax[j + i*number_x] * Hax_y[j + i*number_x] + Hay[j + i*number_x] * Hay_y[j + i*number_x] + Za[j + i*number_x] * Za_y[j + i*number_x]) / Ta;
				Taz = (Hax[j + i*number_x] * Hax_z[j + i*number_x] + Hay[j + i*number_x] * Hay_z[j + i*number_x] + Za[j + i*number_x] * Za_z[j + i*number_x]) / Ta;
				outdata[j + i*number_x] = sqrt(Tax * Tax + Tay * Tay + Taz * Taz);
			}
		}

		////释放内存
		processout_Hax.Delete();
		processout_Hay.Delete();
		processout_Za.Delete();

		processout_Za_z.Delete();
		processout_Za_x.Delete();
		processout_Za_y.Delete();
		processout_Hax_x.Delete();
		processout_Hax_y.Delete();
		processout_Hax_z.Delete();
		processout_Hay_x.Delete();
		processout_Hay_y.Delete();
		processout_Hay_z.Delete();

		return 1;
	}
	int MagComponentTransT2E(CGM_GRDDATA indata, MagneticComponentTransStruct* magdatainfo, CGM_GRDDATA* processdata)
	{
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double* originaldata = indata.GetData();
		double *bounds = indata.GetBounds();
		//首先计算三分量
		//1.1 计算hax
		CGM_GRDDATA processout_Hax;
		processout_Hax.Initialize(number_x, number_y, bounds);
		MagComponentTransT2X(indata, magdatainfo, &processout_Hax);
		//1.2 计算hay
		CGM_GRDDATA processout_Hay;
		processout_Hay.Initialize(number_x, number_y, bounds);
		MagComponentTransT2Y(indata, magdatainfo, &processout_Hay);
		//1.3 计算Za
		CGM_GRDDATA processout_Za;
		processout_Za.Initialize(number_x, number_y, bounds);
		MagComponentTransT2Z(indata, magdatainfo, &processout_Za);
		//1.4 计算三分量的三个方向的导数
		CGM_GRDDATA processout_Hax_x;//Hax的导数
		processout_Hax_x.Initialize(number_x, number_y, bounds);
		DerivativeX(processout_Hax, &processout_Hax_x);
		CGM_GRDDATA processout_Hax_y;
		processout_Hax_y.Initialize(number_x, number_y, bounds);
		DerivativeY(processout_Hax, &processout_Hax_y);
		CGM_GRDDATA processout_Hax_z;
		processout_Hax_z.Initialize(number_x, number_y, bounds);
		DerivativeVertical(processout_Hax, &processout_Hax_z);
		CGM_GRDDATA processout_Hay_x;	//Hay的导数
		processout_Hay_x.Initialize(number_x, number_y, bounds);
		DerivativeX(processout_Hay, &processout_Hay_x);
		CGM_GRDDATA processout_Hay_y;
		processout_Hay_y.Initialize(number_x, number_y, bounds);
		DerivativeY(processout_Hay, &processout_Hay_y);
		CGM_GRDDATA processout_Hay_z;
		processout_Hay_z.Initialize(number_x, number_y, bounds);
		DerivativeVertical(processout_Hay, &processout_Hay_z);
		CGM_GRDDATA processout_Za_x;	//Za的导数
		processout_Za_x.Initialize(number_x, number_y, bounds);
		DerivativeX(processout_Za, &processout_Za_x);
		CGM_GRDDATA processout_Za_y;	//Za的导数
		processout_Za_y.Initialize(number_x, number_y, bounds);
		DerivativeY(processout_Za, &processout_Za_y);
		CGM_GRDDATA processout_Za_z;	//Za的导数
		processout_Za_z.Initialize(number_x, number_y, bounds);
		DerivativeVertical(processout_Za, &processout_Za_z);

		//2.0计算模量
		double* outdata = processdata->GetData();
		double* Hax = processout_Hax.GetData();
		double* Hay = processout_Hay.GetData();
		double* Za = processout_Za.GetData();
		double* Za_z = processout_Za_z.GetData();
		double* Za_x = processout_Za_x.GetData();
		double* Za_y = processout_Za_y.GetData();
		double* Hax_x = processout_Hax_x.GetData();
		double* Hax_y = processout_Hax_y.GetData();
		double* Hax_z = processout_Hax_z.GetData();
		double* Hay_x = processout_Hay_x.GetData();
		double* Hay_y = processout_Hay_y.GetData();
		double* Hay_z = processout_Hay_z.GetData();
		double Theta_Hx, Theta_Hy, Theta_Za;
		for (int i = 0; i < number_y; i++)
		{
			for (int j = 0; j < number_x; j++)
			{
				Theta_Hx = Hax_x[j + i*number_x] * Hax_x[j + i*number_x] + Hax_y[j + i*number_x] * Hax_y[j + i*number_x] + Hax_z[j + i*number_x] * Hax_z[j + i*number_x];
				Theta_Hy = Hay_x[j + i*number_x] * Hay_x[j + i*number_x] + Hay_y[j + i*number_x] * Hay_y[j + i*number_x] + Hay_z[j + i*number_x] * Hay_z[j + i*number_x];
				Theta_Za = Za_x[j + i*number_x] * Za_x[j + i*number_x] + Za_y[j + i*number_x] * Za_y[j + i*number_x] + Za_z[j + i*number_x] * Za_z[j + i*number_x];
				outdata[j + i*number_x] = sqrt((Theta_Hx + Theta_Hy + Theta_Za) / 2.0);
			}
		}

		////释放内存
		processout_Hax.Delete();
		processout_Hay.Delete();
		processout_Za.Delete();

		processout_Za_z.Delete();
		processout_Za_x.Delete();
		processout_Za_y.Delete();
		processout_Hax_x.Delete();
		processout_Hax_y.Delete();
		processout_Hax_z.Delete();
		processout_Hay_x.Delete();
		processout_Hay_y.Delete();
		processout_Hay_z.Delete();

		return 1;
	}
	int MagComponentTransT2Q(CGM_GRDDATA indata, MagneticComponentTransStruct* magdatainfo, CGM_GRDDATA* processdata)
	{
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double* originaldata = indata.GetData();
		double *bounds = indata.GetBounds();
		//首先计算三分量
		//1.1 计算hax
		CGM_GRDDATA processout_Hax;
		processout_Hax.Initialize(number_x, number_y, bounds);
		MagComponentTransT2X(indata, magdatainfo, &processout_Hax);
		//1.2 计算hay
		CGM_GRDDATA processout_Hay;
		processout_Hay.Initialize(number_x, number_y, bounds);
		MagComponentTransT2Y(indata, magdatainfo, &processout_Hay);
		//1.3 计算Za
		CGM_GRDDATA processout_Za;
		processout_Za.Initialize(number_x, number_y, bounds);
		MagComponentTransT2Z(indata, magdatainfo, &processout_Za);
		//1.4 计算三分量的三个方向的导数
		CGM_GRDDATA processout_Hax_x;//Hax的导数
		processout_Hax_x.Initialize(number_x, number_y, bounds);
		DerivativeX(processout_Hax, &processout_Hax_x);
		CGM_GRDDATA processout_Hax_y;
		processout_Hax_y.Initialize(number_x, number_y, bounds);
		DerivativeY(processout_Hax, &processout_Hax_y);
		CGM_GRDDATA processout_Hax_z;
		processout_Hax_z.Initialize(number_x, number_y, bounds);
		DerivativeVertical(processout_Hax, &processout_Hax_z);
		CGM_GRDDATA processout_Hay_x;	//Hay的导数
		processout_Hay_x.Initialize(number_x, number_y, bounds);
		DerivativeX(processout_Hay, &processout_Hay_x);
		CGM_GRDDATA processout_Hay_y;
		processout_Hay_y.Initialize(number_x, number_y, bounds);
		DerivativeY(processout_Hay, &processout_Hay_y);
		CGM_GRDDATA processout_Hay_z;
		processout_Hay_z.Initialize(number_x, number_y, bounds);
		DerivativeVertical(processout_Hay, &processout_Hay_z);
		CGM_GRDDATA processout_Za_x;	//Za的导数
		processout_Za_x.Initialize(number_x, number_y, bounds);
		DerivativeX(processout_Za, &processout_Za_x);
		CGM_GRDDATA processout_Za_y;	//Za的导数
		processout_Za_y.Initialize(number_x, number_y, bounds);
		DerivativeY(processout_Za, &processout_Za_y);
		CGM_GRDDATA processout_Za_z;	//Za的导数
		processout_Za_z.Initialize(number_x, number_y, bounds);
		DerivativeVertical(processout_Za, &processout_Za_z);

		//2.0计算模量
		double* outdata = processdata->GetData();
		double* Hax = processout_Hax.GetData();
		double* Hay = processout_Hay.GetData();
		double* Za = processout_Za.GetData();
		double* Za_z = processout_Za_z.GetData();
		double* Za_x = processout_Za_x.GetData();
		double* Za_y = processout_Za_y.GetData();
		double* Hax_x = processout_Hax_x.GetData();
		double* Hax_y = processout_Hax_y.GetData();
		double* Hax_z = processout_Hax_z.GetData();
		double* Hay_x = processout_Hay_x.GetData();
		double* Hay_y = processout_Hay_y.GetData();
		double* Hay_z = processout_Hay_z.GetData();
		double Theta_Hx, Theta_Hy, Theta_Za;
		double Ta, Tax, Tay, Taz;
		double E, R;
		int index_ij;
		for (int i = 0; i < number_y; i++)
		{
			for (int j = 0; j < number_x; j++)
			{
				index_ij = j + i*number_x;
				Ta = sqrt(Hax[index_ij] * Hax[index_ij] + Hay[index_ij] * Hay[index_ij] + Za[index_ij] * Za[index_ij]);
				Tax = (Hax[index_ij] * Hax_x[index_ij] + Hay[index_ij] * Hay_x[index_ij] + Za[index_ij] * Za_x[index_ij]) / Ta;
				Tay = (Hax[index_ij] * Hax_y[index_ij] + Hay[index_ij] * Hay_y[index_ij] + Za[index_ij] * Za_y[index_ij]) / Ta;
				Taz = (Hax[index_ij] * Hax_z[index_ij] + Hay[index_ij] * Hay_z[index_ij] + Za[index_ij] * Za_z[index_ij]) / Ta;
				R = sqrt(Tax * Tax + Tay * Tay + Taz * Taz);

				Theta_Hx = Hax_x[index_ij] * Hax_x[index_ij] + Hax_y[index_ij] * Hax_y[index_ij] + Hax_z[index_ij] * Hax_z[index_ij];
				Theta_Hy = Hay_x[index_ij] * Hay_x[index_ij] + Hay_y[index_ij] * Hay_y[index_ij] + Hay_z[index_ij] * Hay_z[index_ij];
				Theta_Za = Za_x[index_ij] * Za_x[index_ij] + Za_y[index_ij] * Za_y[index_ij] + Za_z[index_ij] * Za_z[index_ij];
				E = sqrt((Theta_Hx + Theta_Hy + Theta_Za) / 2.0);
				outdata[index_ij] = sqrt(fabs(2 * E * E - R * R));
			}
		}

		////释放内存
		processout_Hax.Delete();
		processout_Hay.Delete();
		processout_Za.Delete();

		processout_Za_z.Delete();
		processout_Za_x.Delete();
		processout_Za_y.Delete();
		processout_Hax_x.Delete();
		processout_Hax_y.Delete();
		processout_Hax_z.Delete();
		processout_Hay_x.Delete();
		processout_Hay_y.Delete();
		processout_Hay_z.Delete();

		return 1;
	}
	int MagComponentTransT2L(CGM_GRDDATA indata, MagneticComponentTransStruct* magdatainfo, CGM_GRDDATA* processdata)
{
	int number_x, number_y;
	indata.GetDimension(number_x, number_y);
	int N = number_x*number_y;
	double* originaldata = indata.GetData();
	double *bounds = indata.GetBounds();
	//首先计算三分量
	//1.1 计算hax
	CGM_GRDDATA processout_Hax;
	processout_Hax.Initialize(number_x, number_y, bounds);
	MagComponentTransT2X(indata, magdatainfo, &processout_Hax);
	//1.2 计算hay
	CGM_GRDDATA processout_Hay;
	processout_Hay.Initialize(number_x, number_y, bounds);
	MagComponentTransT2Y(indata, magdatainfo, &processout_Hay);
	//1.3 计算Za
	CGM_GRDDATA processout_Za;
	processout_Za.Initialize(number_x, number_y, bounds);
	MagComponentTransT2Z(indata, magdatainfo, &processout_Za);
	//1.4 计算三分量的三个方向的导数
	CGM_GRDDATA processout_Hax_x;//Hax的导数
	processout_Hax_x.Initialize(number_x, number_y, bounds);
	DerivativeX(processout_Hax, &processout_Hax_x);
	CGM_GRDDATA processout_Hax_y;
	processout_Hax_y.Initialize(number_x, number_y, bounds);
	DerivativeY(processout_Hax, &processout_Hax_y);
	CGM_GRDDATA processout_Hax_z;
	processout_Hax_z.Initialize(number_x, number_y, bounds);
	DerivativeVertical(processout_Hax, &processout_Hax_z);
	CGM_GRDDATA processout_Hay_x;	//Hay的导数
	processout_Hay_x.Initialize(number_x, number_y, bounds);
	DerivativeX(processout_Hay, &processout_Hay_x);
	CGM_GRDDATA processout_Hay_y;
	processout_Hay_y.Initialize(number_x, number_y, bounds);
	DerivativeY(processout_Hay, &processout_Hay_y);
	CGM_GRDDATA processout_Hay_z;
	processout_Hay_z.Initialize(number_x, number_y, bounds);
	DerivativeVertical(processout_Hay, &processout_Hay_z);
	CGM_GRDDATA processout_Za_x;	//Za的导数
	processout_Za_x.Initialize(number_x, number_y, bounds);
	DerivativeX(processout_Za, &processout_Za_x);
	CGM_GRDDATA processout_Za_y;	//Za的导数
	processout_Za_y.Initialize(number_x, number_y, bounds);
	DerivativeY(processout_Za, &processout_Za_y);
	CGM_GRDDATA processout_Za_z;	//Za的导数
	processout_Za_z.Initialize(number_x, number_y, bounds);
	DerivativeVertical(processout_Za, &processout_Za_z);

	//2.0计算模量
	double* outdata = processdata->GetData();
	double* Hax = processout_Hax.GetData();
	double* Hay = processout_Hay.GetData();
	double* Za = processout_Za.GetData();
	double* Za_z = processout_Za_z.GetData();
	double* Za_x = processout_Za_x.GetData();
	double* Za_y = processout_Za_y.GetData();
	double* Hax_x = processout_Hax_x.GetData();
	double* Hax_y = processout_Hax_y.GetData();
	double* Hax_z = processout_Hax_z.GetData();
	double* Hay_x = processout_Hay_x.GetData();
	double* Hay_y = processout_Hay_y.GetData();
	double* Hay_z = processout_Hay_z.GetData();
	double Theta_Hx, Theta_Hy, Theta_Za;
	double Ta, Tax, Tay, Taz;
	double E, R,Q;
	int index_ij = 0;
	for (int i = 0; i < number_y; i++)
	{
		for (int j = 0; j < number_x; j++)
		{
			index_ij = j + i*number_x;
			Ta = sqrt(Hax[index_ij] * Hax[index_ij] + Hay[index_ij] * Hay[index_ij] + Za[index_ij] * Za[index_ij]);
			Tax = (Hax[index_ij] * Hax_x[index_ij] + Hay[index_ij] * Hay_x[index_ij] + Za[index_ij] * Za_x[index_ij]) / Ta;
			Tay = (Hax[index_ij] * Hax_y[index_ij] + Hay[index_ij] * Hay_y[index_ij] + Za[index_ij] * Za_y[index_ij]) / Ta;
			Taz = (Hax[index_ij] * Hax_z[index_ij] + Hay[index_ij] * Hay_z[index_ij] + Za[index_ij] * Za_z[index_ij]) / Ta;
			R = sqrt(Tax * Tax + Tay * Tay + Taz * Taz);

			Theta_Hx = Hax_x[index_ij] * Hax_x[index_ij] + Hax_y[index_ij] * Hax_y[index_ij] + Hax_z[index_ij] * Hax_z[index_ij];
			Theta_Hy = Hay_x[index_ij] * Hay_x[index_ij] + Hay_y[index_ij] * Hay_y[index_ij] + Hay_z[index_ij] * Hay_z[index_ij];
			Theta_Za = Za_x[index_ij] * Za_x[index_ij] + Za_y[index_ij] * Za_y[index_ij] + Za_z[index_ij] * Za_z[index_ij];
			E = sqrt((Theta_Hx + Theta_Hy + Theta_Za) / 2.0);

			Q = sqrt(fabs(2 * E * E - R * R));
			outdata[index_ij] = Q * Q / Ta;
		}
	}

	////释放内存
	processout_Hax.Delete();
	processout_Hay.Delete();
	processout_Za.Delete();

	processout_Za_z.Delete();
	processout_Za_x.Delete();
	processout_Za_y.Delete();
	processout_Hax_x.Delete();
	processout_Hax_y.Delete();
	processout_Hax_z.Delete();
	processout_Hay_x.Delete();
	processout_Hay_y.Delete();
	processout_Hay_z.Delete();

	return 1;
}

	int EdgeDetective_ImprovedTiltDerivative(CGM_GRDDATA indata, CGM_GRDDATA* processdata)
	{
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double* originaldata = indata.GetData();
		double *bounds = indata.GetBounds();
		//首先计算三分量
	
		//1.4 计算三分量的三个方向的导数
		CGM_GRDDATA processout_x;//Hax的导数
		processout_x.Initialize(number_x, number_y, bounds);
		DerivativeX(indata, &processout_x);
		CGM_GRDDATA processout_y;
		processout_y.Initialize(number_x, number_y, bounds);
		DerivativeY(indata, &processout_y);
		CGM_GRDDATA processout_z;
		processout_z.Initialize(number_x, number_y, bounds);
		DerivativeVertical(indata, &processout_z);
		
		//2.0计算模量
		double* indata_x = processout_x.GetData();
		double* indata_y = processout_y.GetData();
		double* indata_z = processout_z.GetData();
		double* outdata = processdata->GetData();
		double* pIndata = indata.GetData();
		double GradientMoudle,EdgeResult;						//梯度模
		int index_ij = 0;
		for (int i = 0; i < number_y; i++)
		{
			for (int j = 0; j < number_x; j++)
			{
				index_ij = j + i*number_x;
				GradientMoudle = sqrt(indata_x[index_ij] * indata_x[index_ij] + indata_y[index_ij] * indata_y[index_ij] + indata_z[index_ij] * indata_z[index_ij]);
				EdgeResult = atan2(indata_z[index_ij], GradientMoudle);

				outdata[index_ij] = (EdgeResult<0 ? 0 : EdgeResult);
			}
		}

		////释放内存
		processout_x.Delete();
		processout_y.Delete();
		processout_z.Delete();
		return 1;
	}

	int ExtenBoundary(CGM_GRDDATA indata, CGM_GRDDATA* processdata, int method )
	{
		//待扩边数据
		int number_x, number_y;
		indata.GetDimension(number_x, number_y);
		int N = number_x*number_y;
		double* originaldata = indata.GetData();
		double *bounds = indata.GetBounds();
		double dy = fabs(bounds[3] - bounds[2]) / (number_y - 1), dx = fabs(bounds[1] - bounds[0]) / (number_x - 1);
		//扩边后数据
		int number_x_ext, number_y_ext;
		processdata->GetDimension(number_x_ext, number_y_ext);
		int N_ext = number_x_ext*number_y_ext;
		double *extdata = processdata->GetData();
		double* bounds_ext = processdata->GetBounds(); 
		int num_ext_left = (number_x_ext - number_x) / 2;				//扩边后左边多出来的点数
		int num_ext_right = number_x_ext - number_x - num_ext_left;		//扩边后右边多出来的点数
		int num_ext_down = (number_y_ext - number_y) / 2;				//扩编后下面多出的点数
		int num_ext_up = number_y_ext - number_y - num_ext_down;		//扩边后上面多出的点数
		//重新计算坐标范围
		bounds_ext[0] = bounds[0] - num_ext_left*dx;
		bounds_ext[1] = bounds[1] + num_ext_right*dx;
		bounds_ext[2] = bounds[2] - num_ext_down*dy;
		bounds_ext[3] = bounds[3] + num_ext_up*dy;
		int index_ij_ext,index_ij;
		//扩边计算
		//1. 左边部分
		for (int i = 0; i < number_y; i++)
		{
			for (int j = 0; j < num_ext_left;j++)
			{
				index_ij_ext = num_ext_left-1-j + (i + num_ext_down)*number_x_ext;
				index_ij = i*number_x;
				extdata[index_ij_ext] = originaldata[index_ij] * cos(PI*(1 + j) / 2.0 / num_ext_left);
			}
		}
		//2. 中间部分
		for (int i = 0; i < number_y; i++)
		{
			for (int j = 0; j < number_x; j++)
			{
				index_ij = j + i*number_x;
				index_ij_ext = (i + num_ext_down)*number_x_ext + (j + num_ext_left);
				extdata[index_ij_ext] = originaldata[index_ij];
			}
		}
		//3. 右边部分
		for (int i = 0; i < number_y; i++)
		{
			for (int j = 0; j < num_ext_right; j++)
			{
				index_ij_ext = j+num_ext_left+number_x + (i+num_ext_down)*number_x_ext;
				index_ij = i*number_x+number_x-1;
				extdata[index_ij_ext] = originaldata[index_ij] * cos(PI*(j + 1) / 2.0 / num_ext_right);
			}
		}
		//4. 下边部分
		for (int i = 0; i < num_ext_down; i++)
		{
			for (int j = 0; j < number_x_ext; j++)
			{
				index_ij_ext = j + i*number_x_ext;
				index_ij = (num_ext_down)*number_x_ext+j;
				extdata[index_ij_ext] = extdata[index_ij] * cos(PI*(1+i)/2.0/num_ext_down);
			}
		}
		//5. 上边部分
		for (int i = 0; i < num_ext_up; i++)
		{
			for (int j = 0; j < number_x_ext; j++)
			{
				index_ij_ext = j + i*number_x_ext;
				index_ij = (num_ext_down+number_x)*number_x_ext + j;
				extdata[index_ij_ext] = extdata[index_ij] * cos(PI*(1 + i) / 2.0 / num_ext_up);
			}
		}
		return 1;
	}

	
}

