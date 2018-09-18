#include"FFTN.h"
void FFT2d(vector<vector<double>> invector_r, fftw_complex* out_c)
{
	//1. 申明数据输入输出数组
	fftw_complex *in;
	//2. 申明fft执行变量
	fftw_plan p;//,q;
	//3. 给出数据个数
	int rows = invector_r.size();
	int cols = invector_r[0].size();
	int N = rows*cols;
	//4. 根据数据个数开辟动态数组空间
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
	//out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
	//5. 给输入数据数组赋值
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			in[j + i*cols][0] = invector_r[i][j];
			in[j + i*cols][1] = 0.0;
		}
		//printf("%6.2f \n", in[i][0]);
	}
	//6. 调用fft执行函数
	p = fftw_plan_dft_2d(rows, cols, in, out_c, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p); /* repeat as needed*/

	//8. 最后销毁fft相关的变量
	fftw_destroy_plan(p);
	//fftw_destroy_plan(q);
	fftw_free(in);
}
void FFT2d(int rows, int cols, double* realin, fftw_complex* out_c)
{
	//1. 申明数据输入输出数组
	fftw_complex *in;
	//2. 申明fft执行变量
	fftw_plan p;//,q;
	//3. 给出数据个数
	int N = rows*cols;
	//4. 根据数据个数开辟动态数组空间
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
	//out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
	//5. 给输入数据数组赋值
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			in[j + i*cols][0] = realin[i*cols+j];
			in[j + i*cols][1] = 0.0;
		}
		//printf("%6.2f \n", in[i][0]);
	}
	//6. 调用fft执行函数
	p = fftw_plan_dft_2d(rows, cols, in, out_c, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p); /* repeat as needed*/

	//8. 最后销毁fft相关的变量
	fftw_destroy_plan(p);
	//fftw_destroy_plan(q);
	fftw_free(in);
}
void FFT1d(vector<double> invector_r, fftw_complex* out_c)
{
	//1. 申明数据输入输出数组
	fftw_complex *in;
	//2. 申明fft执行变量
	fftw_plan p;//,q;
	//3. 给出数据个数
	int N = invector_r.size();
	//4. 根据数据个数开辟动态数组空间
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
	//out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
	//5. 给输入数据数组赋值
	for (int i = 0; i < N; i++)
	{
		in[i][0] = invector_r[i];
		in[i][1] = 0.0;
		//printf("%6.2f \n", in[i][0]);
	}
	//6. 调用fft执行函数
	p = fftw_plan_dft_1d(N, in, out_c, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p); /* repeat as needed*/

	//8. 最后销毁fft相关的变量
	fftw_destroy_plan(p);
	//fftw_destroy_plan(q);
	fftw_free(in);
}
void FFT1d(int N, double* in, fftw_complex* out_c)
{
	//2. 申明fft执行变量
	fftw_plan p;//,q;
	//6. 调用fft执行函数
	p = fftw_plan_dft_r2c_1d(N, in, out_c, FFTW_ESTIMATE);
	fftw_execute(p); /* repeat as needed*/
	//8. 最后销毁fft相关的变量
	fftw_destroy_plan(p);
	fftw_free(in);
}
void IFFT1d(int N, fftw_complex* fftout, double* out)
{
	fftw_plan p;//,q;
	p = fftw_plan_dft_c2r_1d(N, fftout, out, FFTW_ESTIMATE);
	fftw_execute(p);
	//释放指针
	fftw_destroy_plan(p);
}
void IFFT1d(int N, fftw_complex* fftout, vector<double>& reoutvector)
{
	reoutvector.clear();

	fftw_complex *reout;
	fftw_plan p;//,q;
	reout = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
	p = fftw_plan_dft_1d(N, fftout, reout, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	//赋值结果
	for (int i = 0; i < N; i++)
	{
		reoutvector.push_back(reout[i][0] / N);
	}
	//释放指针
	fftw_destroy_plan(p);
	fftw_free(reout);
}
void IFFT2d(int rows, int cols, fftw_complex* fftout, vector<double>& reoutvector)
{
	reoutvector.clear();

	int N = rows*cols;
	fftw_complex *reout;
	fftw_plan p;//,q;
	reout = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
	p = fftw_plan_dft_2d(rows, cols, fftout, reout, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	//赋值结果
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			reoutvector.push_back(reout[j + i*cols][0] / N);
		}
	}

	//释放指针
	fftw_destroy_plan(p);
	fftw_free(reout);
}
void IFFT2d(int rows, int cols, fftw_complex* fftout, double* out)
{
	int N = rows*cols;
	fftw_complex *reout;
	fftw_plan p;//,q;
	reout = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
	p = fftw_plan_dft_2d(rows, cols, fftout, reout, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	//赋值结果
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			out[i*cols+j]=(reout[j + i*cols][0] / N);
		}
	}

	//释放指针
	fftw_destroy_plan(p);
	fftw_free(reout);
}
void GetSpectrum1d(int N, double dx, fftw_complex* fftout, vector<double>& fvector, vector<double>& SpectrumVector)
{
	fvector.clear();
	SpectrumVector.clear();

	double FS = 1.0 / dx;
	double df = FS / (N - 1);
	int Nhalf = (int)floor(N / 2.0);
	double f = 0, spectrum = 0;
	for (int i = 0; i < Nhalf; i++)
	{
		f = i*df;
		spectrum = sqrt(pow(fftout[i][0], 2.0) + pow(fftout[i][1], 2.0));
		fvector.push_back(f);
		SpectrumVector.push_back(spectrum);
	}
	return;
}