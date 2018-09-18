bool UpWardContinuation(double* OriginalData,double* TransData,const int datanum,const double dx,const double rph);
bool UpWardContinuation_Dct(double* OriginalData,double* TransData,const int datanum,const double dx,const double rph);
bool NormalFullGradient(double* OriginalData,double* TransData,const int datanum,const double dx,const double TRP,const double rph,int NormalizeType);
bool NormalFullGradient(double* OriginalData,double* TransData,const int datanum,const double dx,const double TRP,const int hNum,const double dh,const double h1,int NormalizeType);
bool NormalFullGradient_TSI(double* OriginalData,double* TransData,const int datanum,const double dx,const int N,const int L,const int K,const double rph);
bool NormalFullGradient_TSI(double* OriginalData,double* TransData,const int datanum,const double dx,const int N,const int L,const int K,const int hNum,const double dh,const double h1);
bool DownWardContinuation_Tik(double* OriginalData,double* TransData,const int datanum,const double dx,const double rph,double TRP);
bool Derivative_Z(double* OriginalData,double* TransData,const int datanum,const double dx);
bool Derivative_X(double* OriginalData,double* TransData,const int datanum,const double dx,int Order);
bool Derivative_X_space(double* OriginalData,double* TransData,const int datanum,const double dx,int Order);

bool SourceDepthInversion(double* OriginalData,double* TransData,const int datanum,const double dx,const double P,const double rph);
bool SourceDepthInversion(double* OriginalData,double* TransData,const int datanum,const double dx,const double P,const int hNum,const double dh,const double h1);
bool AmplitudeSpectrum(double* OriginalData,const int datanum,const double dx,char* outfilename);
bool Derivative_Z(double* OriginalData,double* TransData,const int datanum,const double dx,const int hNum,const double dh,const double h1);
bool NormalEAS(double* OriginalData,double* TransData,const int datanum,const double dx,const double TRP,const double rph,int NormalizeType);
bool NormalEAS(double* OriginalData,double* TransData,const int datanum,const double dx,const double TRP,const int hNum,const double dh,const double h1,int NormalizeType);
bool Derivative_Z_DHT(double* OriginalData,double* TransData,const int datanum,const double dx);
bool Derivative_X_DHT(double* OriginalData,double* TransData,const int datanum,const double dx);
bool NormalLocalWaveNumber(double* OriginalData,double* TransData,const int datanum,const double dx,const double TRP,const double rph,int NormalizeType);
bool Separation_EntropyFilter(double* OriginalData,double* TransData,const int datanum,const double dx,int WinSize);
bool NormalizeDownwardConti(double* OriginalData,double* TransData,const int datanum,const double dx,const double TRP,const double rph,int NormalizeType);
