
clear;
clc;
tic;
%{
[x,FS]=wavread('1.wav');
sigLength=length(x);
%h=window('hamming',64);
%[B,f,t]=specgram(x,256,22050,h,60);
t=(0:sigLength-1)/FS;
ts=t;
rdxn=x;
%}
%
%data=load('gphone_velocity(2010.5.24-2010.5.31).dat');
%data=load('Gravity.txt');%加载数据
%data=load('L4Vz.dat');%加载数据,时间信号

data=load('Gravity.txt');
ts=data(:,1);%读取第一列数据
rdxn=data(:,2);

%rdxn=data(:,4);
% %FS=1;%数据的采样频率,若隔0.02s取一个点，则采样频率为1/0.002;
% N=4096;ts=linspace(0,1,N);
% rdxn=(1+0.3*sin(2*pi*15*ts)).*cos(2*pi*50*ts+0.5*sin(2*pi*15*ts))+sin(2*pi*150*ts);

% N=4096;ts=linspace(0,1,N);
% rdxn=0.5*sin(2*pi*400*ts)+0.5*cos(2*pi*100*ts)+0.6*sin(2*pi*500*ts);



%load ecg1410.txt;
%rdxn=ecg1410;
% load eeg0410.mat;
% rdxn=eeg0410;
N=length(rdxn);%数据的长度
%n=1:N;
fs=N/(ts(N)-ts(1));
%t=n/fs;
t_inst=1:length(rdxn); %time vector (drugoj uzhe)
FS=fs;



%绘制原始数据的变化曲线=============================================
figure;
%whitebg(gcf,'black');%%%将图形的绘制窗口的背景色设置为黑色
plot(ts,rdxn);
%title('gphone-data');
title('原始信号');xlabel('时间  t/ms');ylabel('幅值');
%EMD变换
[imf,ort,nbits]=emd(rdxn);
imfsize=size(imf);
i=imfsize(1);

%{
%FFT显示频谱==============(用来显示非平稳信号的频谱，意义不大）==========================================
%
N1=4096;%2的整数次幂,采样点数，数据的采样点数=============处理数据时，根据实际情况需要修改项
n=0:N1-1;%t=1/4096
fs=4096;%采样频率,与数据采样频率一致==========需要修改项
X=fft(rdxn,N1);
% figure;
% plot(abs(X));
figure;
%whitebg(gcf,'black');%%%将图形的绘制窗口的背景色设置为黑色
plot(n*fs/N1,abs(X));%axis([0,1000,0,15000]);
%axis tight;
title('原始信号的频谱图');
%}

%各IMF分量及其残余量====================================================
figure;
%whitebg(gcf,'black');%%%将图形的绘制窗口的背景色设置为黑色
%title('EMD分解');
newrdxn=rdxn';
for j=1:i
    subplot(i/2+1,2,j);plot(ts,imf(j,:));
     if(j<i) Ylabel(sprintf('IMF%d', j));
     else Ylabel('Res');
     end
   %  newrdxn=newrdxn-imf(j,:);
end
% subplot(i/2+1,2,j+1);
% plot(newrdxn);Ylabel('残余量');%最后一次分解的残余量
 %每分解一次后的残余量==================================================
 figure;
% title('残余量');
 newrdxn=rdxn';
 for k=1:i
     newrdxn=newrdxn-imf(k,:);
     subplot(i/2+1,2,k);plot(ts,newrdxn);
     k=k+1;
 end
%{
 %IMF-Fourier 谱===========================================================
  figure;
  %whitebg(gcf,'black');%%%将图形的绘制窗口的背景色设置为黑色
% M = length(imf);
% N = length(x);
% c = linspace(0,(N-1)*Ts,N);   
N1=N; 
for v=1:i
    t=fft(imf(v,:),N1);
    subplot(i/2,2,v),%i/2表示图形的显示形式，可根据需要进行设置
    plot(n*fs/N1,abs(t));axis tight;
    if(v<i)Ylabel(sprintf('IMF%d', v));%title(sprintf('IMF%d-Fourier 谱', v));
    else Ylabel('Res');%title('Res-Fourier 谱');
    end
    %ylabel('IMF-Fourier 谱 ')
    %v=v+1;
end
%}  
%HHT变换===============================================================
freq_resol=300;
[A,f,tt] = hhspectrum(imf(1:end-1,:));
[im,ttt,ff] = toimage(A,f,tt,length(tt),freq_resol);%对每个IMF信号合成求取Hilbert谱，im对应的振幅值，ff:每个网格对应的中心频率，这里横轴是时间，纵轴为频率
%disp_hhs(im,ttt);
disp_hhs(im,[],FS);%FS为数据采样频率
%set(gca,'xtick',[0:441.3108242:length(rdxn)],'xticklabel',[0:0.02:0.1007]);%设置x坐标，前后设置要一致
%length(rdxn)为数据长度，length(rdxn)/250.2与2.388/0.5相等
colorbar('yticklabel',{[0:2:20]})%若色棒出现负值的情况下，可根据其取值范围修改其从零到大发生变化
%colorbar;
hgsave('HH_Spectrum');
%colormap(flipud(gray)); % 黑白显示,打印比较清晰
% axisx=max(ttt);
% axisy=max(ff);
% WriteGrd(im,'time_frqence.grd',axisx,axisy);

%HHT边际谱===========================================================
%%%边际谱从统计意义上表征了整组数据每个频率点的累积幅值分布，边际谱可以处理非平稳信号，如果
%%%信号中存在某一频率的能量出现，就表示一定有该频率的振动波出现，也就是边际谱能比较准确地反
%%%映信号的实际频率成分
fs=1;
%fs=FS;
 %im=flipud(im);%是一个使矩阵上下翻转
for kkk=1:size(im,1)
    aaa(kkk)=sum(im(kkk,:))*1/fs;
end
figure;
%whitebg(gcf,'black');%%%将图形的绘制窗口的背景色设置为黑色
%f=(0:N-3)/N*(fs/2);
plot(FS*ff(1,:),aaa);axis tight;
title('边际谱');%归一化频率乘于采样频率就是真实频率,FS为采样频率
xlabel('频率  f/Hz');ylabel('幅值');
%%==============将分离的数据进行存入文件==============================
%{
%x=0:pi/100000:2*pi;
%N=length(x);
g1=zeros(1,N);%N为信号的长度
g2=zeros(1,N);
g3=zeros(1,N);
g4=zeros(1,N);
diffence=zeros(1,N);
fid1=fopen('高频信号12.dat','wt');
fid2=fopen('低频信号12.dat','wt');
%fprintf(fid,'%s    %s\n','横坐标','坐标值');%用于开始标注数据信息
Q=0;

for j=1:i
    %g(i)=7*sin(x(i))+1;
    if(i<=4)
        g1(j)=Q+imf(j,:);
        fprintf(fid1,'%f %f\n',ts,g1(j));
    end
    if(i>4)
        g2(j)=Q+imf(j,:);
        fprintf(fid2,'%f %f\n',ts,g2(j));
    end
end

for i=1:2
  g1=g1+imf(i,:);%+imf(2,:)+imf(3,:)+imf(4,:);
end
for i=3:4
   g2=g2+imf(i,:);%+imf(6,:)+imf(7,:)+imf(8,:)+;
end
   g3=imf(5,:)+imf(6,:);
   g4=imf(7,:)+imf(8,:);
 fprintf(fid1,'%f  %f \n',ts,g1);
 fprintf(fid2,'%f  %f \n',ts,g2);
%plot(x,g,'g+');
fclose(fid1);
fclose(fid2);
%===================================
%}


%{
g1=zeros(1,N);%N为信号的长度
g2=zeros(1,N);
g3=zeros(1,N);
g4=zeros(1,N);
%g5=zeros(1,N);
g1=imf(1,:)+imf(2,:)+imf(3,:)+imf(4,:);
g2=imf(5,:)+imf(6,:)+imf(7,:)+imf(8,:);
g3=imf(9,:)+imf(10,:)+imf(11,:)+imf(12,:);
g4=imf(13,:)+imf(14,:)+imf(15,:)+imf(16,:)+imf(17,:)+imf(18,:);
%g5=
rdxn=rdxn';
figure
subplot(511)
plot(ts,rdxn);
title('原始信号');
subplot(512)
plot(ts,g1);
title('高频信号')
subplot(513)
plot(ts,g2);
title('次高频信号');
subplot(514)
plot(ts,g3);
title('低频信号');
subplot(515)
plot(ts,g4);
title('次低频信息')

figure
g1=g1';
g2=g2';
g3=g3';
g4=g4'
%g5=g5';
subplot(311)
plot(ts,rdxn);
title('原始信号');
subplot(312)
plot(ts,g1+g2+g3+g4);
title('分解重构后的信号');
subplot(313)
rdxn=rdxn';
diffence=g1+g2+g3+g4-rdxn;
plot(ts,diffence);
title('信号分解重构后的误差');

%}







%plot(f,aaa);

%D=logspace(min(log(ff)),max(log(ff)),size(ff,2));
%loglog(D,aaa);
% [S,freq]=hspec(imf,1000);
%  figure;
%  w1=size(A,2);
%  w2=size(f,2);
  
%[AA,F] = mhs(A,f,0,0.5,1000);

%D=logspace(min(log(F)),max(log(F)),size(F,2));
%loglog(F,AA);
%plot(log(F),log(AA));
toc;
%%%=========基于EMD的小波阈值滤波============================================
%%%=====%%%%实现对一维信号自动消噪===========(效果较好)=============
%{
 q=0;
for j=1:i
    xt = zeros(1,length(t));
    % 将信号nx用小波函数sym5分解到第5层  可对小波函数进行修改如sym3
    % 用minimaxi：极大值极小值阀值选择对系数进行处理，消除噪声信号
    %sqtwolog 取通用阀值，heursure;启发式阀值选择；rigrsure用无偏似然估计
    lev = 2;%小波分解的层数
    xt= wden(imf(j,:),'minimaxi','s','mln',lev,'sym5');  %sym3为所用的小波函数;mln根据不同层的噪声估计来调整阈值，
    %sln 根据第一层的系数来进行噪声层的估计来调整阈值 one 表示不调整
    q=q+xt;
end
%}
%%%%============%%%对小波分解系数进行阈值处理，然后对处理后的小波系数进行重构达到去噪的目的===（效果较好）
%%%%%注意小波基的选择和分解尺度的选择，分解尺度越大，去噪效果越不好======================
%{
q=0;
for j=1:i
    xt = zeros(1,length(t));
   [c,l]=wavedec(imf(j,:),2,'db6');
   %设置尺度向量
   n=[1,2];
   %设置阀值向量
   p=[120,110];%n,p的长度必须相同
   %对高频系数进行阀值处理
   %nc=wthcoef('d',c,l,n,p,'s');
   nc=wthcoef('d',c,l,n);
   %nc=wthcoef('a',c,l);%高频成分
   %对修正后的小波分解结构进行重构，[nc,l]构成一个新的小波分解结构
   xt = waverec(nc,l,'db6');
   q=q+xt;
end
%}
%%%%============利用噪声标准偏差=========(效果较好)============
%{
q=0;
for j=1:i
    xt = zeros(1,length(t));
    %用小波函数 'db6'对信号进行3层分解
   [c,l] = wavedec(imf(j,:),3, 'db6');   
   %估计尺度1的噪声标准偏差
   sigma = wnoisest(c,l,1);      
   alpha = 2; %必须为一个大于1的实数
   %获取消噪过程中的阀值
   thr = wbmpen(c,l,sigma,alpha) ;   %c,l为进行去噪信号小波分解结构
   keepapp = 1;%当取值为1时，低频系数不进行阈值量化，反之低频系数要进行量化
   %对信号进行消噪
   xt= wdencmp('gbl',c,l,'db6',3,thr,'s',keepapp);
   q=q+xt;
end
%}
%%%函数ddencmp的调用==============(效果较好)====================================
%{
q=0;
for j=1:i
    xt = zeros(1,length(t));
    % 获取消噪的阀值
   [thr,sorh,keepapp] = ddencmp('den','wv',imf(j,:));   %den表示进行去噪，wv表示选择小波
   % 对信号进行消噪
   xt = wdencmp('gbl',imf(j,:),'db4',2,thr,sorh,keepapp);%gbl表示选择相同的阈值,3表示小波分解的层数
   
%    %%%函数thselect的调用
%    %获取消噪的阀值
%    THR=thselect(nx,'rigrsure');%自适应阈值选择使用stein的无偏风险估计原理
%    %'heursure'使用启发式阈值选择;'sqtwolog'阈值等于sqrt(2*log(length(x)));'minimaxi'用极大极小原理选择阈值
%    % 对信号进行消噪
%    xd=wdencmp('gbl',nx,'db4',2,THR,'h');%s为软阈2值，h表示硬阈值
   
   q=q+xt;
end
%}
%%%基于提升方案去噪=====去噪效果不好======
%{
q=0;
for j=1:i
    %xt = zeros(1,length(t));
    % 得到sym5小波的提升方案
    lshaar = liftwave('haar');
    % 添加ELS到提升方案中
    els = {'p',[-0.125 0.125],0};
    lsnew = addlift(lshaar,els);
    % 进行提升小波分解
   [cA1,cD1] = lwt(imf(j,:),lsnew);
   [cA2,cD2] = lwt(cA1,lsnew);%%%再对近似信号进行提升小波分解
   length = size(cA2,2);
   c = zeros(1,length*4);
   for i = 1:length;
    c(i) = cA2(i);
   end
   for i = length+1:2*length;
    c(i) = cD2(i-length);
   end;
   for i = length*2+1:4*length;
    c(i) = cD1(i-2*length);
   end;
   l(1) = length;
   l(2) = length;
   l(3) = length*2;
   l(4) = length*4;
   %估计尺度2的噪声标准偏差
   sigma = wnoisest(c,l,2);      
   alpha = 2;%用于处罚的调整参数，它必须是一个大于1的实数，一般取2
   %获取消噪过程中的全局阀值
   thr = wbmpen(c,l,sigma,alpha) ;   
   keepapp = 1;%%%低频系数不进行阈值量化，反之要进行阈值量化
   %对信号进行消噪
   xt= wdencmp('gbl',c,l,'db4',2,thr,'s',keepapp);
   %%%gbl全局阈值，lvb每一层用不同的阈值进行处理；2表示分解的层数
   %%%s表示软阈值，h表示硬阈值

   q=q+xt;
end
%}
%{
figure;
subplot(411);
plot(ts,rdxn);title('带噪声原始信号');axis tight;
subplot(412);
plot(ts,q);title('滤波合成后的信号');axis tight;
subplot(413);
%err=q-x;
q=q';%对维数进行转置，使其与rdxn的维数一直
plot(ts,rdxn,'b',ts,q,'r');legend('带噪声原始信号','滤波合成后的信号');
title('滤波效果比较');axis tight;
diffence=q-ts;%滤波前后的差值
subplot(414);
plot(ts,diffence,'g');title('滤波前后的差值');
axis tight;
%}
%%%==============利用小波变换进行时频分析===============================================
%{
%fs=4096;  %采样频率，可根据需要设置采样频率，一般为2的整数次幂
fs=FS;
%{
R=xlsread('Gravity');%读取数据
t=R(:,1);
s=R(:,2);
%}
%a=load('Gravity.txt');
%{
t=data(:,1);%读取第一列数据
s=rdxn;%读取第二列数据
%}
t=ts;
s=rdxn;
%读取第一行，第二行表示方法： t=a(1,:);s=a(2,:);
%======小波变换时频图的绘制========================================================
%wavename='cmor3-3';%复小波
wavename='db4';%带宽参数为4，中心频率为2 小波的带宽参数和中心频率取得越大，时频图上反映的时频聚集性越好
totalscal=256;                    %尺度序列的长度，即scal的长度
wcf=centfrq(wavename);            %小波的中心频率
cparam=2*wcf*totalscal;           %为得到合适的尺度所求出的参数
a=totalscal:-1:1;  
scal=cparam./a;                   %得到各个尺度，以使转换得到频率序列为等差序列
coefs=cwt(s,scal,wavename);       %得到小波系数
f=scal2frq(scal,wavename,1/fs);   %将尺度转换为频率
figure;
imagesc(t,f,abs(coefs));          %绘制色谱图
axis('xy');
colorbar;
%colormap(flipud(gray));% 黑白显示,打印比较清晰
xlabel('时间  t/ms');
ylabel('频率  f/Hz');%频率轴的最大值为采样频率的1/2
title('小波变换（WT）时频图');
%}
%==================================================================================