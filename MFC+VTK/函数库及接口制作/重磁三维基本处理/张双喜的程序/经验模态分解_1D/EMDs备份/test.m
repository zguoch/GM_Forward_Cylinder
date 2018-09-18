clear all
close all
clc

% t1=0:.001:2.5;
% t2=2.5:.001:5;
% t=0:.001:(5+.001);
% a=sin(20*pi*t1+0.5*sin(5*pi*t1));
% b=sin(10*pi*t2+.5*sin(5*pi*t2));
% c=0.5*cos(3*pi*t);
% x=[a(:,:) b(:,:)]+c;
%stop = [0.1,0.1,0.1];
%Ts=0.001;
%Fs=1/Ts;%采样频率
%t=0:Ts:1;

% N=1024;t=linspace(0,1,N);
% signal=sin(2*pi*50*t)+sin(2*pi*10*t);
% x=signal';

N=512;t=linspace(0,1,N);
%N=1024;n=0:N-1;t=Ts*n;
x=sin(2*pi*50*t)+sin(2*pi*10*t);
x=x';
%}
%x1=6*sin(2*pi*100*t)+5*cos(2*pi*400*t)+10*sin(2*pi*1000*t);
%x1=(1+0.3*sin(2*pi*15*ts))*cos(2*pi*50*ts+0.5*sin(2*pi*15*ts))+sin(2*pi*150*ts);
%x1=cos(tan(pi*t));
%x=awgn(x1,40,'measured');%y = awgn(x,SNR,SIGPOWER);x为信号，SNR为信噪比，若SIGPOWER为‘measured’
%则函数将在加入噪声之前测定信号强度
% x=sin(2*pi*90*t)+sin(2*pi*30*t);
% x1=sin(2*pi*90*t)+sin(2*pi*30*t);
%x=sin(2*pi*30*t)+1*randn(1, length(t)); %用于去噪
%{
y=zeros(1,length(t));
for i=1:500
    x(i)=2*sin(2*pi*240*t(i));
end
for i=501:1000
    x(i)=3*cos(2*pi*180*t(i));
end
%}
[imf,ort,nbits]=emd(x);
%原始信号的绘制==========================================================
figure;
plot(t,x);title('模拟信号');ylabel('幅值');xlabel('时间  t/s');
%IMF各分量================================================================
figure;
imfsize=size(imf);
i=imfsize(1);
for j=1:i
    subplot(i/2+1,2,j);plot(t,imf(j,:));
    if(j<i) Ylabel(sprintf('IMF%d', j));
    else Ylabel('Res');
    end
end

%个IMF各分量的合成与原始信号的比较=========================================
figure;
 newrdxn=x;
 for k=1:i
    newrdxn=newrdxn-imf(k,:);
    subplot(i/2+1,2,k);plot(t,newrdxn);title(sprintf('IMF%d的剩余信号', k));
    %k=k+1;
 end 
figure;
m=0;
for n=1:i
    m=m+imf(n,:);
end
%使用循环对各个IMF分量进行相加
m=m+newrdxn;%构造原始信号，IMF各分量和剩余信号
%m=imf(3,:)+imf(4,:)+imf(5,:)+imf(6,:)+newrdxn;%lk生成原始信号，可实现去噪
subplot(3,1,1);plot(t,x);%plot(t,x,'b',t,m,'r');
title('原始信号');%legend('原始信号','合成后的信号');
subplot(3,1,2);plot(t,m);title('合成后的信号');
diffence=m-x;
subplot(3,1,3);plot(t,diffence);title('分解合成后的信号与原始信号的误差');
%axis tight;
%绘制原始信号的频谱========================================================
figure;
[A,f,tt] = hhspectrum(imf(1:end-1,:));
%plot(fftshift(abs(fft(x,1024))));
plot(abs(fft(x)));
title('频谱图');ylabel('幅值');xlabel('时间');
axis([0,500,0,1000]);%设置坐标轴的取值范围
%axis tight;%把数值变化范围设置为刻度
% HHT变换==============================================================
freq_resol=300;
[A,f,tt] = hhspectrum(imf(1:end,:));%HHT谱画出除最后一条的IMF分量，hhspectrum（imf）画出所有分量
%im=toimage(A,f,tt);
%[A,f,tt] = hhspectrum(imf(1:end-1,:));
[im,ttt,ff] = toimage(A,f,tt,length(tt),freq_resol);
%disp_hhs(im,t);%绘制归一化频率
disp_hhs(im,[],N);%绘制真实频率，N为采样频率
set(gca,'xtick',[0:409.6:4096],'xticklabel',[0:0.1:1]);%设置x坐标，前后设置要一致
%hgsave('HH_Spectrum');%以名字'HH_Spectrum'保存图片
colorbar('yticklabel',{[0:2:20]})%若色棒出现负值的情况下，可根据其取值范围修改其从零到大发生变化
%colorbar;
%HHT边际谱============================================================
 %im=flipud(im);%是一个使矩阵上下翻转
 fs=1;
for kkk=1:size(im,1)
    aaa(kkk)=sum(im(kkk,:))*1/fs;
end
figure;
%f=(0:N-3)/N*(fs/2);
plot(N*ff(1,:),aaa);%归一化频率乘于采样频率就是真实频率,N为采样频率
%plot(ff(1,:),aaa);%求取归一化频率
axis tight;
title('边际谱');%显示每个频率的总振幅
xlabel('频率');ylabel('幅值');
axis([0,500,0,1000]);%设置坐标轴的取值范围
%横坐标是归一化的频率
%%%边际谱表示信号中某一频率在各个时刻的幅值值之和，更能真实的反应频率在信号中是否存在
%将所有时刻某一瞬时频率的能量（幅值）加起来就是信号中该频率的总能量（总幅值），即边际谱线的高度
%%%从边际谱上可以看出信号频率分布的真实情况
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

%%%%===========基于EMD(EEMD)的小波阈值去噪=============================================
%方法原理： (1)对原始数据序列x(t)进行EEMD分解，设分解阶数为x;
%          (2)由于第一阶IMF几乎是噪声，可直接将其去除
%          (3)对IMF2-IMFk进行阈值降噪处理
%          (4)对阈值量化后的各IMF分量和剩余分量进行重构，得到降噪后的信号
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
%%%函数ddencmp的调用==============(效果较好)======
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
plot(t,x1);title('不带噪声原始信号');
subplot(412);
plot(t,x);title('带噪声原始信号');axis tight;
subplot(413);
plot(t,q);title('滤波合成后的信号');axis tight;
subplot(414);
%err=q-x;
plot(t,x,'b',t,q,'r');legend('带噪声原始信号','滤波合成后的信号');
title('滤波效果比较');axis tight;
%}
%================利用小波包进行阈值滤波=====================================
%{
figure
%%%使用db2小波包对信号x进行s层分解
%%%%使用shannon熵
s=3;%可对其进行修改
wpt=wpdec(x,s,'db2');
%plot(wpt);
%%%计算结点的系数
for m=1:s
    for n=1:2^m
        cfs(m,n)=wpcoef(wpt,[m,n-1]);
        figure
        plot(cfs(m,n));
    end
end
%}
%{
%%%函数ddencmp的调用==============(效果较好)======
q=0;
for m=1:s
    for n=1:2^m
    xt = zeros(1,length(t));
    % 获取消噪的阀值
   [thr,sorh,keepapp] = ddencmp('den','wv',cfs(m,n));   %den表示进行去噪，wv表示选择小波
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
end
plot(q);
%}
%%%====观察去噪效果的好坏========================================
%{
m='rm';%也可以设置为“rv”
%x1为原始信号，x为重构信号，
value=SNR(x1,x,m)    %获得未处理的信号的信噪比
%%%函数ddencmp的调用，实现阈值去噪====================
    xt = zeros(1,length(t));
    % 获取消噪的阀值
   [thr,sorh,keepapp] = ddencmp('den','wv',x);   %den表示进行去噪，wv表示选择小波
   % 对信号进行消噪
   xt = wdencmp('gbl',x,'db4',3,thr,sorh,keepapp);%gbl表示选择相同的阈值,3表示小波分解的层数
%    %%%函数thselect的调用
%    %获取消噪的阀值
%    THR=thselect(nx,'rigrsure');%自适应阈值选择使用stein的无偏风险估计原理
%    %'heursure'使用启发式阈值选择;'sqtwolog'阈值等于sqrt(2*log(length(x)));'minimaxi'用极大极小原理选择阈值
%    % 对信号进行消噪
%    xd=wdencmp('gbl',nx,'db4',2,THR,'h');%s为软阈2值，h表示硬阈
value1=SNR(x1,xt,m)   %获得经过小波阈值处理后的信噪比
%%%====基于EMD的小波阈值滤波====================================
q=0;
for j=1:i
    xt = zeros(1,length(t));
    % 获取消噪的阀值
   [thr,sorh,keepapp] = ddencmp('den','wv',imf(j,:));   %den表示进行去噪，wv表示选择小波
   % 对信号进行消噪
   xt = wdencmp('gbl',imf(j,:),'db4',3,thr,sorh,keepapp);%gbl表示选择相同的阈值,3表示小波分解的层数 
%    %%%函数thselect的调用
%    %获取消噪的阀值
%    THR=thselect(nx,'rigrsure');%自适应阈值选择使用stein的无偏风险估计原理
%    %'heursure'使用启发式阈值选择;'sqtwolog'阈值等于sqrt(2*log(length(x)));'minimaxi'用极大极小原理选择阈值
%    % 对信号进行消噪
%    xd=wdencmp('gbl',nx,'db4',2,THR,'h');%s为软阈2值，h表示硬阈值
   q=q+xt;
end
%m='rm';%也可以设置为“rv”
value2=SNR(x1,q,m)   %获得经过基于EMD小波阈值处理后的信噪比
%%%当信噪比较高时，小波分解层数可设置为2,；当信噪比较低时，小波分解层数可设置大些，
%%%当信噪比较大的时候，小波去噪效果较好；当信噪比较小的时候，基于EMD小波阈值去噪处理效果好
%}
%{
figure
subplot(311)
plot(t,x1);title('带噪声信号');
xlabel('时间');ylabel('幅值');
subplot(312)
plot(t,xt);title('小波阈值去噪信号');
xlabel('时间');ylabel('幅值');
subplot(313)
plot(t,q);title('基于EMD小波阈值去噪信号');
xlabel('时间');ylabel('幅值');
%}
%%====================小波变换用于时频分析=================================
%{
%%interpreting the cwt coefficients
－、绘制原理

    1.需要用到的小波工具箱中的三个函数

    COEFS = cwt(S,SCALES,'wname') 
    说明：该函数能实现连续小波变换，其中S为输入信号，SCALES为尺度，wname为小波名称。
          COEFS为进行连续小波变换后返回的系数矩阵
        
    FREQ = centfrq('wname')
    说明：该函数能求出以wname命名的母小波的中心频率。

    F = scal2frq(A,'wname',DELTA) 
    说明：该函数能将尺度转换为实际频率，其中A为尺度，wname为小波名称，DELTA为采样周期。

    注：这三个函数还有其它格式，具体可参阅matlab的帮助文档。

    2.尺度与频率之间的关系

    设a为尺度，fs为采样频率，Fc为小波中心频率，则a对应的实际频率Fa为
                      
                      Fa＝Fc×fs/a                                     (1)

显然，为使小波尺度图的频率范围为(0,fs/2)，尺度范围应为(2*Fc,inf),其中inf表示为无穷大。在实际应用中，只需取尺度足够大即可。
   
    3.尺度序列的确定 

    由式(1)可以看出，为使转换后的频率序列是一等差序列，尺度序列必须取为以下形式：
    
                 c/totalscal,...,c/(totalscal-1),c/4,c/2,c        (2)  

其中，totalscal是对信号进行小波变换时所用尺度序列的长度(通常需要预先设定好)，c为一常数。

   下面讲讲c的求法。
     
    根据式(1)容易看出，尺度c/totalscal所对应的实际频率应为fs/2，于是可得

                      c=2×Fc/totalscal                               (3)

将式(3)代入式(2)便得到了所需的尺度序列。
    
    4.时频图的绘制

    确定了小波基和尺度后，就可以用cwt求小波系数coefs（系数是复数时要取模），然后用scal2frq将尺度序列转换为实际频率序列f，
最后结合时间序列t，用imagesc(t,f,abs(coefs))便能画出小波时频图。

    注意：直接将尺度序列取为等差序列，例如1:1:64，将只能得到正确的尺度－时间－小波系数图，而无法将其转换为频率－时间－小波系数图。这是因为此时的频率间隔不为常数。
此时，可通过查表的方法将尺度转化为频率或直接修改尺度轴标注。同理，利用本帖所介绍的方法只能得到频率－时间－小波系数图，不能得到正确的尺度－时间－小波系数图

说明：(1)应用时只须改变wavename和totalscal两个参数即可。
              (2)在这个例子中，最好选用复的morlet小波，其它小波的分析效果不好，而且morlet小波的带宽参数和中心频率取得越大，时频图上反映的时频聚集性越好。

%}
%{
% clear;
% clc;
fs=4096;  %采样频率，可根据需要设置采样频率，一般为2的整数次幂
% f1=100;                         
% f2=200;
% t=0:1/fs:1; %读取数据需将其注释
%{
R=xlsread('Gravity');%读取数据
t=R(:,1);
s=R(:,2);
%}
%s=10*sin(2*pi*f1*t)+15*sin(2*pi*f2*t);  %两个不同频率正弦信号合成的仿真信号
%s=sin(2*pi*f1*t)/2;  
%%%%分段函数的时频分析================================================
figure
s=x1;
%%%小波时频图绘制=======================================================
wavename='cmor3-3';   %小波名称 复小波cgau8，shan1-1.5,bior4.4,coif5,sym6,db4
%wavename='bior3.5';
totalscal=256;                    %尺度序列的长度，即scal的长度,可根据实际需要对其进行更改
wcf=centfrq(wavename);            %小波的中心频率
cparam=2*wcf*totalscal;           %为得到合适的尺度所求出的参数
a=totalscal:-1:1; 
scal=cparam./a;                   %得到各个尺度，以使转换得到频率序列为等差序列
coefs=cwt(s,scal,wavename);       %得到小波系数
f=scal2frq(scal,wavename,1/fs);   %将尺度转换为实际频率
%尺度的倒数对应于频率：尺度大，频率低，频率分辨率高，时间分辨率低；尺度小，频率高，频率分辨率低，时间分辨率高
imagesc(t,f,abs(coefs));          %绘制色谱图
axis xy;colormap jet;
colorbar;
%colormap(flipud(gray));% 黑白显示,打印比较清晰
xlabel('时间 t/s');
ylabel('频率 f/Hz');%频率轴的最大值为采样频率的1/2
title('小波变换时频图');
figure;
mesh(t,f,abs(coefs));colormap jet;
xlabel('时间 t/s');
ylabel('频率 f/Hz');
zlabel('幅值');
title('小波变换时频图');
axis tight;
%}
%==========EMD用于二维重力异常的分离======================================
%{
N=1500;
D=600;%大圆柱体的埋深
%D1=1000;%更深大圆柱体的埋深
d=150;%小圆柱体的埋深
%R1=80;%更深大圆柱体的半径
R=30;%大圆柱体的半径
r=15;%小圆柱体的半径
X=1;%剩余密度
z1=X*pi*(R.^2);%单位长度大圆柱体的剩余质量
z2=X*pi*(r.^2);%单位长度小圆柱体的剩余质量
%z3=X*pi*(R1.^2);%单位长度更深大圆柱体的剩余质量
G=6.67;%万有引力常量
x=zeros(1,2*N);
for i=1:(2*N)
   x(i)=i-1-N;
end
M=zeros(1,2*N);
F=zeros(1,2*N);
K=zeros(1,2*N);
%Q=zeros(1,2*N);
for j=1:(2*N)
    g0=2*G*z2*r/((x(j)+1200).^2+d.^2);%从左向右依次表示
    g1=2*G*z1*R/((x(j)+800).^2+D.^2);
    g2=2*G*z2*r/((x(j)+500).^2+d.^2);
    
    g3=2*G*z2*r/((x(j)-500).^2+d.^2);
    g4=2*G*z1*R/((x(j)-800).^2+D.^2);
    g5=2*G*z2*r/((x(j)-1200).^2+d.^2);
   %g6=2*G*z3*R1/((x(j)).^2+D1.^2);
    M(j)=g0+g1+g2+g3+g4+g5;  %总异常
    F(j)=g0+g2+g3+g5;
    K(j)=g1+g4;
    %Q(j)=g6;
end
figure
curve=plot(x,M,'r',x,F,'g',x,K,'k');
%curve=plot(x,M,'r',x,F,'g',x,K,'k',x,Q,'c');
ylabel('重力值  Δg/g.u');xlabel('坐标值  X/m');
title('模型重力异常');
%legend('模型叠加异常','模型局部异常');
%set(lsline,'LineWidth',10);
set(curve(1),'linewidth',3);%设置曲线1的粗细
set(curve(2),'linewidth',1.5);%设置曲线1的粗细
set(curve(3),'linewidth',1.5);%设置曲线1的粗细
%set(curve(4),'linewidth',1.5);%设置曲线1的粗细
 %legend('模型叠加异常','模型局部异常','模型区域异常');
%legend('模型叠加异常','模型局部异常','模型区域异常');
% rdxn=zeros(1,2*N);%获取行向量
% for i=1:(2*N)
%     rdxn(i)=M(i);
% end
[imf,ort,nbits]=emd(M);
%IMF各分量================================================================
figure;
imfsize=size(imf);
i=imfsize(1);
for j=1:i
    subplot(i,1,j);plot(x,imf(j,:));
    if(j<i) Ylabel(sprintf('IMF%d', j));
    else Ylabel('Res');
    end
end
%对IMF进行重构，观察重构误差========================
w=zeros(1,2*N);
w(1)=zeros(1,2*N);
for j=1:i
   w(j+1)= w(j)+imf(j,:)
end
diffence=w-M;
plot(x,diffence);title('重力重构误差');
ylabel('重力重构误差值  Δg/g.u');xlabel('坐标值  X/m');
%}
%{
%%%====
%用db1小波函数对信号进行多尺度小波分解
[C,L]=wavedec(M,3,'sym1');
%提取尺度1的低频系数
cA1=appcoef(C,L,'sym1',1);
%提取尺度2的低频系数
cA2=appcoef(C,L,'sym1',3);
%提取尺度1的高频系数
cD1=detcoef(C,L,1);
%提取尺度2的高频系数
cD2=detcoef(C,L,3);
figure
subplot(411)
plot(cA1);
subplot(412)
plot(cD1);
subplot(413)
plot(cA2);
subplot(414)
plot(cD2);
%%%===================径向平均对数功率谱=========================================
%}
%{

dt = 1;%时间采样间隔dt，ms
fn = 1000/(2*dt);%奈奎斯特频率，hz
fmin=0;%最小频率从0Hz开始
fmax=fn %最大频率为奈奎斯特频率

%-----利用不同频率的Ricker子波制作一个合成记录--------------------------------
fp1=150;
fp2=90;
fp3=60;
fp4=30;
wave_t=120; %子波的时间长度，ms
w_t=-wave_t/(2*1000):dt/1000:wave_t/(2*1000);
rick1=(1-2*(pi*fp1*(w_t)).^2).*exp(-(pi*fp1*(w_t)).^2);
w_t=-wave_t/(2*1000):dt/1000:wave_t/(2*1000);
rick2=-1*(1-2*(pi*fp2*(w_t)).^2).*exp(-(pi*fp2*(w_t)).^2);
w_t=-wave_t/(2*1000):dt/1000:wave_t/(2*1000);
rick3=(1-2*(pi*fp3*(w_t)).^2).*exp(-(pi*fp3*(w_t)).^2);
w_t=-wave_t/(2*1000):dt/1000:wave_t/(2*1000);
rick4=-1*(1-2*(pi*fp4*(w_t)).^2).*exp(-(pi*fp4*(w_t)).^2);
syn1=[rick1 rick2 rick3 rick4];
wave_t=100; %子波的时间长度，ms
w_t=-wave_t/(2*1000):dt/1000:wave_t/(2*1000);
rick1=(1-2*(pi*fp1*(w_t)).^2).*exp(-(pi*fp1*(w_t)).^2);
w_t=-wave_t/(2*1000):dt/1000:wave_t/(2*1000);
rick2=-1*(1-2*(pi*fp2*(w_t)).^2).*exp(-(pi*fp2*(w_t)).^2);
w_t=-wave_t/(2*1000):dt/1000:wave_t/(2*1000);
rick3=(1-2*(pi*fp3*(w_t)).^2).*exp(-(pi*fp3*(w_t)).^2);
w_t=-wave_t/(2*1000):dt/1000:wave_t/(2*1000);
rick4=-1*(1-2*(pi*fp4*(w_t)).^2).*exp(-(pi*fp4*(w_t)).^2);
syn2=[rick1 rick2 rick3 rick4 zeros(1,80)];
syn=syn1+syn2;%合成一个信号

N=1000;%采样频率
[imf,ort,nbits]=emd(syn);
w_t=0:dt:length(syn)-1;
%{
figure
plot(w_t,syn);title('模拟信号');xlabel('时间  t/ms');
ylabel('幅值');
figure
N=512;n=0:N-1;t=0.001*n;q=n*1000/N;
y=fft(syn,N);plot(q,abs(y));title('经FFT得到的模拟信号的频谱');
xlabel('频率  f/Hz');ylabel('幅值');
%}
%IMF各分量================================================================
%原始信号的绘制==========================================================

%IMF各分量================================================================
figure;
imfsize=size(imf);
i=imfsize(1);
%个IMF各分量的合成与原始信号的比较=========================================

% HHT变换==============================================================
figure
freq_resol=300;
[A,f,tt] = hhspectrum(imf(1:end-1,:));%HHT谱画出除最后一条的IMF分量，hhspectrum（imf）画出所有分量
%im=toimage(A,f,tt);
%[A,f,tt] = hhspectrum(imf(1:end-1,:));
[im,ttt,ff] = toimage(A,f,tt,length(tt),freq_resol);
%disp_hhs(im,t);%绘制归一化频率
disp_hhs(im,[],N);%绘制真实频率，N为采样频率
set(gca,'xtick',[0:50:length(syn)],'xticklabel',[0:50:length(syn)]);%设置x坐标，前后设置要一致
%hgsave('HH_Spectrum');%以名字'HH_Spectrum'保存图片
%colorbar('yticklabel',{[0:2:20]})%若色棒出现负值的情况下，可根据其取值范围修改其从零到大发生变化
colorbar;
%HHT边际谱============================================================
 %im=flipud(im);%是一个使矩阵上下翻转
 N=1000;%采样频率
 fs=1;
for kkk=1:size(im,1)
    aaa(kkk)=sum(im(kkk,:))*1/fs;
end
figure;
%f=(0:N-3)/N*(fs/2);
plot(N*ff(1,:),aaa);%归一化频率乘于采样频率就是真实频率,N为采样频率
%plot(ff(1,:),aaa);%求取归一化频率
axis tight;
title('边际谱');%显示每个频率的总振幅
xlabel('频率  f/Hz');ylabel('幅值');
%}
