
clear;
clc;

%计算信号的信噪比
% load LorenzSNR1.mat;
% [DATAfile DATApath]=uigetfile('*.txt','输入信号');
% FILENAME=[DATApath,DATAfile];
% x=load(FILENAME);
% sig1=x';
% sig2=Y(1:length(sig1));
% % load Lorenz.mat;
% % sig1=x;
% % sig2=Y;
%{
f=30;
fs=1000;
number=200;
t=-number/2+1:number/2;
a=(1-2*(pi*f*t/fs).^2).*exp(-(pi*f*t/fs).^2);
s=awgn(a,10,'measured'); 
%}
N=4096;t=linspace(0,1,N);
%N=1024;n=0:N-1;t=Ts*n;
x=3*sin(2*pi*90*t)+2*cos(2*pi*180*t);
x1=3*sin(2*pi*90*t)+2*cos(2*pi*180*t)+0.1*randn(1,length(t));
sig1=x;   % sig1 原始信号
sig2=x1;   % sig2 重构信号
m='rm';
value=SNR(sig1,sig2,m)
%%%函数ddencmp的调用，实现阈值去噪====================
    xt = zeros(1,length(t));
    % 获取消噪的阀值
   [thr,sorh,keepapp] = ddencmp('den','wv',x1);   %den表示进行去噪，wv表示选择小波
   % 对信号进行消噪
   xt = wdencmp('gbl',x1,'db4',2,thr,sorh,keepapp);%gbl表示选择相同的阈值,3表示小波分解的层数
   
%    %%%函数thselect的调用
%    %获取消噪的阀值
%    THR=thselect(nx,'rigrsure');%自适应阈值选择使用stein的无偏风险估计原理
%    %'heursure'使用启发式阈值选择;'sqtwolog'阈值等于sqrt(2*log(length(x)));'minimaxi'用极大极小原理选择阈值
%    % 对信号进行消噪
%    xd=wdencmp('gbl',nx,'db4',2,THR,'h');%s为软阈2值，h表示硬阈
value1=SNR(sig1,xt,m)



