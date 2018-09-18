
clear ;
clc;

%模拟信号
N=2048;%采样点数
t=linspace(0,1,N);
sampFreq=N/(1-0);%采样频率
n=0:N-1;T=n*sampFreq/N;

% 模拟信号
x01=chirp(t,50,1,500);
x02=chirp(t,10,1,200);
x0=x01+x02;

x0_length=length(x0); % 原始信号长度

% 信号扩边
K=x0_length+1024; % 扩边后的信号长度
[x,k1,k2]=Signal_Extension_1D(x0,x0_length,K);

[imf,ort,nbits]=emd(x);

% 信号缩边
imfsize=size(imf);
imf_length=imfsize(1);
for i=1:imf_length
    [imf_cut(i,:)]=Signal_Cut_1D(imf(i,:),k1,k2);  
end
imf=imf_cut;

figure;
subplot(311)
plot(t,x01);title('y1');
subplot(312)
plot(t,x02);title('y2');
subplot(313)
plot(t,x0);title('y1+y2');

%{
%原始信号的绘制
figure;
plot(t,x);title('模拟信号');ylabel('幅值');xlabel('时间  t/s');

%IMF各分量
figure;
imfsize=size(imf);
i=imfsize(1);
for j=1:i
    [imf_cut(j,:)]=Signal_Cut_1D(imf(j,:),k1,k2);
    subplot(i/2+1,2,j);plot(t,imf(j,:));
    if(j<i) ylabel(sprintf('IMF%d', j));
    else ylabel('Res');
    end
end

imf=imf_cut
%分解后IMF剩余分量
figure;
newrdxn=x;
for k=1:i
    newrdxn=newrdxn-imf(k,:);
    subplot(i/2+1,2,k);plot(t,newrdxn);title(sprintf('IMF%d的剩余信号', k));
end 

%绘制原始信号的频谱
figure;
plot(T,abs(fft(x)));
title('频谱图');ylabel('幅值');xlabel('频率/HZ');
axis([0,sampFreq/2,0,1000]);%设置坐标轴的取值范围
%}
%希尔伯特谱  在图形的右上角添加标号如（a）
freq_resol=500;
[A,f,tt] = hhspectrum(imf(1:end,:));%修改imf(1:end-2,:)，可显示某一阶频率
[im,ttt,ff] = toimage(A,f,tt,length(tt),freq_resol);
[T,IM]=disp_hhs(im,[],N);%绘制真实频率，N为采样频率  ====IM为返回的参数====
colormap(flipud(gray));% 黑白显示,打印比较清晰
set(gca,'ytick',(0:100:N/2),'yticklabel',(0:100:N/2));%设置y坐标，前后设置要一致
set(gca,'FontName','Times New Roman');

%{
%边际谱
%%%边际谱表示信号中某一频率在各个时刻的幅值值之和，更能真实的反应频率在信号中是否存在
%将所有时刻某一瞬时频率的能量（幅值）加起来就是信号中该频率的总能量（总幅值），即边际谱线的高度
%%%从边际谱上可以看出信号频率分布的真实情况

fs=1;
aaa=zeros(length(size(im,1)));
for kkk=1:size(im,1)
    aaa(kkk)=sum(im(kkk,:))*1/fs;
end

figure;
plot(N*ff(1,:),aaa);%归一化频率乘于采样频率就是真实频率,N为采样频率
axis tight;
title('边际谱');%显示每个频率的总振幅
xlabel('频率');ylabel('幅值');
axis([0,sampFreq/4,0,1000]);%设置坐标轴的取值范围
%}
%{
%==================================================================================
% 保存成grapher、surfer用格式，用于绘图
Sig_length=length(x);
disp('正在保存本征模态函数数据……');

% 保存每一阶本征模态函数（IMF）成.dat格式=======
filename1='第一阶本征模态函数(IMF1).dat';  % 根据需要需改一下文件名
rdxn1=imf(1,:); % 根据需要需改一下文件名
fidD1=fopen(filename1,'wt');

for i=1:Sig_length
    fprintf(fidD1,'%.4f    %.6f',t(i),rdxn1(i)); 
    fprintf(fidD1,'\n');  % 换行
end
fclose(fidD1);
disp('本征模态函数数据保存完成！');

% 边界谱保存===================
disp('正保存边际谱数据……');

filename2='边界谱.dat';  
rdxn2=ff(1,:); 
fidD2=fopen(filename2,'wt');

bianji_length=length(ff(1,:));
for j=1:bianji_length
    fprintf(fidD2,'%.4f    %.6f',Sig_length*rdxn2(j),aaa(j)); % Sig_length为原始信号长度，也就是N
    fprintf(fidD2,'\n');  % 换行
end
fclose(fidD2);
disp('边际谱数据保存完成！');
%}
%{
% 保存hilber谱===============
% 方法一

disp('正保存希尔伯特谱数据……');

[yn,xn]=size(IM);

YN=yn/4;  % xiugai
XN=xn;
% 查找Z值极大值、极小值
M=zeros(XN*YN,1);
k=1;
for i0=1:YN
    for j0=1:XN
        if IM(i0,j0)==-Inf
            IM(i0,j0)=0;
        end
        M(k)=IM(i0,j0);
        k=k+1;
    end
end
Z_MAX=max(M);
Z_MIN=min(M);

%XX=[0,(1/sampFreq)*N];
X_MAX=max(T);
X_MIN=min(T);

YY=[0,sampFreq]; % sampFreq 采样频率  xiugai
Y_MAX=max(YY);
Y_MIN=min(YY);

filename3='希尔伯特谱.grd';
fidD3=fopen(filename3,'wt');
fprintf(fidD3,'%s','DSAA');
fprintf(fidD3,'\n');
fprintf(fidD3,'%d  %d',XN,YN);
fprintf(fidD3,'\n');
fprintf(fidD3,'%f  %f',X_MIN,X_MAX);
fprintf(fidD3,'\n');
fprintf(fidD3,'%f  %f',Y_MIN,Y_MAX);
fprintf(fidD3,'\n');
fprintf(fidD3,'%f  ',Z_MIN,Z_MAX);
fprintf(fidD3,'\n');

for i1=1:YN   %
    for j1=1:XN
        fprintf(fidD3,'%f  ',IM(i1,j1));
    end
    fprintf(fidD3,'\n');
end
fclose(fidD3);
disp('希尔伯特谱数据保存完成！');
%}

%{
% 方法二  保存成dat格式，再网格化
[YN,XN]=size(IM);

% 查找Z值极大值、极小值
M=zeros(XN*YN,1);
for i0=1:YN
    for j0=1:XN
        if IM(i0,j0)==-Inf
            IM(i0,j0)=0;
        end
    end
end
% x,y坐标矢量化
XX=linspace(0,(1/sampFreq)*N,XN);
YY=linspace(0,0.5*sampFreq,YN);
xx=XX;
yy=YY;

filename3='希尔伯特谱.dat';  
fidD3=fopen(filename3,'wt');


for i1=1:XN
    for j1=1:YN
        fprintf(fidD3,'%.4f    %.6f    %.6f',xx(i1),yy(j1),IM(j1,i1)); % Sig_length为原始信号长度，也就是N
    end
    fprintf(fidD3,'\n');  % 换行
end
fclose(fidD3);
%}










