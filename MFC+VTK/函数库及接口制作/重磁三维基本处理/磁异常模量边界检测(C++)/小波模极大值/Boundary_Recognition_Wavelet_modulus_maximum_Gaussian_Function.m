
clear; 
clc;

%利用小波模极大值进行位场边界识别，提取线性结构
% 1 什么是模极大值？一般信号的主要信息，由拐点（二阶导数为零的点）确定，而由于噪声的影响，直接求拐点显然困难。于是，我们求一阶导数的模的极大值。
% 
% 2 什么是小波模极大值？就是先将小波函数和原信号卷积（连续小波变换），然后对结果取模，最后找到极大值。上述步骤，也就等价于：先把某一光滑函数求导（求导后满足积分为零的条件成为小波函数），然后卷积源信号，接着取模，最后发现极大值。
% 
% 3 图像处理的操作。
% 
%     a、给定某一尺度，求出二维高斯函Phi_y。数沿x和沿y方向的导数Phi_x,这两个函数就等价于小波函数。
% 
%     b、用Phi_x,Phi_y分别与图像卷积得到Gx，Gy。
% 
%     c、求出每一个像素点的梯度大小G=(Gx*Gx+Gy*Gy).^(1/2)，用反正切求梯度方向或者称幅角atan(Gy/Gx)。这里，注意的是反正切只能求出一、四象限的角度，其它象限要分别处理。且Gx为一个很小的数值时，也要处理。
% 
%     d、把求得幅角，分成四种方向。第一种0或180方向（水平），第二种90或270方向（垂直），第三种45或225方向（正对角线），第四种135或315方向（负对角线）。也就是说，看看你求出幅角的大小与上面的哪个方向最接近。
% 
%     e、依次检测每一个像素点，看看在它对应“幅角最接近的方向上”是否是极大值。如果是，纪录该梯度值。若不是，把梯度值置零。
% 
%     f、找到记录梯度值中的最大值，然后以该值做归一化。比较每一个像素归一化的梯度值，当该梯度值大于某个阈值的时候，就是真正边缘，否则认为是伪边缘。
% 
% 4 实际上这个算法和canny算子本质上等价的。让我们再来回顾canny本人经典的原话，来体会边缘提取的目标到底是什么。
% 
%     a、好的检测性能。不漏检真实边缘，也不把非边缘点作为边缘点检出，使输出的信噪比最大。
% 
%     b、好的定位性能。检测到的边缘点与实际边缘点位置最近。
% 
%     c、唯一性。对于单个边缘点仅有一个响应。

tic;
filename='局部异常.grd'; %读取文件
fname=filename(1:size(filename,2)-4);

% 打开并读取GRD文件
[XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,Z_MIN,Z_MAX,data]=open_grd(filename);

[SIZE_Y,SIZE_X]=size(data);% 获取信号大小
X_grd= (X_MAX-X_MIN)/(XN-1); % x方向点距
Y_grd= (Y_MAX-Y_MIN)/(YN-1); % y方向点距

% 待修改参数：
% m ; 滤波器长度N; 阈值：threshold

% 在位场边界识别中，m,N可设定某一固定的数值，只需改变阈值threshold的数值

% 多尺度
m=0;               % 阶次,可经过调试，此数值可保持不变
delta=2^m;

% 构造高斯函数x,y方向偏导
N=20;  % 滤波器长度（需要调整，必须是偶数）不宜太长
A=-1/sqrt(2*pi);  % 幅度
phi_x=zeros(N,N);
phi_y=zeros(N,N);
for index_y=1:N;
    for index_x=1:N;
        x=index_x-(N+1)/2;
        y=index_y-(N+1)/2;
        phi_x(index_x,index_y)=A*(x/delta^2).*exp(-(x.*x+y.*y)/(2*delta^2));%对x偏导
        phi_y(index_x,index_y)=A*(y/delta^2).*exp(-(x.*x+y.*y)/(2*delta^2));%对y偏导
    end
end;

phi_x=phi_x/norm(phi_x);  %  能量归一化
phi_y=phi_y/norm(phi_y);  %  能量归一化

%  对二维信号做行列卷积
Gx=conv2(data,phi_x,'same');
Gy=conv2(data,phi_y,'same');% 参数"same"返回二维卷积结果中与data大小相同的部分
%"full"返回二维卷积的全部结果
%"valid"返回在卷积过程中，未使用边缘补 0 部分进行计算的卷积结果部分，
%当 size(A)>size(B) 时，size(C)=[Ma-Mb+1,Na-Nb+1].

% 求梯度
Grads=zeros(SIZE_Y,SIZE_X);
for i0=1:SIZE_Y
    for j0=1:SIZE_X
          Grads(i0,j0)=sqrt((Gx(i0,j0)*Gx(i0,j0))+(Gy(i0,j0)*Gy(i0,j0)));
    end
end
M=fix(Grads);% 用于获得直方图的参数，进而得到最大最小阈值

% 求幅角（梯度方向）
angle_array=zeros(SIZE_X,SIZE_Y);  %  角度

for i=1:SIZE_Y;
    for j=1:SIZE_X
        if (abs(Gx(i,j))>eps*100)  %  x的绝对值足够大
            p=atan(Gy(i,j)/Gx(i,j))*180/pi;  %  反正切求角度值(1,4象限)
            if (p<0)        %  负的幅角（4象限）
                p=p+360;
            end;
            if (Gx(i,j)<0 && p>180)     %  2象限的特殊处理
                p=p-180;
            elseif (Gx(i,j)<0 && p<180) %  3象限的特殊处理
                p=p+180;
            end
        else  %  90或270度
            p=90;
        end
        angle_array(i,j)=p;  %  幅角赋值
    end
end;

% 找边缘
edge_array=zeros(SIZE_Y,SIZE_X);

% 遍历,沿着梯度方向检测小波变换系数模的局部极大值点
for i=2:SIZE_Y-1
    for j=2:SIZE_X-1
        if ((angle_array(i,j)>=(-22.5) && angle_array(i,j)<=22.5) || ...
            (angle_array(i,j)>=(180-22.5) && angle_array(i,j)<=(180+22.5)))     %  0/180
            if (Grads(i,j)>Grads(i+1,j) && Grads(i,j)>Grads(i-1,j))
                edge_array(i,j)=Grads(i,j);
            end
        elseif ((angle_array(i,j)>=(90-22.5) && angle_array(i,j)<=(90+22.5)) || ...
                (angle_array(i,j)>=(270-22.5) && angle_array(i,j)<=(270+22.5))) %  90/270
            if (Grads(i,j)>Grads(i,j+1) && Grads(i,j)>Grads(i,j-1))
                edge_array(i,j)=Grads(i,j);
            end
        elseif ((angle_array(i,j)>=(45-22.5) && angle_array(i,j)<=(45+22.5)) || ...
                (angle_array(i,j)>=(225-22.5) && angle_array(i,j)<=(225+22.5))) %  45/225
            if (Grads(i,j)>Grads(i+1,j+1) && Grads(i,j)>Grads(i-1,j-1))
                edge_array(i,j)=Grads(i,j);
            end
        else  %  135/215
            if (Grads(i,j)>Grads(i+1,j-1) && Grads(i,j)>Grads(i-1,j+1))
                edge_array(i,j)=Grads(i,j);
            end
        end
    end
end

% 去除伪边缘
MAX_E=max(max(edge_array).');     % 最大幅度值
edge_array=edge_array/MAX_E;      % 归一化

% 根据直方图获取高低阈值 同canny算子相似
count=0;% 累计直方图统计可能是边缘点的个数
NUM=2048;
hist=zeros(1,NUM);
for k1=1:SIZE_Y
    for k2=1:SIZE_X
        if (edge_array(k1,k2)~=0)
            hist(1,M(k1,k2)+1)=hist(1,M(k1,k2)+1)+1;
            count=count+1;
        end
    end
end
for k3=1:NUM
    if hist(1,k3)~=0
        nmaxmag=k3;
    end
end
p=0.7;
dot=ceil(p*count);
k4=1;
nedgenb=hist(1,1);
while (k4<(nmaxmag-1)&&(nedgenb<dot))
    k4=k4+1;
    nedgenb=nedgenb+hist(1,k4);
end
k41=k4/50; % 要获得更多细节信息，可将分母数值增大
MAX_threshold=k41;                                       % 高阈值
MIN_threshold=0.4*MAX_threshold;                         % 低阈值,一般高为低的2.5倍

% 遍历边缘点判定
for i5=1:SIZE_Y
    for j5=1:SIZE_X
        if (edge_array(i5,j5)>=MAX_threshold)                 % 高于高阈值一定是边缘
            edge_array(i5,j5)=1;
        elseif (edge_array(i5,j5)<=MIN_threshold)              % 低于低阈值一定不是边缘
            edge_array(i5,j5)=0;
        else                                   % 高低阈值之间判断8邻域是否高于阈值否则不是边缘
            if((edge_array(i5-1,j5)>MAX_threshold)&&(edge_array(i5,j5-1)>MAX_threshold)&&(edge_array(i5+1,j5)>=MAX_threshold)&&(edge_array(i5,j5+1)>=MAX_threshold)...
                    &&(edge_array(i5-1,j5-1)>MAX_threshold)&&(edge_array(i5+1,j5-1)>MAX_threshold)&&(edge_array(i5-1,j5+1)>MAX_threshold)&&(edge_array(i5+1,j5+1)>MAX_threshold))
               edge_array(i5,j5)=1;
            else
               edge_array(i5,j5)=0;
            end
        end
    end
end

% 图像上、下翻转
edge_array_new=zeros(SIZE_Y,SIZE_X);
for i0=1:SIZE_Y
    for j0=1:SIZE_X
         edge_array_new(i0,j0)=edge_array(SIZE_Y-i0+1,j0);         
    end        
end

%图像显示
figure;
imshow(edge_array_new);
title(strcat('小波模极大值边界识别结果','（','阈值','threshold=',num2str(MIN_threshold),'-',num2str(MAX_threshold),'）'));

%{
% 保存GRD数据
Z_MIN_NEW=min(min(edge_array_new));
Z_MAX_NEW=max(max(edge_array_new));
%filename1=strcat(fname,'_Boundary_Recognition_Wavelet_modulus_maximum','_threshold=',num2str(MIN_threshold),'-',num2str(MAX_threshold),'.grd');
filename1=strcat(fname,'_threshold=',num2str(MIN_threshold),'-',num2str(MAX_threshold),'.grd');
save_grd(filename1,XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,Z_MIN_NEW,Z_MAX_NEW,edge_array_new);
%}
% 保存.dat数据（把数值为1的边界点全部保存）
filename2=strcat(fname,'_边界位置','.dat'); 
fidD=fopen(filename2,'wt'); 
for m=1:SIZE_Y
    for n=1:SIZE_X
        if edge_array(m,n)==1
            fprintf(fidD,'%.6f   %.6f   %.2f',(X_MIN+(n-1)*X_grd),(Y_MIN+(m-1)*Y_grd),edge_array(m,n));
            fprintf(fidD,'\n');
        end
    end
end      
fclose(fidD);
%}
toc;
