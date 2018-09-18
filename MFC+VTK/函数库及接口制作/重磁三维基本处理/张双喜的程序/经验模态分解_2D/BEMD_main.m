
clear ;
clc;

% 利用二维经验模态分解进行多尺度分解、位场分离
tic;      % tic toc可以用来计算运算时间
filename='长方体组合模型_with_noise_2.5%.grd'; % 读取文件
fname=filename(1:size(filename,2)-4);% size(filename,2)表示2维数据的长度，size(filename,2)-4就是出去.grd的文件名
[XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,Z_MIN,Z_MAX,data]=open_grd(filename);% 读取网格数据
[ imf_matrix residue YuZhi] = bemd(data);% bemd分解
M=size(imf_matrix);% 返回值为[57,100,5];行数YN为57，列数XN为100，分解阶次K为5
N=M(3);% 分解阶次

%  分解的各阶次细节（固有模态函数BIMF）
for i=1:(N-1)
filename=strcat(fname,'_ε=',num2str(YuZhi),'_BIMF_','D',int2str(i),'.grd'); % 'D'表示细节Datail
z=(imf_matrix(:,:,i));
zmax=max(z(:));
zmin=min(z(:));
save_grd(filename,XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,zmin,zmax ,imf_matrix(:,:,i));
end
%  分解的各阶次逼近
for j=2:N  % 分解一次后才有剩余信号，                    
filename=strcat(fname,'_ε=',num2str(YuZhi),'_BIMF_','A',int2str(j-1),'.grd'); %‘A’表示Approximate
z1=(residue(:,:,j));
zmax1=max(z1(:));
zmin1=min(z1(:));
save_grd(filename,XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,zmin1,zmax1 ,residue(:,:,j));
end  

%
% 对各模态函数叠加重构实现分场
% 将分解得到的前K阶异常进行叠加，叠加结果作为局部异常
Num1=3; 
filename1=strcat(fname,'_','ε=',num2str(YuZhi),'局部异常（','BIMF1','～','BIMF',int2str(Num1),'）','.grd');
Z1=zeros(YN,XN);% 创建YN行XN列的零向量

for i1=1:Num1
    Z1=Z1+(imf_matrix(:,:,i1));
end

%{
% 将剩余分量叠加重构，叠加结果作为区域异常
zmax=max(Z1(:));
zmin=min(Z1(:));
save_grd(filename1,XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,zmin,zmax,Z1);

Num2=Num1+1;
filename2=strcat(fname,'_','ε=',num2str(YuZhi),'区域异常（','BIMF',int2str(Num2),'～','BIMF',int2str(N),'）','.grd');
Z2=zeros(YN,XN);%创建YN行XN列的零向量
for j1=Num2:N
     Z2=Z2+(imf_matrix(:,:,j1));
end    
zmax=max(Z2(:));
zmin=min(Z2(:));
save_grd(filename2,XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,zmin,zmax,Z2);
%}
toc;