

clear;
clc;

filename='立方体磁异常.grd';
% [h,l,ymin,ymax,xmin,xmax,zmin,zmax,z,dx,dy]=opengrd(filename);%打开DeltaT文件

I0=45;%给出地磁场参数
A0=0;
%调用函数
savepath='E:\Study\重磁数据处理软件系统开发\程序模块及系统\函数库及接口制作\重磁三维基本处理\石达三分量转换及求梯度';
cal_mag_trans_new_ext(filename,I0,A0,2,savepath,'grd',1);