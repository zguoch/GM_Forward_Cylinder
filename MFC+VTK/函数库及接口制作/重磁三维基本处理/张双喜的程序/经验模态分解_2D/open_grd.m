function [XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,Z_MIN,Z_MAX,data]=open_grd(filename)
%打开并读取GRD文件
fid = fopen(filename,'r');
ID=fgets(fid); 
XN=fscanf(fid,'%d',1);%网格文件列数
YN=fscanf(fid,'%d',1);%网格文件行数
X_MIN=fscanf(fid,'%f',1);%x,y,z的变换范围
X_MAX=fscanf(fid,'%f',1);
Y_MIN=fscanf(fid,'%f',1);
Y_MAX=fscanf(fid,'%f',1);
Z_MIN=fscanf(fid,'%f',1);
Z_MAX=fscanf(fid,'%f',1);
data=fscanf(fid,'%f',[XN,YN]);%至少读取XN*YN个元素
data=data';
fclose(fid);