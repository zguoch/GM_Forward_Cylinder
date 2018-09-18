function [XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,Z_MIN,Z_MAX,data]=open_grd(filename)
%打开并读取GRD文件
fid = fopen(filename,'r');
ID=fgets(fid); 
XN=fscanf(fid,'%d',1);
YN=fscanf(fid,'%d',1);
X_MIN=fscanf(fid,'%f',1);
X_MAX=fscanf(fid,'%f',1);
Y_MIN=fscanf(fid,'%f',1);
Y_MAX=fscanf(fid,'%f',1);
Z_MIN=fscanf(fid,'%f',1);
Z_MAX=fscanf(fid,'%f',1);
data=fscanf(fid,'%f',[XN,YN]);
data=data';
fclose(fid);