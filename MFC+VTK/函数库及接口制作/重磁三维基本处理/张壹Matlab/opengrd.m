function [h,l,xmin,xmax,ymin,ymax,zmin,zmax,z,dx,dy]=opengrd(filename)
%函数功能：打开grd文件，返回相关信息
%函数返回参数为：行数，列数，x坐标最大值，最小值，y坐标最大值，最小值，z最小值，最大值，z，dx,dy
%[h,l,xmin,xmax,ymin,ymax,zmin,zmax,z,dx,dy]=opengrd(filename)
str=strcat(filename,'.grd');
fidx=fopen(str,'r'); %读入二维平面grd数据
cdumx=fread(fidx,4,'uint8=>char')'; %读取文件头DSAA
fscanf(fidx,'\n');
l=fscanf(fidx,'%d',1);%读取格网数据的x坐标个数（列数）,y坐标个数（行数）
fp=fseek(fidx,1,0);
h=fscanf(fidx,'%d',1);
fscanf(fidx,'\n');
xmin=fscanf(fidx,'%f',1);%读取x坐标的最小值和最大值
fseek(fidx,4,1);
xmax=fscanf(fidx,'%f',1);
fscanf(fidx,'\n');
ymin=fscanf(fidx,'%f',1);%读取y坐标的最小值和最大值
fseek(fidx,4,1);
ymax=fscanf(fidx,'%f',1);
fscanf(fidx,'\n');
zmin=fscanf(fidx,'%f',1);%读取z的最小值和最大值
fseek(fidx,4,1);
zmax=fscanf(fidx,'%f',1);
fscanf(fidx,'\n');
z=zeros(h,l);

fin=['读入',filename,'.grd'];
w=waitbar(0,'1','name',fin);
for i=1:h
    str=['已完成：',num2str(i*100/h),'%'];
    waitbar(i/h,w,str);
   for m=1:l
     z(i,m)=fscanf(fidx,'%f',1);
     fseek(fidx,4,1);
   end
end
fclose(fidx);
dx=((xmax-xmin)/(l-1));
dy=((ymax-ymin)/(h-1));
close(w);