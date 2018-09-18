function save_grd(filename,XN,YN,X_MIN,X_MAX,Y_MIN,Y_MAX,Z_MIN,Z_MAX,data)
%打开并读取GRD文件
fid = fopen(filename,'wt');
fprintf(fid, '%s', 'DSAA'); 
fprintf(fid, '\n');
fprintf(fid,'%d ',XN);
fprintf(fid,'%d\n',YN);
fprintf(fid,'%f ',X_MIN);
fprintf(fid,'%f\n',X_MAX);
fprintf(fid,'%f ',Y_MIN);
fprintf(fid,'%f\n',Y_MAX);
fprintf(fid,'%.5f ',Z_MIN);
fprintf(fid,'%.5f\n',Z_MAX);
for i=1:YN;
    for j=1:XN;
        fprintf(fid,'%f ',data(i,j));
    end
    fprintf(fid,'\n');
end    
fclose(fid);