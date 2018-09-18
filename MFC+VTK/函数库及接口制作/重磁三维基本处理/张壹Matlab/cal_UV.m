function [U,V]=cal_UV(Yn,Xn,dx,dy)
Xn2=Xn/2+1;
du=1/((Xn2-1)*dx);
U=zeros(Yn,Xn);
for i=1:1:Yn
  for row=1:1:Xn2
    u=du*(row-1);
    U(i,row)=u;
    if row~=1&&row~=Xn2
        k=Xn+2-row;
        U(i,k)=-u;
    end
  end
end

Yn2=Yn/2+1;
dv=1/((Yn2-1)*dy);
V=zeros(Yn,Xn);
for i=1:1:Xn
  for row=1:1:Yn2
    v=dv*(row-1);
    V(row,i)=v;
    if row~=1&&row~=Yn2
        k=Yn+2-row;
        V(k,i)=-v;
    end
  end
end
end
