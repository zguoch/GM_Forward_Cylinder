function [T,k1,k2,k3,k4]=kuobian_duoxiangshi(t,H,L,op2,op,n)
%函数功能：局部多项式扩边，一维或二维
%函数返回值：扩边后数据（矩阵形式）
%函数参数：t待扩边数据；H、L扩边后行列数，一维是只需要H；op2=1计算一维
%op2=2计算二维；n多项式阶数
%op=1端点值为0，op=2端点值为原数据端点值平均值
%[T,k1,k2,k3,k4]=kuobian_duoxiangshi(t,H,L,op2,op,n)
if op2==1
    m=numel(t);
    M=H;
    T=zeros(M,1);

    k1=floor((M-m)/2);
    k2=ceil((M-m)/2);

    K1=ones(n+1,n+1);
    for i=2:n+1
        for j=2:n+1
            K1(i,j)=(k1+i-1)^(j-1);
        end
    end

    G1=zeros(n+1,1);
    if op==1
        G1(1)=0;
        else if op==2
             G1(1)=(t(1)+t(m))/2;
            end
    end

    for i=2:n+1
        G1(i)=t(i-1);
    end
    A1=K1\G1;

    T1=ones(k1,n+1);
    for i=2:k1
        for j=2:n+1
            T1(i,j)=i^(j-1);
        end
    end
    T2=T1*A1;
    for i=1:k1
        T(i)=T2(i);
    end

    K2=ones(n+1,n+1);
    for i=2:n+1
        for j=2:n+1
            K2(i,j)=(k2+i-1)^(j-1);
        end
    end
    
    G2=zeros(n+1,1);
    if op==1
        G2(1)=0;
    else if op==2
            G2(1)=(t(1)+t(m))/2;
        end
    end

    for i=2:n+1
        G2(i)=t(m-i+2);
    end
    A2=K2\G2;

    T3=ones(k2,n+1);
    for i=2:k2
        for j=2:n+1
            T3(i,j)=i^(j-1);
        end
    end
    T4=T3*A2;
    for i=1:k2
        T(i+k1+m)=T4(k2-i+1);
    end

    for i=k1+1:k1+m
        T(i,1)=t(i-k1);
    end
else if op2==2
        sizet=size(t);
ht=sizet(1,1);
lt=sizet(1,2);
Xn=L;%扩边到2^n个数据
Yn=H;

T=zeros(Yn,Xn);

k1=floor((Xn-lt)/2);
k2=ceil((Xn-lt)/2);
k3=floor((Yn-ht)/2);
k4=ceil((Yn-ht)/2);

for h=1:ht %向左扩边
    K1=ones(n+1,n+1);
    for i=2:n+1
        for j=2:n+1
            K1(i,j)=(k1+i-1)^(j-1);
        end
    end
    
    G1=zeros(n+1,1);
    if op==1
        G1(1)=0;
    else if op==2
            G1(1)=(t(h,1)+t(h,lt))/2;
        end
    end
    
    for i=2:n+1
        G1(i)=t(h,i-1);
    end
    A1=K1\G1;
    T1=ones(k1,n+1);
    for i=2:k1
        for j=2:n+1
            T1(i,j)=i^(j-1);
        end
    end
    T2=T1*A1;
    for i=1:k1
        T(h+k4,i)=T2(i);
    end
end

for h=1:ht %向右扩边
    K2=ones(n+1,n+1);
    for i=2:n+1
        for j=2:n+1
            K2(i,j)=(k2+i-1)^(j-1);
        end
    end
    
    G2=zeros(n+1,1);
    if op==1
        G2(1)=0;
    else if op==2
            G2(1)=(t(h,1)+t(h,lt))/2;
        end
    end
    
    for i=2:n+1
        G2(i)=t(h,lt-i+2);
    end
    A2=K2\G2;
    T3=ones(k2,n+1);
    for i=2:k2
        for j=2:n+1
            T3(i,j)=i^(j-1);
        end
    end
    T4=T3*A2;
    for i=1:k2
        T(h+k4,i+k1+lt)=T4(k2-i+1);
    end
end

for h=1:ht
    for l=1:lt
        T(h+k4,l+k1)=t(h,l);
    end
end

for l=1:Xn %向上扩边
    K3=ones(n+1,n+1);
    for i=2:n+1
        for j=2:n+1
            K3(i,j)=(k4+i-1)^(j-1);
        end
    end
    
    G3=zeros(n+1,1);
    if op==0
        G3(1)=0;
    else if op==2
            G3(1)=(T(1+k4,l)+T(ht+k4,l))/2;
        end
    end
    
    for i=2:n+1
        G3(i)=T(k4+i-1,l);
    end
    A3=K3\G3;
    T5=ones(k4,n+1);
    for i=2:k4
        for j=2:n+1
            T5(i,j)=i^(j-1);
        end
    end
    T6=T5*A3;
    for i=1:k4
        T(i,l)=T6(i);
    end
end

for l=1:Xn %向下扩边
    K4=ones(n+1,n+1);
    for i=2:n+1
        for j=2:n+1
            K4(i,j)=(k3+i-1)^(j-1);
        end
    end
    
    G4=zeros(n+1,1);
    if op==0
        G4(1)=0;
    else if op==2
            G4(1)=(T(1+k4,l)+T(ht+k4,l))/2;
        end
    end
    
    for i=2:n+1
        G4(i)=T(k4+ht-i+2,l);
    end
    A4=K4\G4;
    T7=ones(k3,n+1);
    for i=2:k3
        for j=2:n+1
            T7(i,j)=i^(j-1);
        end
    end
    T8=T7*A4;
    for i=1:k3
        T(i+k4+ht,l)=T8(k3-i+1);
    end
end
    end
end