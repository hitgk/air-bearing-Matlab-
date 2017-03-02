R1=25*10^(-3);%轴承的外半径
R2=50*10^(-3);%轴承的内半径
h1=60*10^(-6);
h2=10*10^(-6);
h11=h1/h2;
belta=90;
b=0.3;
bbelta=b*belta;
m=21;
n=41;
mm=7;
deltar=(1-R1/R2)/(n-1);
deltatheta=belta*pi/90*(m-1);
H=ones(m,n);
P=zeros(m,n);
P(1,:)=1;
p(m,:)=1;
P(:,1)=1;
P(:,n)=1;
%确定每个节点处的半径值r(i,j)
r1=R1/R2;
r=zeros(m,n);
for i=1:m
    for j=1:n
        r2=r1+(j-1)*deltar;
        r(i,j)=r2;
    end 
end 
r
%间隙分布初值
alpha=0.0429;
chi=13.44;
for i=1:1:m
    theta=(i-1)*belta/(m-1);
    for j=1:1:n
        if theta<bbelta
            g=(h11-1)*(1-theta/bbelta);
            H(i,j)=1+g;
        else
            H(i,j)=1;
         if P(i,j)<1
             H(i,j)=H(i,j);
         else
             H(i,j)=H(i,j)+alpha*(P(i,j)-1);
         end
        end
    end
end
H
figure(1);%初始间隙分布
[x,y]=meshgrid((r1:deltar:1),(0:belta/(m-1):belta));
mesh(x,y,H)
axis([0.5 1 0 90 0 6])
%计算压力分布
%迭代次数
S=0;T=0;
wucha=1;
error=10^(-6);
count=0;
w=1.8;
A1=ones(m,n);
A2=ones(m,n);
A3=ones(m,n);
A4=ones(m,n);
A5=ones(m,n);
A6=ones(m,n);
PP=ones(m,n);
P2=ones(m,n);
P3=ones(m,n);
H2=H;
while(wucha>=error)
    PP=P;
    for i=2:m-1
        for j=2:n-1
           P2(i,j)=P(i,j);
          A1(i,j)=(r(i,j)+deltar/2)*(H(i,j)^3+H(i,j+1)^3)/2*r(i,j)*deltar^2;
          A2(i,j)=(r(i,j)-deltar/2)*(H(i,j)^3+H(i,j-1)^3)/2*r(i,j)*deltar^2;
          A3(i,j)=(H(i+1,j)^3+H(i,j)^3)/2*r(i,j)^2*deltatheta^2;
          A4(i,j)=(H(i-1,j)^3+H(i,j)^3)/2*r(i,j)^2*deltatheta^2;
          p1=(P(i+1,j)+P(i,j))/2;
         p2=(P(i-1,j)+P(i,j))/2;
          p3=(P(i,j+1)+P(i,j))/2;
          p4=(P(i,j-1)+P(i,j))/2;
         A5(i,j)=chi*p1*(H(i+1,j)+H(i,j))/2*deltatheta-chi*p2*(H(i-1,j)+H(i,j))/2*deltatheta;
         P3(i,j)=A1(i,j)*p3*P(i,j+1)+A2(i,j)*p4*P(i,j-1)+A3(i,j)*p1*P(i+1,j)+A4(i,j)*p2*P(i-1,j)-A5(i,j);
         A6(i,j)=A1(i,j)*p3+A2(i,j)*p4+A3(i,j)*p1+A4(i,j)*p2;
         P3(i,j)=P3(i,j)/A6(i,j);
          P(i,j)=(1-w)*P2(i,j)+w*P3(i,j);
            if P(i,j)<0
                P(i,j)=0;
            end
            if P(i,j)>1
             H(i,j)=H(i,j)+alpha*(P(i,j)-1);
            else
                H(i,j)=H2(i,j);
            end
        end
    end
    for i=2:m-1
        for j=2:n-1
            S=S+P(i,j)-PP(i,j);
            T=T+P(i,j);
        end
    end
    wucha=S/T;
    count=count+1
end
 P
 figure(2);
[x,y]=meshgrid((r1:deltar:1),(0:belta/(m-1):belta));
mesh(x,y,P)
axis([0.5 1 0 90 0 100])
H
figure(3);
[x,y]=meshgrid((r1:deltar:1),(0:belta/(m-1):belta));
mesh(x,y,H)
axis([0.5 1 0 90 0 6000000])
         
           
          
           
          
