epison=0.4;
nx=180;
ny=100;
deltax=2*180/nx;
deltay=2/ny;
d=75;
L=50;
lamta=d/L;
lam=(lamta*deltax*pi/180*deltay)^2;
P=1.0*ones(nx+1,ny+1);
P(1,:)=1;
P(nx+1,:)=1;
P(:,1)=1;
P(:,ny+1)=1;
H=ones(nx+1,ny+1);
A=ones(nx+1,ny+1);
B=ones(nx+1,ny+1);
C=ones(nx+1,ny+1);
D=ones(nx+1,ny+1);
E=ones(nx+1,ny+1);
F=ones(nx+1,ny+1);
P2=ones(nx+1,ny+1);
P3=ones(nx+1,ny+1);
%计算间隙初始值
for i=1:1:nx+1
     theta=(i-1)*deltax*pi/180;
        for j=1:1:ny+1
            H(i,j)=1+epison*cos(theta);
        end 
end
H
figure(1);%初始间隙分布
[x,y]=meshgrid((0:deltay:2),(0:deltax:2*180));
mesh(x,y,H)
axis([0 2 0 2*180 0 1.5])
S=0;T=0;
wucha=1;
error=10^(-3);
count=0;
w=1.75;
while(wucha>=error)
    PP=P;
    for i=2:nx
        for j=2:ny
            P2(i,j)=P(i,j);
            A(i,j)=H(i+1,j)^3;
            B(i,j)=H(i-1,j)^3;
            C(i,j)=lam*H(i,j+1)^3;
            D(i,j)=lam*H(i,j-1)^3;
            E(i,j)=A(i,j)+B(i,j)+C(i,j)+D(i,j);
            F(i,j)=3*(deltax*pi/180)*(H(i+1,j)-H(i-1,j));
            P3(i,j)=(A(i,j)*P(i+1,j)+B(i,j)*P(i-1,j)+C(i,j)*P(i,j+1)+D(i,j)*P(i,j-1)-F(i,j))/E(i,j);
            P(i,j)=(1-w)*P2(i,j)+w*P3(i,j);
            if P(i,j)<1
                P(i,j)=1;
            end
        end
    end
    for i=2:nx
        for j=2:ny
            S=S+P(i,j)-PP(i,j);
            T=T+P(i,j);
        end
    end
    wucha=S/T;
    count=count+1
end 
P
figure(2);%压力分布图
[x,y]=meshgrid((0:deltay:2),(0:deltax:2*180));
mesh(x,y,P)
axis([0 2 0 2*180 0 3])


    
            