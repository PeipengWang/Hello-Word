clear
clc
tic%计算运行时间
L=1e-2;%系统长
R=5e-3;%系统半径
I = -1/3*1.0e-1;%入射的每个宏电子电流
R0=2.0e-3;%最内圈宏电子注入半径
Nz=201;%纵向网格数
Nr=101;%径向网格数
Ntha=36;%每一圈电子选择的有限大小粒子数
dz=5.0e-5;%纵向空间步长
dr=5.0e-5;%径向空间步长
dtha=pi/18;%角向空间步长
dt=1.0e-13;%时间步长
e0=8.85e-12;%自由空间介电常数
q=-1.6e-19;%电子电荷量
qBm=-1.76e11;%电子的荷质比
Vol0=100;%左端面电压
T=6001;%总时间步数

numPar=0;%粒子个数
N=5;%五圈电子
num=180;%每次发射宏电子总数
vz0(1:num,1)=1.0e7;%纵向初始速度
vtha0(1:num,1)=0;%角向初始速度
vr0(1:num,1)=0;%径向初始速度 
r0(1:num,1)=0;%径向初始位置
tha0(1:num,1)=0;%径向初始位置
z0(1:num,1)=0;%纵向初始位置
for i=1:N
    for j=1:Ntha
        r0(j+(i-1)*Ntha,1)=R0+R0/10*(i-1);%径向初始位置
        tha0(j+(i-1)*Ntha,1)=dtha*(j-1);%角向初始位置 
    end
end
z0=z0+vz0*dt/2;

%初始化电荷密度
NthaAdd2=Ntha+2;
Rou=zeros(Nr,NthaAdd2,Nz);

%初始化电势
Voltage=zeros(Nr,NthaAdd2,Nz);
for k=1:Nz
    for j=1:NthaAdd2
        dV=Vol0/(Nz-1);
        Voltage(:,:,k)=(k-1)*dV;
    end
end
VoltageNew=Voltage;

figNum=0;
r_zR=zeros(T,N);
r_zZ=zeros(T,N);
hold off
for t=1:T 
    t %显示当前计算时间步数
    %========发射新粒子========
    if mod(t,10)==1
        %初始化粒子速度及位置
        vz(numPar+1:numPar+num,1)=vz0;
        vtha(numPar+1:numPar+num,1)=vtha0;
        vr(numPar+1:numPar+num,1)=vr0;
        r(numPar+1:numPar+num,1)=r0;
        tha(numPar+1:numPar+num,1)=tha0;
        z(numPar+1:numPar+num,1)=z0;        
        numPar=numPar+num;
        %为电场分量提供空间
        Ez=zeros(numPar,1);
        Er=zeros(numPar,1);
        Etha=zeros(numPar,1); 
    end
    %========发射新粒子—END========
    

    for numRZ=1:5
        r_zR(t,numRZ)=r(1+(numRZ-1)*36);
        r_zZ(t,numRZ)=z(1+(numRZ-1)*36);
    end

    %========更新电荷密度========

    Rou(:,:,:)=0;
    for n=1:numPar
        i=floor(r(n)/dr)+1;
        j=floor(tha(n)/dtha)+2;
        k=floor(z(n)/dz)+1;
        if k<0||j<0||i<0
            break;
        end  
        if i==1
            drou=I*dt/dtha/dz^3/(2*i-1)/3;
            Rou(i,j,k)=Rou(1,2,k)+drou;
            Rou(i+1,j,k)=Rou(i+1,j,k)+drou;
            Rou(i+1,j+1,k)=Rou(i+1,j+1,k)+drou;
            Rou(i,j,k+1)=Rou(1,2,k+1)+drou;
            Rou(i+1,j,k+1)=Rou(i+1,j,k+1)+drou;
            Rou(i+1,j+1,k+1)=Rou(i+1,j+1,k+1)+drou;              
        else
            drou=I*dt/dtha/dz^3/(2*i-1)/4;
            Rou(i,j,k)=Rou(i,j,k)+drou;
            Rou(i+1,j,k)=Rou(i+1,j,k)+drou;
            Rou(i,j+1,k)=Rou(i,j+1,k)+drou;
            Rou(i+1,j+1,k)=Rou(i+1,j+1,k)+drou;
            Rou(i,j,k+1)=Rou(i,j,k+1)+drou;
            Rou(i+1,j,k+1)=Rou(i+1,j,k+1)+drou;
            Rou(i,j+1,k+1)=Rou(i,j+1,k+1)+drou;
            Rou(i+1,j+1,k+1)=Rou(i+1,j+1,k+1)+drou;              
        end
   
    end
    Rou(:,2,:)=Rou(:,2,:)+Rou(:,NthaAdd2,:);
    Rou(:,NthaAdd2,:)=Rou(:,2,:);    

    %========更新电荷密度—END========
        
    %========更新电势========
    for k=2:Nz-1
        i=1;
        for j=1:NthaAdd2
            VoltageNew(i,j,k)=2/3*Voltage(2,2,k)+1/6*Voltage(1,2,k+1)+1/6*...
            Voltage(1,2,k-1)+dz^2/e0/6*Rou(1,2,k);
        end
        for i=2:Nr-1
            a0=1/(4+2/(i-1)^2/dtha^2);
            a1=a0*(1-1/(i-1)/2);
            a2=a0*(1+1/(i-1)/2);        
            a3=a0/(i-1)^2/dtha^2;       
            for j=2:Ntha+1
                VoltageNew(i,j,k)=a3*(Voltage(i,j-1,k)+Voltage(i,j+1,k))...
                +a1*Voltage(i-1,j,k)+a2*Voltage(i+1,j,k)+a0*...
                (Voltage(i,j,k-1)+Voltage(i,j,k+1))+dz^2*a0/e0*Rou(i,j,k);
            end
        end

        i=Nr;
        a0=1/(4+2/(i-1)^2/dtha^2);
        a1=a0*(1-1/(i-1)/2);
        a2=a0*(1+1/(i-1)/2);        
        a3=a0/(i-1)^2/dtha^2;       
        for j=2:Ntha+1
            VoltageNew(i,j,k)=a3*(Voltage(i,j-1,k)+Voltage(i,j+1,k))+...
            2*a0*Voltage(i-1,j,k)+a0*(Voltage(i,j,k-1)+Voltage(i,j,k+1))+...
            dz^2*a0/e0*Rou(i,j,k);
        end   
    end
    for k=1:Nz
        for i=1:Nr
            VoltageNew(i,1,k)=VoltageNew(i,Ntha+1,k);
            VoltageNew(i,NthaAdd2,k)=VoltageNew(i,2,k);
        end
    end
    Voltage=VoltageNew;
    Voltage(40,2,:);
    %========更新电势-END========
           
    %========更新电场========
    for n=1:numPar
        i=floor(r(n)/dr)+1;
        j=floor(tha(n)/dtha)+2;
        k=floor(z(n)/dz)+1;
        if i==1
            Ez(n)=(Voltage(i,j,k)+Voltage(i+1,j,k)+Voltage(i+1,j+1,k)-...
            Voltage(i,j,k+1)-Voltage(i+1,j,k+1)-Voltage(i+1,j+1,k+1))/dz/3;
            Er(n)=(Voltage(i,j,k)+Voltage(i,j,k+1)+Voltage(i,j+1,k)+...
            Voltage(i,j+1,k+1)-Voltage(i+1,j,k)-Voltage(i+1,j,k+1)-...
            Voltage(i+1,j+1,k)-Voltage(i+1,j+1,k+1))/dz/4;
            Etha(n)=(Voltage(i+1,j,k)-Voltage(i+1,j+1,k)+...
            Voltage(i+1,j,k+1)-Voltage(i+1,j+1,k+1))/dtha/i/dr/2;        
        else
            Ez(n)=(Voltage(i,j,k)+Voltage(i,j+1,k)+Voltage(i+1,j,k)+...
            Voltage(i+1,j+1,k)-Voltage(i,j,k+1)-Voltage(i,j+1,k+1)-...
            Voltage(i+1,j,k+1)-Voltage(i+1,j+1,k+1))/dz/4;
            Er(n)=(Voltage(i,j,k)+Voltage(i,j,k+1)+Voltage(i,j+1,k)+...
            Voltage(i,j+1,k+1)-Voltage(i+1,j,k)-Voltage(i+1,j,k+1)-...
            Voltage(i+1,j+1,k)-Voltage(i+1,j+1,k+1))/dz/4;
            Etha(n)=(Voltage(i,j,k)-Voltage(i,j+1,k)+Voltage(i,j,k+1)-...
            Voltage(i,j+1,k+1))/dtha/(i-1)/dr/4+...
            (Voltage(i+1,j,k)-Voltage(i+1,j+1,k)+Voltage(i+1,j,k+1)-...
            Voltage(i+1,j+1,k+1))/dtha/i/dr/4;
        end
    end
    
    %========更新粒子速度========
    vz=vz+qBm*dt*Ez;
    vrOld=vr;
    vr=vr+qBm*dt*Er+dt*vtha.^2./r;
    vtha=vtha+qBm*dt*Etha-dt*vrOld.*vtha./r;
    %=========更新粒子速度—END=========
    
    %========更新粒子位置========
    z=z+vz*dt;
    rOld=r;
    r=r+vr*dt;
    rOld=(rOld+r)/2;
    tha=tha+vtha*dt./rOld;
    tha=mod(tha,2*pi);%使角向位置位于0-2pi内
    %========更新粒子位置—END========
    
    %========更新图片显示========
    if rem(t,10)==0
    xx=rOld.*sin(tha);
    yy=rOld.*cos(tha);
    zz=z;
    pause(0.000001)
    subplot(2,2,1); plot3(xx,zz,yy,'.','MarkerSize',1), 
                    xlabel('x');ylabel('y'), zlabel('z'),  title('三维图') ;
    subplot(2,2,2);plot3(xx,yy,zz,'.','MarkerSize',1);
                    xlabel('x');ylabel('y'); zlabel('z'); title('三维图');
    subplot(2,2,3);plot(z,rOld,'.','MarkerSize',1);
                   xlabel('z');ylabel('r'); title('R-Z局部轨迹图');
    subplot(2,2,4);plot(z,rOld,'.','MarkerSize',1);
                   xlabel('z');ylabel('r');axis([0 0.01 0 0.005]); 
                   title('R-Z整体轨迹图');
    end 
    %========更新图片显示—END========
end

