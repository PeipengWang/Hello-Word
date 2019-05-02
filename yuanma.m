clear
clc
tic%��������ʱ��
L=1e-2;%ϵͳ��
R=5e-3;%ϵͳ�뾶
I = -1/3*1.0e-1;%�����ÿ������ӵ���
R0=2.0e-3;%����Ȧ�����ע��뾶
Nz=201;%����������
Nr=101;%����������
Ntha=36;%ÿһȦ����ѡ������޴�С������
dz=5.0e-5;%����ռ䲽��
dr=5.0e-5;%����ռ䲽��
dtha=pi/18;%����ռ䲽��
dt=1.0e-13;%ʱ�䲽��
e0=8.85e-12;%���ɿռ��糣��
q=-1.6e-19;%���ӵ����
qBm=-1.76e11;%���ӵĺ��ʱ�
Vol0=100;%������ѹ
T=6001;%��ʱ�䲽��

numPar=0;%���Ӹ���
N=5;%��Ȧ����
num=180;%ÿ�η�����������
vz0(1:num,1)=1.0e7;%�����ʼ�ٶ�
vtha0(1:num,1)=0;%�����ʼ�ٶ�
vr0(1:num,1)=0;%�����ʼ�ٶ� 
r0(1:num,1)=0;%�����ʼλ��
tha0(1:num,1)=0;%�����ʼλ��
z0(1:num,1)=0;%�����ʼλ��
for i=1:N
    for j=1:Ntha
        r0(j+(i-1)*Ntha,1)=R0+R0/10*(i-1);%�����ʼλ��
        tha0(j+(i-1)*Ntha,1)=dtha*(j-1);%�����ʼλ�� 
    end
end
z0=z0+vz0*dt/2;

%��ʼ������ܶ�
NthaAdd2=Ntha+2;
Rou=zeros(Nr,NthaAdd2,Nz);

%��ʼ������
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
    t %��ʾ��ǰ����ʱ�䲽��
    %========����������========
    if mod(t,10)==1
        %��ʼ�������ٶȼ�λ��
        vz(numPar+1:numPar+num,1)=vz0;
        vtha(numPar+1:numPar+num,1)=vtha0;
        vr(numPar+1:numPar+num,1)=vr0;
        r(numPar+1:numPar+num,1)=r0;
        tha(numPar+1:numPar+num,1)=tha0;
        z(numPar+1:numPar+num,1)=z0;        
        numPar=numPar+num;
        %Ϊ�糡�����ṩ�ռ�
        Ez=zeros(numPar,1);
        Er=zeros(numPar,1);
        Etha=zeros(numPar,1); 
    end
    %========���������ӡ�END========
    

    for numRZ=1:5
        r_zR(t,numRZ)=r(1+(numRZ-1)*36);
        r_zZ(t,numRZ)=z(1+(numRZ-1)*36);
    end

    %========���µ���ܶ�========

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

    %========���µ���ܶȡ�END========
        
    %========���µ���========
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
    %========���µ���-END========
           
    %========���µ糡========
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
    
    %========���������ٶ�========
    vz=vz+qBm*dt*Ez;
    vrOld=vr;
    vr=vr+qBm*dt*Er+dt*vtha.^2./r;
    vtha=vtha+qBm*dt*Etha-dt*vrOld.*vtha./r;
    %=========���������ٶȡ�END=========
    
    %========��������λ��========
    z=z+vz*dt;
    rOld=r;
    r=r+vr*dt;
    rOld=(rOld+r)/2;
    tha=tha+vtha*dt./rOld;
    tha=mod(tha,2*pi);%ʹ����λ��λ��0-2pi��
    %========��������λ�á�END========
    
    %========����ͼƬ��ʾ========
    if rem(t,10)==0
    xx=rOld.*sin(tha);
    yy=rOld.*cos(tha);
    zz=z;
    pause(0.000001)
    subplot(2,2,1); plot3(xx,zz,yy,'.','MarkerSize',1), 
                    xlabel('x');ylabel('y'), zlabel('z'),  title('��άͼ') ;
    subplot(2,2,2);plot3(xx,yy,zz,'.','MarkerSize',1);
                    xlabel('x');ylabel('y'); zlabel('z'); title('��άͼ');
    subplot(2,2,3);plot(z,rOld,'.','MarkerSize',1);
                   xlabel('z');ylabel('r'); title('R-Z�ֲ��켣ͼ');
    subplot(2,2,4);plot(z,rOld,'.','MarkerSize',1);
                   xlabel('z');ylabel('r');axis([0 0.01 0 0.005]); 
                   title('R-Z����켣ͼ');
    end 
    %========����ͼƬ��ʾ��END========
end

