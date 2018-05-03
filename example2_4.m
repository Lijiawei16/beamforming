%**************example2.4***********************************************
close all
clear all
clc


N=10;
lambda=1;
noise = 1;
sig1 = 10^1.5;
theta1=-10;  
v1=[-sind(theta1),cosd(theta1)];
k1=2*pi/lambda*v1;

sig2 = 10^3;
theta2=10;
v2=[-sind(theta2),cosd(theta2)];
k2=2*pi/lambda*v2;

%*************************基阵位置******************************************
position=zeros(2,N);
position(2,:)=0;
for mm=1:N
    position(1,mm)=(mm-1)*lambda./2;   
    as1(mm)=exp(-i*(k1*position(:,mm)));   %信号1 15dB -10度基阵响应向量
    as2(mm)=exp(-i*(k2*position(:,mm)));   %信号2 30dB 10度基阵响应向量
end
as1=as1';%列向量
as2=as2';
as=as1+as2;
Rx0 =sig1*as1*as1'+sig2*as2*as2'+noise*eye(N);%+ sigma_i*a_i*a_i';  %接收信号协方差矩阵

w_beamforming=inv(Rx0)*as/(as'*inv(Rx0)*as);  %MVDR波束加权向量
w2_beamforming=as/N;  %常规波束加权向量

%*******************************MVDR方位谱**********************************
theta_0=-90:1:90;
P=zeros(1,length(theta_0));
a=zeros(N,1);
for mm=1:length(theta_0)
    vd=[-sind(theta_0(mm)),cosd(theta_0(mm))];
    kd=2*pi./lambda*vd;
    for nn=1:N
        a(nn)=exp(-i*(kd*position(:,nn)));
    end
%     Rx = sig1*as*as'+noise*eye(N);
   w_beamforming=inv(Rx0)*a/(a'*inv(Rx0)*a); 
  ww(mm,:)=inv(Rx0)*a/(a'*inv(Rx0)*a); 
   P(mm)=w_beamforming'*Rx0*w_beamforming;
end
figure(2)
plot(theta_0,10*log10(abs(P)),'LineWidth',1.5);
grid on;
% xlim([-90,90]);ylim([-15,40])
xlabel('方位/(^o)');ylabel('方位谱/dB')
title('方位谱')
%**************************************************************************
%*******************************常规方位谱**********************************
 w_beamforming= w2_beamforming;
theta_0=-90:1:90;
% figure
P=zeros(1,length(theta_0));
a=zeros(N,1);
for mm=1:length(theta_0)
    vd=[-sind(theta_0(mm)),cosd(theta_0(mm))];
    kd=2*pi./lambda*vd;
    for nn=1:N
        a(nn)=exp(-i*(kd*position(:,nn)));
    end
%     Rx = sig1*as*as'+noise*eye(N);
   w_beamforming=a/N; 
  ww(mm,:)=a/N; 
   P(mm)=w_beamforming'*Rx0*w_beamforming;
end
figure()
plot(theta_0,10*log10(abs(P)),'LineWidth',1.5);
grid on;
% xlim([-90,90]);ylim([-15,40])
xlabel('方位/(^o)');ylabel('方位谱/dB')
title('方位谱')
