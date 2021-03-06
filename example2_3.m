%******************example2.1*************
close all
clear all
clc

N=10;
lambda=1;  %波长
theta0=10;  %单频信号入射方向
v0=[-sind(theta0),cosd(theta0)];  %波束方向
k0=2*pi./lambda*v0;
sig=1000;
noise=1;

%*************************基阵位置******************************************
position=zeros(2,N);
position(2,:)=0;
for mm=1:N
    position(1,mm)=(mm-1)*lambda./2;   
    as(mm)=exp(-i*(k0*position(:,mm)));   %基阵响应向量
end
as=as';%列向量
Rx0 =sig*as*as'+noise*eye(N);%+ sigma_i*a_i*a_i';  %接收信号协方差矩阵
w_beamforming=inv(Rx0)*as/(as'*inv(Rx0)*as);  %MVDR波束加权向量
% w_beamforming=as/N;  %常规波束加权向量

%***************************波束扫描方向************************************
theta = (-90:0.1:90).';
ptheta=zeros(1,length(theta));
v = [-sind(theta) cosd(theta)];
k = 2*pi/lambda*v;   
for nn=1:N
    a_theta(:,nn)=exp(i*(k*position(:,nn))) ;
end
%**************************************************************************
ptheta= w_beamforming'*a_theta.';  %波束响应  pthata=加权矢量*阵列流形向量   即是有阵决定的，与该方向是否有信号无关
figure
plot(theta,20*log10(abs(ptheta)),'LineWidth',1.5);%%%%期望信号10度方向的波束图
grid on
% xlim([-90 90]);ylim([-90 10])
xlabel('波束图  方位');ylabel('波束/dB');
% title('10元阵列MVDR波束形成')



%***********方位普**********************
theta_0=-90:1:90;
figure
P=zeros(1,length(theta_0));
a=zeros(N,1);
for mm=1:length(theta_0)
    vd=[-sind(theta_0(mm)),cosd(theta_0(mm))];
    kd=2*pi./lambda*vd;
    for nn=1:N
        a(nn)=exp(-i*(kd*position(:,nn)));
    end
%     Rx = sig*as*as'+noise*eye(N);
   w_beamforming=inv(Rx0)*a/(a'*inv(Rx0)*a); 
  ww(mm,:)=inv(Rx0)*a/(a'*inv(Rx0)*a); 
   P(mm)=w_beamforming'*Rx0*w_beamforming;
end
figure(2)
plot(theta_0,10*log10(abs(P)),'LineWidth',1.5);
grid on;
xlim([-90,90]);ylim([-15,40])
xlabel('方位/(^o)');ylabel('方位谱/dB')
title('方位谱')
figure
 plot(theta_0,10*log10(ww(1:181,:).^2),'LineWidth',1.5);figure(gcf)
 grid on;
 xlabel('加权向量范数');ylabel('l0lg(||w||^2)')
% plot(theta_0,10*log10(abs(w_beamforming.^2)),'LineWidth',1.5)