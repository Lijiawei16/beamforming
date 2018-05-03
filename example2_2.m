%*************************example2.2*************************************
close all
clear all
clc

N=10;
lambda=1;
theta0=10;
v0=[-sind(theta0),cosd(theta0)];  %波束方向
k0=2*pi./lambda*v0;
sig=1;
%*************************基阵位置******************************************
position=zeros(2,N);
position(2,:)=0;
for mm=1:N
    position(1,mm)=(mm-1)*lambda./2;   
    as(mm)=exp(-i*(k0*position(:,mm)));   %基阵响应向量
end
as=as';%列向量
Rx0 =sig*as*as' ;%+ sigma_i*a_i*a_i';  %接收信号协方差矩阵
% w_beamforming=as/(as'*as);  %常规波束加权向量
w_beamforming=as/N;  %常规波束加权向量

%***************************波束扫描方向**********************************
theta = (-90:0.1:90).';
ptheta=zeros(1,length(theta));
v = [-sind(theta) cosd(theta)];
k = 2*pi/lambda*v;   
for nn=1:N
    a_theta(:,nn)=exp(i*(k*position(:,nn))) ;
end
% a_theta = [exp(-(k*p1)*i),exp(-(k*p2)*i),exp(-(k*p3)*i),exp(-(k*p4)*i),exp(-(k*p5)*i),exp(-(k*p6)*i),exp(-(k*p7)*i),exp(-(k*p8)*i),exp(-(k*p9)*i),exp(-(k*p10)*i)];
ptheta= w_beamforming'*a_theta.';  %波束响应  pthata=加权矢量*阵列流形向量   即是有阵决定的，与该方向是否有信号无关
figure
plot(theta,20*log10(abs(ptheta)),'LineWidth',1.5);%%%%期望信号10度方向的波束图
grid on;xlabel('波束图  方位');ylabel('波束/dB');

