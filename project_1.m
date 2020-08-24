clear all; clc; close all;
%% 参数设置
load('ball_256_32_1result.mat');
fc = 77e9; %载频77GHz
PRT0 = 138e-6; %脉冲的频率
PRF = 1 / PRT0;
fs = 10e6;%采样率
k = 105e6*1e6;  %斜率
c = 3e8;   %光速
N1 = 32;   %帧数
N2 = 256;  %单周期PRT快时间点数
N3 = 255;  %每帧的PRT数
N_total = size(adcData,2); %adc采集的全部数据
t = (0:N2-1)/fs;

%% 单列每一帧的数据
adcData_frame = zeros();
for i = 1:4
    adcData_frame = reshape(adcData(i,:),N2*N3,N1)+adcData_frame;
end
adcData_frame = adcData_frame';

%% 对每一帧信号进行RD处理
% Data = zeros();
Rx = zeros(1,32);
final_Data = zeros(N1,N3,N2);
for i = 1:N1
    Data0 = reshape(adcData_frame(i,:),[N2,N3]); 
    Data0 = Data0';
    Data1 = padarray(Data0,[1,0],'symmetric','post');
    Data1(1,:) = [];
    Data = -Data0+Data1;    %MTI
    fft_Data0 = (fft2(Data));
    [cfar_Data,x] = cfar(abs(fft_Data0)); %调用CFAR检测函数
    Rx(i) = x;
    cfar_Data1 = cfar_Data;
    final_Data(i,:,:) = reshape(cfar_Data1,[1,N3,N2]);
end
Distance = zeros(1,32);
for i = 1:N1
    f0 = ((Rx(i)-16)/N2)*fs;
    Distance(i) = f0*c/2/k;
end
%% 距离-帧数图
figure(1);
scatter(1:N1,Distance,'*');
grid on;
title('单摆距离-帧数图');
ylabel('距离/m');
xlabel('帧数/个');
%% 单摆RD图像
final_Data1 = [];
Xlabel = (-128:127)/N2*fs*c/2/k;
Ylabel = (-127:127)/N3*PRF*c/2/fc;
for i = 1:N1
    final_Data1 = reshape(final_Data(i,:,:),[N3,N2]);
    figure(2);
    imagesc(Xlabel,Ylabel,fftshift(final_Data1));
    xlabel('距离/m');
    ylabel('速度/(m/s)');
    title('单摆RD图像');
    pause(0.5);
end
