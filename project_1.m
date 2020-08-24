clear all; clc; close all;
tic;
%% ��������
load('ball_256_32_1result.mat');
fc = 77e9; %��Ƶ77GHz
PRT0 = 138e-6; %�����Ƶ��
PRF = 1 / PRT0;
fs = 10e6;%������
k = 105e6*1e6;  %б��
c = 3e8;   %����
N1 = 32;   %֡��
N2 = 256;  %������PRT��ʱ�����
N3 = 255;  %ÿ֡��PRT��
N_total = size(adcData,2); %adc�ɼ���ȫ������
t = (0:N2-1)/fs;

%% ����4ͨ������
adcData_frame = zeros();
for i = 1:4
    adcData_frame = reshape(adcData(i,:),N2*N3,N1)+adcData_frame;
end
adcData_frame = adcData_frame'/4;

%% ��ÿһ֡�źŽ���RD����
Rx = [];
final_Data = zeros(N1,N3,N2);
for i = 1:N1
    Data0 = reshape(adcData_frame(i,:),[N2,N3]); 
    Data0 = Data0';
    Data1 = padarray(Data0,[1,0],'symmetric','post');
    Data1(1,:) = [];
    Data = -Data0+Data1;    %MTI
    fft_Data0 = (fft2(Data));
    [cfar_Data,x] = cfar(abs(fft_Data0)); %����CFAR��⺯��
    Rx(i) = x;
    cfar_Data1 = cfar_Data;
    final_Data(i,:,:) = reshape(cfar_Data1,[1,N3,N2]);
end
Distance = [];
for i = 1:N1
    f0 = ((Rx(i)-17)/N2)*fs;
    Distance(i) = f0*c/2/k;
end
%% ����-֡��ͼ
figure(1);
scatter(1:N1,Distance,'*');
grid on;
title('���ھ���-֡��ͼ');
ylabel('����/m');
xlabel('֡��/��');
%% ����RDͼ��
final_Data1 = [];
Xlabel = (-128:127)/N2*fs*c/2/k;
Ylabel = (-127:127)/N3*PRF*c/2/fc;
for i = 1:N1
    final_Data1 = reshape(final_Data(i,:,:),[N3,N2]);
    figure(2);
    final_Data1_shift = fftshift(final_Data1);
    imagesc(Xlabel(100:180),Ylabel(50:200),final_Data1_shift(50:200,100:180));
    xlabel('����/m');
    ylabel('�ٶ�/(m/s)');
    title('����RDͼ��');
    pause(0.3)
end
toc;