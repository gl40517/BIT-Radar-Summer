%% 设置参数
c = 3e8;                 %光速
fc = 10e9;           %载频
Tp = 10e-6;             %脉冲宽度（s）
B =  10e6;               %发射信号带宽（Hz）
K = B/Tp;                %调频率
fs = 100e6;              %采样率（Hz）
R = 3000;                %初始采样距离

Nr = (2*R/c + Tp) * fs; %距离向点数
tr = (0:Nr-1)*(1/fs); %距离向时间
Na = Tp * fs;         %纯信号部分点数
% ta = (0:Na-1)*(1/fs);   %及其时间
ta = (-Na/2:Na/2-1)*(1/fs);   %及其时间
f = (0:Na-1)*fs/Na - fs/2; %及其频率

tdelay = 2*R/c;
%% 回波信号生成
% sr = rectpuls(tr-tdelay-Tp/2,Tp) .* exp(1j*pi*K*(tr-tdelay).^2) .* exp(-1*1j*2*pi*fc*tdelay); %带初始0值
sr = rectpuls(tr-tdelay-Tp/2,Tp) .* exp(1j*pi*K*(tr-tdelay-Tp/2).^2) .* exp(-1*1j*2*pi*fc*tdelay); %带初始0值
s = exp(1j*pi*K*ta.^2) .* exp(-1*1j*2*pi*fc*tdelay);  %纯信号
fft_s = fftshift(fft(s));

figure;
subplot(221)
plot(tr*c/2, real(sr));
xlabel('I(t)距离单位/m')
ylabel('幅度')
subplot(222)
plot(tr*c/2, imag(sr));
xlabel('Q(t)距离单位/m')
ylabel('幅度')
suptitle('回波信号')
subplot(223)
plot(f/1e6,20*log10(abs(fft_s)/max(abs(fft_s))));
xlabel('频率/MHz')
ylabel('幅度/dB')
subplot(224)
plot(f/1e6, phase(s));
xlabel('频率/MHz')
ylabel('弧度/rad')
%% 脉冲压缩
S_RC = fft(sr);                                     %对回波信号进行fft（距离向） 
h_RC = ifft(conj(fft(s,round(Nr))));    %匹配滤波函数    
figure
suptitle('参考信号')
subplot(221)
plot(tr*c/2, real(h_RC))
xlabel('距离单位/m')
title('滤波器实部')
subplot(222)
plot(tr*c/2, imag(h_RC))
xlabel('距离单位/m')
title('滤波器虚部')
subplot(223)
plot(f/1e6, 20*log10(abs(fftshift(fft(h_RC(2001:3000))))/max(abs( fft( h_RC(2000:2999) ) ))))
xlabel('频率/MHz')
title('滤波器频谱')
subplot(224)
plot(f/1e6, phase(h_RC(2001:3000))/pi*180)
xlabel('频率/MHz')
title('滤波器相位')

H_RC = fft(h_RC);                                   %对匹配滤波函数进行fft（距离向）
Sout_RC = S_RC.*H_RC;                              %匹配滤波（距离向）
s_RC = ifft(Sout_RC);    %距离向脉冲压缩结果
s_RC_norm = s_RC/max(s_RC);
figure()
plot(tr*c/2, 20*log10(abs(s_RC_norm)));
hold on
plot([0,4500], [-13.26 -13.26],'r')
ylim([-150,0])
xlabel('距离单位/m')
ylabel('幅度/dB')
title('显示副瓣高度的脉冲压缩后信号');
figure;
fft_s_out = fftshift(fft(s_RC_norm));
plot(20*log10(abs(fft_s_out)));
