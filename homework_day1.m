%% ���ò���
c = 3e8;                 %����
fc = 10e9;           %��Ƶ
Tp = 10e-6;             %�����ȣ�s��
B =  10e6;               %�����źŴ���Hz��
K = B/Tp;                %��Ƶ��
fs = 100e6;              %�����ʣ�Hz��
R = 3000;                %��ʼ��������

Nr = (2*R/c + Tp) * fs; %���������
tr = (0:Nr-1)*(1/fs); %������ʱ��
Na = Tp * fs;         %���źŲ��ֵ���
% ta = (0:Na-1)*(1/fs);   %����ʱ��
ta = (-Na/2:Na/2-1)*(1/fs);   %����ʱ��
f = (0:Na-1)*fs/Na - fs/2; %����Ƶ��

tdelay = 2*R/c;
%% �ز��ź�����
% sr = rectpuls(tr-tdelay-Tp/2,Tp) .* exp(1j*pi*K*(tr-tdelay).^2) .* exp(-1*1j*2*pi*fc*tdelay); %����ʼ0ֵ
sr = rectpuls(tr-tdelay-Tp/2,Tp) .* exp(1j*pi*K*(tr-tdelay-Tp/2).^2) .* exp(-1*1j*2*pi*fc*tdelay); %����ʼ0ֵ
s = exp(1j*pi*K*ta.^2) .* exp(-1*1j*2*pi*fc*tdelay);  %���ź�
fft_s = fftshift(fft(s));

figure;
subplot(221)
plot(tr*c/2, real(sr));
xlabel('I(t)���뵥λ/m')
ylabel('����')
subplot(222)
plot(tr*c/2, imag(sr));
xlabel('Q(t)���뵥λ/m')
ylabel('����')
suptitle('�ز��ź�')
subplot(223)
plot(f/1e6,20*log10(abs(fft_s)/max(abs(fft_s))));
xlabel('Ƶ��/MHz')
ylabel('����/dB')
subplot(224)
plot(f/1e6, phase(s));
xlabel('Ƶ��/MHz')
ylabel('����/rad')
%% ����ѹ��
S_RC = fft(sr);                                     %�Իز��źŽ���fft�������� 
h_RC = ifft(conj(fft(s,round(Nr))));    %ƥ���˲�����    
figure
suptitle('�ο��ź�')
subplot(221)
plot(tr*c/2, real(h_RC))
xlabel('���뵥λ/m')
title('�˲���ʵ��')
subplot(222)
plot(tr*c/2, imag(h_RC))
xlabel('���뵥λ/m')
title('�˲����鲿')
subplot(223)
plot(f/1e6, 20*log10(abs(fftshift(fft(h_RC(2001:3000))))/max(abs( fft( h_RC(2000:2999) ) ))))
xlabel('Ƶ��/MHz')
title('�˲���Ƶ��')
subplot(224)
plot(f/1e6, phase(h_RC(2001:3000))/pi*180)
xlabel('Ƶ��/MHz')
title('�˲�����λ')

H_RC = fft(h_RC);                                   %��ƥ���˲���������fft��������
Sout_RC = S_RC.*H_RC;                              %ƥ���˲���������
s_RC = ifft(Sout_RC);    %����������ѹ�����
s_RC_norm = s_RC/max(s_RC);
figure()
plot(tr*c/2, 20*log10(abs(s_RC_norm)));
hold on
plot([0,4500], [-13.26 -13.26],'r')
ylim([-150,0])
xlabel('���뵥λ/m')
ylabel('����/dB')
title('��ʾ����߶ȵ�����ѹ�����ź�');
figure;
fft_s_out = fftshift(fft(s_RC_norm));
plot(20*log10(abs(fft_s_out)));
