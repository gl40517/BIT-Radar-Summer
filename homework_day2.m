%% ���ò���
c = 3e8;                                                              %����
fc = 10e9;                                                           %��Ƶ
lambda = c/fc;
Tp = 10e-6;                                                        %�����ȣ�s��
B =  10e6;                                                          %�����źŴ���Hz��
K = B/Tp;                                                            %��Ƶ��
fs = 100e6;                                                         %�����ʣ�Hz��
Rmax = 3000;                                                     %��ʼ��������
PRT = 100e-6;                                                    %�����ظ�����
PRF = 1/PRT;                                                      %�����ظ�Ƶ�ʣ�Hz��
Nr = round( (2*Rmax/c + Tp) * fs );                                       %���������
tr = (0:Nr-1)*(1/fs);                                             %������ʱ��
N_signal = Tp * fs;                                                       %���źŲ��ֵ���
t_signal = (-N_signal/2:N_signal/2-1)*(1/fs);                                 %����ʱ��
v = 60;
Na = 64;
ta = (0:Na-1)/(PRF);       %��λ��ʱ��
x_target = Rmax - ta*v;
%% �����ٶ�ģ���ͷֱ���
delta_f = PRF / Na;
delta_v = lambda / 2 * delta_f;
fprintf('�ٶȷֱ���Ϊ%fm/s\n', delta_v);
max_v = c* PRF/2 / 2 / fc;
fprintf('���ģ���ٶ�Ϊ����%fm/s\n', max_v);
%% �ز��ź�����
sr = zeros(Na, Nr);

for n = 1:Na
    tdelay = 2*x_target(n) / c;
    sr(n,:) = rectpuls(tr-tdelay-Tp/2,Tp) .* exp(1i*pi*K*(tr-tdelay-Tp/2).^2) .* exp(-1i*4*pi*x_target(n)/lambda); 
end
%% �ο��ź�����
% mesh(real(sr))
s_signal = exp(1j*pi*K*t_signal.^2) .* exp(-1*1j*2*pi*fc*tdelay);  %���ź�
h_RC = repmat(s_signal, Na,1);
H_RC = conj(fft( h_RC,  Nr, 2 ));    %ƥ���˲����� 
%%��������գ�PD������
S_RC = fft(sr,[],2);%�Իز��źŽ���fft�������� 
S_out = S_RC .* H_RC;
s_out = ifft(S_out, [], 2);
S_out_1 = fftshift(fft(s_out), 1);
f = (-Na/2:Na/2-1)/Na*PRF;
figure
[X,Y] = meshgrid(tr*c/2, f/1000);
mesh(X, Y, abs(S_out_1))
xlabel('���뵥λ/m')
ylabel('Ƶ��/kHz')
title('����-����������άͼ')
figure;
imagesc(tr*c/2, f/1000, abs(S_out_1))
axis xy
xlabel('���뵥λ/m')
ylabel('Ƶ��/kHz')
title('����-���������άͼ')
figure
plot(f/1000, 20*log10(abs(S_out_1(:,2000))))
xlabel('Ƶ��/kHz')
title('ѡȡ����=3km���Ķ����������ͼ')

%% ���ķ������ٶȼ����루�в���ͼ��룩
sr = zeros(Na, Nr, 4);

for n = 1:Na
    tdelay = 2*x_target(n) / c;
    sr(n, :, 1) = rectpuls(tr-tdelay-Tp/2,Tp) .* exp(1i*pi*K*(tr-tdelay-Tp/2).^2) .* exp(-1i*4*pi*x_target(n)/lambda); 
    sr(n, :, 2) = rectpuls(tr-tdelay-Tp/2,Tp) .* exp(1i*pi*K*(tr-tdelay-Tp/2).^2) .* exp(-1i*4*pi*x_target(n)/lambda); 
    sr(n, :, 3) = rectpuls(tr-tdelay-Tp/2,Tp) .* exp(1i*pi*K*(tr-tdelay-Tp/2).^2) .* exp(-1i*4*pi*x_target(n)/lambda);
    sr(n, :, 4) = rectpuls(tr-tdelay-Tp/2,Tp) .* exp(1i*pi*K*(tr-tdelay-Tp/2).^2) .* exp(-1i*4*pi*x_target(n)/lambda);
end
sr(:, :, 2) = awgn(sr(:, :, 1), 0);
sr(:, :, 3) = awgn(sr(:, :, 1), 10);
sr(:, :, 4) = awgn(sr(:, :, 1), 20);
%% �ο��ź�����
s_signal = exp(1j*pi*K*t_signal.^2) .* exp(-1*1j*2*pi*fc*tdelay);  %���ź�
h_RC = repmat(s_signal, Na,1);
H_RC = conj(fft( h_RC,  Nr, 2 ));    %ƥ���˲����� 
%% ���ķ����������ٶȣ��в��㣩
 f = (-Na/2:Na/2-1)/Na*PRF;
 N_process = 8;                                         %ѡȡ8������м�Ȩ����
 Final_distance = zeros(1, 4);
 Final_velocity = zeros(1, 4);
 Tr = (0:4*Nr-1) * (1/(4*fs));                      %������ʱ���ṹ�죨����ԭ����4����
 f = (0: 8*Na-1) / (8*Na) * PRF;                 %�����Ķ�����Ƶ���ṹ�죨����ԭ����8����
 ready_for_distance = zeros(1, 4*Nr);
 ready_for_velocity = zeros(8*Na, 1);
 
for i=1:4
    S_RC = sr(:, :, i);
    S_RC = fft(S_RC, [], 2);                     
    S_out = S_RC .* H_RC;
    s_out = ifft(S_out, [], 2);
    
    max_value = max(max(abs(s_out)));
    [velocity_index_before, distance_index_before] = find(abs(s_out) == max_value);
    ready_for_distance = [S_out(velocity_index_before, :), zeros(1,3*Nr)];
    ready_for_distance = ifft(ready_for_distance, [], 2);
    [~, distance_index] = max(ready_for_distance);
    ready_for_velocity = [s_out(:, distance_index_before); zeros(7*Na, 1)];
    ready_for_velocity = fft(ready_for_velocity);
    [~, velocity_index] = max(ready_for_velocity);

    count_distance = 0;
    count_distance_magnitude = 0;
    for j = -N_process/2: N_process/2 - 1
        R = Tr(distance_index + j)*c / 2;
        count_distance = count_distance + abs(ready_for_distance(distance_index+j)) * R;
        count_distance_magnitude = count_distance_magnitude + abs(ready_for_distance(distance_index+j));
    end
    distance = count_distance / count_distance_magnitude;
    Final_distance(i) = distance;
    %-------------------�������ķ�����------------------
    count_velocity = 0;
    count_velocity_magnitude = 0;
    for j = -N_process/2: N_process/2 - 1
        V = f(velocity_index + j)*lambda/2;
        count_velocity = count_velocity + abs(ready_for_velocity(velocity_index+j)) * V;
        count_velocity_magnitude = count_velocity_magnitude + abs(ready_for_velocity(velocity_index+j));
    end
    velocity = count_velocity / count_velocity_magnitude;
    Final_velocity(i) = velocity;
    
end
Final_distance
Final_velocity


