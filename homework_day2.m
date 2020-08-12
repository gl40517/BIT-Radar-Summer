%% 设置参数
c = 3e8;                                                              %光速
fc = 10e9;                                                           %载频
lambda = c/fc;
Tp = 10e-6;                                                        %脉冲宽度（s）
B =  10e6;                                                          %发射信号带宽（Hz）
K = B/Tp;                                                            %调频率
fs = 100e6;                                                         %采样率（Hz）
Rmax = 3000;                                                     %初始采样距离
PRT = 100e-6;                                                    %脉冲重复周期
PRF = 1/PRT;                                                      %脉冲重复频率（Hz）
Nr = round( (2*Rmax/c + Tp) * fs );                                       %距离向点数
tr = (0:Nr-1)*(1/fs);                                             %距离向时间
N_signal = Tp * fs;                                                       %纯信号部分点数
t_signal = (-N_signal/2:N_signal/2-1)*(1/fs);                                 %及其时间
v = 60;
Na = 64;
ta = (0:Na-1)/(PRF);       %方位向时间
x_target = Rmax - ta*v;
%% 计算速度模糊和分辨率
delta_f = PRF / Na;
delta_v = lambda / 2 * delta_f;
fprintf('速度分辨率为%fm/s\n', delta_v);
max_v = c* PRF/2 / 2 / fc;
fprintf('最大不模糊速度为正负%fm/s\n', max_v);
%% 回波信号生成
sr = zeros(Na, Nr);

for n = 1:Na
    tdelay = 2*x_target(n) / c;
    sr(n,:) = rectpuls(tr-tdelay-Tp/2,Tp) .* exp(1i*pi*K*(tr-tdelay-Tp/2).^2) .* exp(-1i*4*pi*x_target(n)/lambda); 
end
%% 参考信号生成
% mesh(real(sr))
s_signal = exp(1j*pi*K*t_signal.^2) .* exp(-1*1j*2*pi*fc*tdelay);  %纯信号
h_RC = repmat(s_signal, Na,1);
H_RC = conj(fft( h_RC,  Nr, 2 ));    %匹配滤波函数 
%%脉冲多普勒（PD）处理
S_RC = fft(sr,[],2);%对回波信号进行fft（距离向） 
S_out = S_RC .* H_RC;
s_out = ifft(S_out, [], 2);
S_out_1 = fftshift(fft(s_out), 1);
f = (-Na/2:Na/2-1)/Na*PRF;
figure
[X,Y] = meshgrid(tr*c/2, f/1000);
mesh(X, Y, abs(S_out_1))
xlabel('距离单位/m')
ylabel('频率/kHz')
title('距离-多普勒域三维图')
figure;
imagesc(tr*c/2, f/1000, abs(S_out_1))
axis xy
xlabel('距离单位/m')
ylabel('频率/kHz')
title('距离-多普勒域二维图')
figure
plot(f/1000, 20*log10(abs(S_out_1(:,2000))))
xlabel('频率/kHz')
title('选取距离=3km处的多普勒域截面图')

%% 形心法计算速度及距离（有补零和加噪）
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
%% 参考信号生成
s_signal = exp(1j*pi*K*t_signal.^2) .* exp(-1*1j*2*pi*fc*tdelay);  %纯信号
h_RC = repmat(s_signal, Na,1);
H_RC = conj(fft( h_RC,  Nr, 2 ));    %匹配滤波函数 
%% 形心法计算距离和速度（有补零）
 f = (-Na/2:Na/2-1)/Na*PRF;
 N_process = 8;                                         %选取8个点进行加权计算
 Final_distance = zeros(1, 4);
 Final_velocity = zeros(1, 4);
 Tr = (0:4*Nr-1) * (1/(4*fs));                      %补零后的时间轴构造（补到原长的4倍）
 f = (0: 8*Na-1) / (8*Na) * PRF;                 %补零后的多普勒频率轴构造（补到原长的8倍）
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
    %-------------------距离形心法计算------------------
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


