% 32个串电脑爆炸了，故只仿真了8个脉冲序列
%随机码采用的是m序列伪随机码
%仿真了模糊函数、其等高线和两个轴向的截面图
%% 
f0 = 1e9;%单频信号频率
fs = 1e6;%采样率
Tp = 100e-6;%脉冲时间
PRT = 1e-3;
PRF = 1/PRT;
Np = 8;
Nr = PRT * fs;
N_total = Np*PRT*fs;
N_signal = Tp*fs;
tr = (0: Nr-1)/fs;            %信号的时间
delta_fd=10*PRF / 1000;%频偏的增量
fd = -10*PRF/2: delta_fd: 10*PRF/2-delta_fd;%频偏的范围
tr_total = (0: N_total-1)/fs;%时延的范围
tau = -PRT*Np : 1/fs: PRT*Np-2/fs;
x = zeros(1000, N_total*2-1);
for i = 1:2
    if i ==1 
        st = rectpuls(tr - Tp/2, Tp);%cw信号
    else
        st = (idinput(127, 'prbs'))';     %Nr改了这个也要改
        st = st(1:N_signal);
        st = [st, zeros(1, Nr-N_signal)];
    end
    st = repmat(st, 1, Np);
    for a=1:length(fd)
        stao=st.*exp(1j*2*pi*fd(a)*tr_total);%频偏后的信号
        x(a,:) = xcorr(st,stao);%对原始信号跟频偏后的信号作相关
    end
    figure;
    mesh(tau*1e6,fd,abs(x))
    xlabel('时延/us');
    ylabel('频偏/Hz');
    if i ==1
        title('脉冲串模糊函数三维图')
    else
        title('加随机码(m序列)的模糊函数三维图')
    end
    figure;
    contour(tau*1e6,fd,abs(x));
    xlabel('时延/us');
    ylabel('频偏/Hz');
        if i ==1
        title('脉冲串模糊函数等高线图')
    else
        title('加随机码(m序列)的模糊函数等高线图')
        end
    figure;
    plot(abs(x(500,:)))
    if i ==1
        title('脉冲串\tau轴剖面图')
    else
        title('加随机码(m序列)的\tau轴剖面图')
    end
    figure;
    plot(abs(x(:,4000)))
    if i ==1
        title('脉冲串多普勒轴剖面图')
    else
        title('加随机码(m序列)的多普勒轴剖面图')
    end
    
end