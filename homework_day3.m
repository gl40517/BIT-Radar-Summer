% 32�������Ա�ը�ˣ���ֻ������8����������
%�������õ���m����α�����
%������ģ����������ȸ��ߺ���������Ľ���ͼ
%% 
f0 = 1e9;%��Ƶ�ź�Ƶ��
fs = 1e6;%������
Tp = 100e-6;%����ʱ��
PRT = 1e-3;
PRF = 1/PRT;
Np = 8;
Nr = PRT * fs;
N_total = Np*PRT*fs;
N_signal = Tp*fs;
tr = (0: Nr-1)/fs;            %�źŵ�ʱ��
delta_fd=10*PRF / 1000;%Ƶƫ������
fd = -10*PRF/2: delta_fd: 10*PRF/2-delta_fd;%Ƶƫ�ķ�Χ
tr_total = (0: N_total-1)/fs;%ʱ�ӵķ�Χ
tau = -PRT*Np : 1/fs: PRT*Np-2/fs;
x = zeros(1000, N_total*2-1);
for i = 1:2
    if i ==1 
        st = rectpuls(tr - Tp/2, Tp);%cw�ź�
    else
        st = (idinput(127, 'prbs'))';     %Nr�������ҲҪ��
        st = st(1:N_signal);
        st = [st, zeros(1, Nr-N_signal)];
    end
    st = repmat(st, 1, Np);
    for a=1:length(fd)
        stao=st.*exp(1j*2*pi*fd(a)*tr_total);%Ƶƫ����ź�
        x(a,:) = xcorr(st,stao);%��ԭʼ�źŸ�Ƶƫ����ź������
    end
    figure;
    mesh(tau*1e6,fd,abs(x))
    xlabel('ʱ��/us');
    ylabel('Ƶƫ/Hz');
    if i ==1
        title('���崮ģ��������άͼ')
    else
        title('�������(m����)��ģ��������άͼ')
    end
    figure;
    contour(tau*1e6,fd,abs(x));
    xlabel('ʱ��/us');
    ylabel('Ƶƫ/Hz');
        if i ==1
        title('���崮ģ�������ȸ���ͼ')
    else
        title('�������(m����)��ģ�������ȸ���ͼ')
        end
    figure;
    plot(abs(x(500,:)))
    if i ==1
        title('���崮\tau������ͼ')
    else
        title('�������(m����)��\tau������ͼ')
    end
    figure;
    plot(abs(x(:,4000)))
    if i ==1
        title('���崮������������ͼ')
    else
        title('�������(m����)�Ķ�����������ͼ')
    end
    
end