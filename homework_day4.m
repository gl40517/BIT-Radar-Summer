%% 
%运行时间约25秒
%利用RD算法成像
%% 设置参数
tic
c = 3e8;                 %光速
vs = 100;               %卫星速度 100 m/s
fc = 15e9;           %载频
lambda = c/fc;           %载波波长（m
Tp = 0.5e-6;             %脉冲宽度（s）
B =  200e6;               %发射信号带宽（Hz）
K = B/Tp;                %调频率
fs = 240e6;              %采样率（Hz）
PRF = 500;               %脉冲重复频率（Hz）
PRT = 1/PRF;             %脉冲重复周期
A = 0.02;               %方位向波束宽度（rad）
Rmin = 9600;   %初始采样距离
x_target = [ 0 ]; %目标坐标
r_target = [10000];
Ls = A*Rmin;             %合成孔径长度

Nr =round(2*fs*(r_target(1)+400 - Rmin)/c + Tp * fs);                   %距离向点数
tr_rd = 2*Rmin/c+(0:Nr-1)*(1/(fs)); %距离向时间
Na = ceil(Ls*PRF/vs);                  %方位向点数
ta_rd = (-Na/2:Na/2-1)/(PRF);       %方位向时间
%% 计算参数
tar_num = length(x_target);
x_radar = ta_rd*vs;     %雷达方位向坐标
r_radar = 0;         %雷达距离向坐标
%% 理论计算
resR = 0.886*c/(2*B);                                           % 距离向理论分辨率                                        
resA = lambda*Rmin/2/Ls;                                                   % 方位向理论分辨率
disp(['距离向理论分辨率为(带0.886的系数)：',num2str(resR),' m'])
disp(['方位向理论分辨率为：',num2str(resA),'m'])
%% 回波信号
sr = zeros(Na,Nr);
Rn1 = zeros(Na,tar_num); %存放斜距

for n = 1:Na
    for k = 1:tar_num
      Rn1(n,k) = sqrt( ( x_radar(n)-x_target(k) ).^2+( r_radar-r_target(k) ).^2 );
         tdelay = 2*(Rn1(n,k))/c;
         s1 = rectpuls(tr_rd-tdelay-Tp/2,Tp);
         s2 = rectpuls(x_radar(n)-x_target(k),Ls);
         s3 = exp(1i*pi*K*(tr_rd-tdelay-Tp/2).^2);
         s4 = exp(-1i*4*pi*Rn1(n,k)/lambda);
         sr(n,:) = sr(n,:)+s1.*s2.*s3.*s4;
     end
end

figure;
imagesc(abs(sr));
title('RD算法点目标回波信号');
%% 距离徙动计算
Rf = zeros(Na,Nr);   %Rf为距离徙动计算值
for n = 1:Na
    for t = 1:Nr
        Rf(n,t) = (1/8) * ((lambda/vs).^2).*(Rmin).*(ta_rd(n)/Na*(PRF)^2).^2; %距离徙动计算
        
    end
end
figure;
imagesc(abs(Rf));
title('RD算法距离徙动大小');
%% 距离向脉冲压缩
S_RC = fft(sr,[],2);                                     %对回波信号进行fft（距离向）
h_RC = rectpuls(tr_rd-2*Rmin/c-Tp/2,Tp).*exp(1i*K*pi*(tr_rd-2*Rmin/c-Tp/2).^2);  %匹配滤波函数
h_RC = repmat(h_RC,Na,1);
H_RC = fft(h_RC,[],2);                                   %对匹配滤波函数进行fft（距离向）
Sout_RC = S_RC.*conj(H_RC);                              %匹配滤波（距离向）
s_RC = ifft(Sout_RC,[],2);                               %距离向脉冲压缩结果
figure;
imagesc(abs(s_RC));
title('RD算法距离向脉冲压缩后信号');
%% 变换到方位频域
S_AC = fftshift(fft(s_RC),1);            %变换到方位频域
figure;
imagesc(abs(S_AC));
title('RD算法变换到方位频域后的信号');
%% 距离徙动校正1

sj = zeros(Na,Nr);                       %sj为距离徙动校正后信号
cha = zeros(Na,Nr);                      %sinc插值核
Rj = 2*Rf*fs/c;                          %计算得到当前Rf对应的点数，便于sinc插值
Rjz = fix(Rj);                           %整数部分
Rjd = Rj - Rjz;                          %小数部分
for n = 1:Na
    
    for t = 1:Nr
        for i = -4:3                     %进行8点sinc插值计算
        cha(n,t) = sinc(Rjd(n,t)-i);     %构造sinc插值核
        if t+Rjz(n,t)+i<Nr+1 && t+Rjz(n,t)+i>0
            sj(n,t) = S_AC(n,t+Rjz(n,t)+i).*cha(n,t)+sj(n,t);
        else
            sj(n,t) = S_AC(n,t);
        end
        end
    end
end

figure(5)
imagesc(abs(sj));
title('RD算法距离徙动校正结果');        
    
%% 距离徙动校正2
% sj = zeros(Na,Nr);                       %sj为距离徙动校正后信号
% Rj_used_to_correct = zeros(Na,Nr);
% Rj = 2*Rf*fs/c;            %计算得到当前Rf对应的点数，便于sinc插值
% k = 1:Nr;
% S_RAC = fft(S_AC,[],2);  
% for i = 1:Na
%     Rj_used_to_correct(i,:) = exp(1j * 2*pi/Nr .* k .* Rj(i,1));
%     sj(i,:) = S_RAC(i,:) .* Rj_used_to_correct(i,:);
%     
% end
% sj = ifft(sj ,[],2);  
% figure;
% imagesc(abs(sj));
% title('距离徙动校正结果');        
%% 方位向脉冲压缩
h_AC = zeros(Na,Nr);                   %h_AC为方位向匹配滤波函数                          
Rline = Rmin + (0:Nr-1)*c/(2*fs);      %不同的距离门
for t = 1:Nr
    ka = 2*vs^2/(lambda.*Rline(t));    %调频斜率
    h_AC(:,t) = rectpuls(vs*ta_rd,Ls).*exp(-1i*pi*ka*ta_rd.^2);
end
H_AC = fftshift(fft(h_AC),1);
Sout_AC = sj.*conj(H_AC);              %匹配滤波（方位向）
%% 方位向逆fft
s_RAC = fftshift(ifft(fftshift(Sout_AC,1)),1);  %s_RAC为RD算法后的结果
figure;
% imagesc(c*tr/2, x_radar, abs(s_RAC));
imagesc(tr_rd*c/2000, ta_rd*vs, abs(s_RAC));
title('RD算法处理后结果');
xlabel('距离向距离/km')
ylabel('方位向距离/m')
%% 升采样观察点扩展函数（8倍升采样）
interp = 8; %升采样倍数
figure;
lx = 101;
ly = 201;
lx1 = interp*lx;
lx2 = lx1-lx;
tr_cut = tr_rd(600:700);
ta_cut = ta_rd(400:600);
AziComp = s_RAC(400:600, 600:700);
fAziComp = fft(AziComp,[],2);
fAziComp1 = [ fAziComp(:,1:round(lx/2)) , zeros(ly,lx2) , fAziComp(:,round(lx/2)+1:lx)]; %频域补零
AziComp1 = ifft(fAziComp1,[],2);  %距离向频域补零后的数据
ly1 = interp*ly;
ly2 = ly1-ly;
fAziComp1 = fft(AziComp1);
fAziComp2 = [ fAziComp1(1:round(ly/2),:) ; zeros(ly2,lx1) ; fAziComp1(round(ly/2)+1:ly,:)]; %频域补零
AziComp2 = ifft(fAziComp2);  %方位向频域补零后的数据

tr_up = (tr_cut(1)*(interp*fs): tr_cut(101)*(interp*fs))/(interp*fs);
ta_up = (ta_cut(1)*(interp*PRF) : ta_cut(201)*(interp*PRF))/(interp*PRF);
imagesc(tr_up*c/2000, ta_up*vs, abs(AziComp2))
xlabel('距离向升采样/km')
ylabel('方位向升采样/m')
title('RD算法点扩展函数')

%% BP算法仿真
%% 参数设置
tr_bp = (0:Nr-1)*(1/fs);
% 目标坐标：单点
distance_target = [r_target(1)] ;          
position_target = [0];
height_target = [0];

tar_num = length(x_target);                                                % 目标个数

% 雷达初始坐标
distance_radar = 0;        
position_radar = 0;
height_radar = 0;
echo = zeros(Na,Nr);
for n = 1:Na
    for tar = 1:tar_num
        Rn = sqrt((distance_radar-distance_target(tar)).^2 + (position_radar+ta_rd(n)*vs-position_target(tar)).^2 ...
            + (height_radar-height_target(tar)).^2 );
        tdelay = 2*(Rn-Rmin)/c;
        r1 = rectpuls(tr_bp-tdelay-Tp/2,Tp);
        r2 = exp(1i*pi*K*(tr_bp-tdelay-Tp/2).^2);
        r3 = exp(-1i*4*pi*Rn/lambda);
        echo(n,:) = echo(n,:)+r1.*r2.*r3;
       
    end
end
figure;
imagesc(tr_bp*c/2+Rmin, ta_rd, abs(echo));
ylabel('方位向时刻（s）')
xlabel('距离向时刻（s）')
title('BP算法回波信号幅度值')
%% 距离向脉冲压缩
echo_cmp_temp = fft(echo,[],2);
match_filter = rectpuls(tr_bp-Tp/2,Tp).*exp(1i*pi*K*(tr_bp-Tp/2).^2);
match_filter = repmat(match_filter,Na,1);
match_filter_f = fft(match_filter,[],2);
echo_cmp_f = echo_cmp_temp.*conj(match_filter_f);
interp = 8;                                                                % 升采样倍数
Nr_new = interp*Nr;                                                        % 补零后距离向点数
Nr_res = Nr_new-Nr;
echo_cmp_f_add = [echo_cmp_f(:,1:Nr/2),zeros(Na,Nr_res),echo_cmp_f(:,Nr/2:Nr)]; %频域补零
echo_cmp = ifft(echo_cmp_f_add,[],2);
figure;
imagesc(abs(echo_cmp));
xlabel('方位向点数')
ylabel('距离向点数')
title('BP算法距离向脉冲压缩结果')

%% 成像区域确定（输出图像）
deltaR = (1/fs)*c/2;
deltaA = (1/PRF)*vs;
% 确定场景大小
xleft = distance_target-100;
xright = distance_target+100;
ydown = position_target-50;
yup = position_target+50 ;
%x、y分别为横坐标和纵坐标
x_line = xleft : deltaR : xright;
y_line = ydown : deltaA : yup;
lx = length(x_line);
ly = length(y_line);
x = ones(ly,1)*x_line;
y = y_line'*ones(1,lx);
%% BP成像
R = zeros(ly,lx);                                                          % 存储斜距历程                                                
BP_Image = zeros(ly,lx);                                                   % 存储BP成像结果
dtr = c/(2*fs*interp);                                                     % 补零后的距离向采样间隔
for index = 1:Na
    yradar = position_radar+vs*ta_rd(index);
    R = sqrt((x-distance_radar).^2 + (y-yradar).^2);                         % 计算到各像素点的斜距历程
    for ii = 1:ly
        for jj = 1:lx
            r = R(ii,jj);          
            flag = round((r-Rmin)/dtr)+1; 
            if flag > 0 && flag <= Nr_new
               BP_Image(ii,jj) = BP_Image(ii,jj) +echo_cmp(index,flag).*exp(1j*4*pi*r/lambda);
      
            else
            end
        end
    end
end    
figure;
imagesc(x_line,y_line,abs(BP_Image));
xlabel('横向距离（m）')
ylabel('纵向距离（m）')
title('BP成像结果')
%% 二维零阶插值（补）
% 首先进行补零，保证频谱在原本的位置，进行保相处理
BP_Image_add = zeros(ly,lx);                                         
for ii = 1:ly
    for jj = 1:lx
        BP_Image_add(ii,jj) = BP_Image(ii,jj).*exp(-1i*4*pi*R(ii,jj)/lambda);
    end
end
%补零操作 
interp_new = 8;    
% 对距离向进行频域补零的操作
lx_new = interp_new*lx;
lx_res = lx_new-lx;
fBP_Image_add = fft(BP_Image_add,[],2);
fBP_Image_add1 = [ fBP_Image_add(:,1:round(lx/2)) , zeros(ly,lx_res) , fBP_Image_add(:,round(lx/2)+1:lx)]; %频域补零
BP_Image_add1 = ifft(fBP_Image_add1,[],2);  %频域补零后的数据
% 对方位向进行频域补零的操作
ly_new = interp_new*ly;
ly_res = ly_new-ly;
fBP_Image_add1 = fft(BP_Image_add1);
fBP_Image_add2 = [ fBP_Image_add1(1:round(ly/2),:) ; zeros(ly_res,lx_new) ; fBP_Image_add1(round(ly/2)+1:ly,:)]; %频域补零
BP_Image_final = ifft(fBP_Image_add2);  %方位向频域补零后的数据
figure;
imagesc(abs(BP_Image_final)) 
xlabel('横向点数')
ylabel('纵向点数')
title('BP成像结果（升采样）')
figure;
mesh(abs(BP_Image_final));
xlabel('横向点数')
ylabel('纵向点数')
title('BP算法点扩展函数')
toc