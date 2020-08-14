%% 
%����ʱ��Լ25��
%����RD�㷨����
%% ���ò���
tic
c = 3e8;                 %����
vs = 100;               %�����ٶ� 100 m/s
fc = 15e9;           %��Ƶ
lambda = c/fc;           %�ز�������m
Tp = 0.5e-6;             %�����ȣ�s��
B =  200e6;               %�����źŴ���Hz��
K = B/Tp;                %��Ƶ��
fs = 240e6;              %�����ʣ�Hz��
PRF = 500;               %�����ظ�Ƶ�ʣ�Hz��
PRT = 1/PRF;             %�����ظ�����
A = 0.02;               %��λ������ȣ�rad��
Rmin = 9600;   %��ʼ��������
x_target = [ 0 ]; %Ŀ������
r_target = [10000];
Ls = A*Rmin;             %�ϳɿ׾�����

Nr =round(2*fs*(r_target(1)+400 - Rmin)/c + Tp * fs);                   %���������
tr_rd = 2*Rmin/c+(0:Nr-1)*(1/(fs)); %������ʱ��
Na = ceil(Ls*PRF/vs);                  %��λ�����
ta_rd = (-Na/2:Na/2-1)/(PRF);       %��λ��ʱ��
%% �������
tar_num = length(x_target);
x_radar = ta_rd*vs;     %�״﷽λ������
r_radar = 0;         %�״����������
%% ���ۼ���
resR = 0.886*c/(2*B);                                           % ���������۷ֱ���                                        
resA = lambda*Rmin/2/Ls;                                                   % ��λ�����۷ֱ���
disp(['���������۷ֱ���Ϊ(��0.886��ϵ��)��',num2str(resR),' m'])
disp(['��λ�����۷ֱ���Ϊ��',num2str(resA),'m'])
%% �ز��ź�
sr = zeros(Na,Nr);
Rn1 = zeros(Na,tar_num); %���б��

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
title('RD�㷨��Ŀ��ز��ź�');
%% �����㶯����
Rf = zeros(Na,Nr);   %RfΪ�����㶯����ֵ
for n = 1:Na
    for t = 1:Nr
        Rf(n,t) = (1/8) * ((lambda/vs).^2).*(Rmin).*(ta_rd(n)/Na*(PRF)^2).^2; %�����㶯����
        
    end
end
figure;
imagesc(abs(Rf));
title('RD�㷨�����㶯��С');
%% ����������ѹ��
S_RC = fft(sr,[],2);                                     %�Իز��źŽ���fft��������
h_RC = rectpuls(tr_rd-2*Rmin/c-Tp/2,Tp).*exp(1i*K*pi*(tr_rd-2*Rmin/c-Tp/2).^2);  %ƥ���˲�����
h_RC = repmat(h_RC,Na,1);
H_RC = fft(h_RC,[],2);                                   %��ƥ���˲���������fft��������
Sout_RC = S_RC.*conj(H_RC);                              %ƥ���˲���������
s_RC = ifft(Sout_RC,[],2);                               %����������ѹ�����
figure;
imagesc(abs(s_RC));
title('RD�㷨����������ѹ�����ź�');
%% �任����λƵ��
S_AC = fftshift(fft(s_RC),1);            %�任����λƵ��
figure;
imagesc(abs(S_AC));
title('RD�㷨�任����λƵ�����ź�');
%% �����㶯У��1

sj = zeros(Na,Nr);                       %sjΪ�����㶯У�����ź�
cha = zeros(Na,Nr);                      %sinc��ֵ��
Rj = 2*Rf*fs/c;                          %����õ���ǰRf��Ӧ�ĵ���������sinc��ֵ
Rjz = fix(Rj);                           %��������
Rjd = Rj - Rjz;                          %С������
for n = 1:Na
    
    for t = 1:Nr
        for i = -4:3                     %����8��sinc��ֵ����
        cha(n,t) = sinc(Rjd(n,t)-i);     %����sinc��ֵ��
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
title('RD�㷨�����㶯У�����');        
    
%% �����㶯У��2
% sj = zeros(Na,Nr);                       %sjΪ�����㶯У�����ź�
% Rj_used_to_correct = zeros(Na,Nr);
% Rj = 2*Rf*fs/c;            %����õ���ǰRf��Ӧ�ĵ���������sinc��ֵ
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
% title('�����㶯У�����');        
%% ��λ������ѹ��
h_AC = zeros(Na,Nr);                   %h_ACΪ��λ��ƥ���˲�����                          
Rline = Rmin + (0:Nr-1)*c/(2*fs);      %��ͬ�ľ�����
for t = 1:Nr
    ka = 2*vs^2/(lambda.*Rline(t));    %��Ƶб��
    h_AC(:,t) = rectpuls(vs*ta_rd,Ls).*exp(-1i*pi*ka*ta_rd.^2);
end
H_AC = fftshift(fft(h_AC),1);
Sout_AC = sj.*conj(H_AC);              %ƥ���˲�����λ��
%% ��λ����fft
s_RAC = fftshift(ifft(fftshift(Sout_AC,1)),1);  %s_RACΪRD�㷨��Ľ��
figure;
% imagesc(c*tr/2, x_radar, abs(s_RAC));
imagesc(tr_rd*c/2000, ta_rd*vs, abs(s_RAC));
title('RD�㷨�������');
xlabel('���������/km')
ylabel('��λ�����/m')
%% �������۲����չ������8����������
interp = 8; %����������
figure;
lx = 101;
ly = 201;
lx1 = interp*lx;
lx2 = lx1-lx;
tr_cut = tr_rd(600:700);
ta_cut = ta_rd(400:600);
AziComp = s_RAC(400:600, 600:700);
fAziComp = fft(AziComp,[],2);
fAziComp1 = [ fAziComp(:,1:round(lx/2)) , zeros(ly,lx2) , fAziComp(:,round(lx/2)+1:lx)]; %Ƶ����
AziComp1 = ifft(fAziComp1,[],2);  %������Ƶ����������
ly1 = interp*ly;
ly2 = ly1-ly;
fAziComp1 = fft(AziComp1);
fAziComp2 = [ fAziComp1(1:round(ly/2),:) ; zeros(ly2,lx1) ; fAziComp1(round(ly/2)+1:ly,:)]; %Ƶ����
AziComp2 = ifft(fAziComp2);  %��λ��Ƶ����������

tr_up = (tr_cut(1)*(interp*fs): tr_cut(101)*(interp*fs))/(interp*fs);
ta_up = (ta_cut(1)*(interp*PRF) : ta_cut(201)*(interp*PRF))/(interp*PRF);
imagesc(tr_up*c/2000, ta_up*vs, abs(AziComp2))
xlabel('������������/km')
ylabel('��λ��������/m')
title('RD�㷨����չ����')

%% BP�㷨����
%% ��������
tr_bp = (0:Nr-1)*(1/fs);
% Ŀ�����꣺����
distance_target = [r_target(1)] ;          
position_target = [0];
height_target = [0];

tar_num = length(x_target);                                                % Ŀ�����

% �״��ʼ����
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
ylabel('��λ��ʱ�̣�s��')
xlabel('������ʱ�̣�s��')
title('BP�㷨�ز��źŷ���ֵ')
%% ����������ѹ��
echo_cmp_temp = fft(echo,[],2);
match_filter = rectpuls(tr_bp-Tp/2,Tp).*exp(1i*pi*K*(tr_bp-Tp/2).^2);
match_filter = repmat(match_filter,Na,1);
match_filter_f = fft(match_filter,[],2);
echo_cmp_f = echo_cmp_temp.*conj(match_filter_f);
interp = 8;                                                                % ����������
Nr_new = interp*Nr;                                                        % �������������
Nr_res = Nr_new-Nr;
echo_cmp_f_add = [echo_cmp_f(:,1:Nr/2),zeros(Na,Nr_res),echo_cmp_f(:,Nr/2:Nr)]; %Ƶ����
echo_cmp = ifft(echo_cmp_f_add,[],2);
figure;
imagesc(abs(echo_cmp));
xlabel('��λ�����')
ylabel('���������')
title('BP�㷨����������ѹ�����')

%% ��������ȷ�������ͼ��
deltaR = (1/fs)*c/2;
deltaA = (1/PRF)*vs;
% ȷ��������С
xleft = distance_target-100;
xright = distance_target+100;
ydown = position_target-50;
yup = position_target+50 ;
%x��y�ֱ�Ϊ�������������
x_line = xleft : deltaR : xright;
y_line = ydown : deltaA : yup;
lx = length(x_line);
ly = length(y_line);
x = ones(ly,1)*x_line;
y = y_line'*ones(1,lx);
%% BP����
R = zeros(ly,lx);                                                          % �洢б������                                                
BP_Image = zeros(ly,lx);                                                   % �洢BP������
dtr = c/(2*fs*interp);                                                     % �����ľ�����������
for index = 1:Na
    yradar = position_radar+vs*ta_rd(index);
    R = sqrt((x-distance_radar).^2 + (y-yradar).^2);                         % ���㵽�����ص��б������
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
xlabel('������루m��')
ylabel('������루m��')
title('BP������')
%% ��ά��ײ�ֵ������
% ���Ƚ��в��㣬��֤Ƶ����ԭ����λ�ã����б��ദ��
BP_Image_add = zeros(ly,lx);                                         
for ii = 1:ly
    for jj = 1:lx
        BP_Image_add(ii,jj) = BP_Image(ii,jj).*exp(-1i*4*pi*R(ii,jj)/lambda);
    end
end
%������� 
interp_new = 8;    
% �Ծ��������Ƶ����Ĳ���
lx_new = interp_new*lx;
lx_res = lx_new-lx;
fBP_Image_add = fft(BP_Image_add,[],2);
fBP_Image_add1 = [ fBP_Image_add(:,1:round(lx/2)) , zeros(ly,lx_res) , fBP_Image_add(:,round(lx/2)+1:lx)]; %Ƶ����
BP_Image_add1 = ifft(fBP_Image_add1,[],2);  %Ƶ����������
% �Է�λ�����Ƶ����Ĳ���
ly_new = interp_new*ly;
ly_res = ly_new-ly;
fBP_Image_add1 = fft(BP_Image_add1);
fBP_Image_add2 = [ fBP_Image_add1(1:round(ly/2),:) ; zeros(ly_res,lx_new) ; fBP_Image_add1(round(ly/2)+1:ly,:)]; %Ƶ����
BP_Image_final = ifft(fBP_Image_add2);  %��λ��Ƶ����������
figure;
imagesc(abs(BP_Image_final)) 
xlabel('�������')
ylabel('�������')
title('BP����������������')
figure;
mesh(abs(BP_Image_final));
xlabel('�������')
ylabel('�������')
title('BP�㷨����չ����')
toc