function [OUT, x] = cfar( Img )
%% CFAR������Ҫ����Ŀ����
%   PfaΪ�龯����
%   ImgΪ�����ͼ��
Pfa = 1e-6;
[m, n] = size(Img);
Img = Img/max(max(Img));
%% �趨����Ӧ���ڴ�С    mΪ����nΪ��
Len_BackWin = floor(n/16);
Len_AimWin = floor(0.15*Len_BackWin); %Ŀ�괰Length�ֽ׶α���Ϊ����
if mod(Len_AimWin,2) == 0
    Len_AimWin = Len_AimWin+1;
end
Width_Win = floor(m/15);
Img = padarray(Img,[0,Len_BackWin],'both');   %��0
[m, n] = size(Img);
Len_Img2 = m;
Len_Img1 = n;
%% ����Ŀ�괰�ͱ������Ĵ�С��Ԥ����󣬴���С��
G_Aim = zeros(Width_Win, Len_AimWin);
G_Back1 =  zeros(Width_Win, Len_BackWin);
G_Back2 = G_Back1;
G_Back = [G_Back1,G_Back2];
%�������
OUT1 = Img;
%����ѭ������
step = Width_Win;

%% ca-cfar
for i = Len_BackWin+floor(Len_AimWin/2)+1 : Len_Img1-(Len_BackWin+floor(Len_AimWin/2))%�����ң���ѭ����
    for j = 1:step:Len_Img2-Width_Win+1 %���ϵ��£���ѭ����
        G_Aim = Img(j:j+Width_Win-1, i-floor(Len_AimWin/2):i+floor(Len_AimWin/2)); %Ŀ�괰(������)
        G_Back1 = Img(j:j+Width_Win-1, i-floor(Len_AimWin/2)-Len_BackWin:i-floor(Len_AimWin/2)-1); %ǰ������(��Ⱥ�Ŀ�괰һ��)
        G_Back2 = Img(j:j+Width_Win-1, i+floor(Len_AimWin/2)+1:i+floor(Len_AimWin/2)+Len_BackWin); %�󱳾���(��Ⱥ�Ŀ�괰һ��)
        G_Back = [G_Back1,G_Back2];
        Mean_Back = mean(mean(G_Back)); %��������ֵ
        [m1,n1] = size(G_Back);
        [m2,n2] = size(G_Aim);
        BackSize = m1*n1;
        AimSize = m2*n2;
        R = BackSize+AimSize;              %�ܵ���
        TN = BackSize*(Pfa^(-1/R)-1).*Mean_Back;       %�о�����
        %�ж��Ƿ�ΪĿ��
        for a = 0:Len_AimWin-1
            for b = 0:Width_Win-1
                if G_Aim(b+1, a+1) < TN
                    OUT1(j+b, i-floor(Len_AimWin/2)+a) = 0; 
                else
                     OUT1(j+b, i-floor(Len_AimWin/2)+a) = G_Aim(b+1, a+1);
                end
            end
        end
    end
end

%% �㼣����
[x0,y0] = find(  abs(OUT1) == max(max(abs(OUT1))) );
OUT1(x0,y0) = 1;
for i = 1:Len_Img2
    for j = 1:Len_Img1
        if OUT1(i,j) ~= 1
            OUT1(i,j) = 0;
        end
    end
end
OUT1(:,Len_Img1-Len_BackWin+1:Len_Img1) = [];
OUT1(:,1:Len_BackWin) = [];
x = y0;
OUT = OUT1;
