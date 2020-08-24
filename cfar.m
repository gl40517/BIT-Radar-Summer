function [OUT, x] = cfar( Img )
%% CFAR函数主要用于目标检测
%   Pfa为虚警概率
%   Img为待检测图像
Pfa = 1e-6;
[m, n] = size(Img);
Img = Img/max(max(Img));
%% 设定自适应窗口大小    m为长，n为宽
Len_BackWin = floor(n/16);
Len_AimWin = floor(0.15*Len_BackWin); %目标窗Length现阶段必须为奇数
if mod(Len_AimWin,2) == 0
    Len_AimWin = Len_AimWin+1;
end
Width_Win = floor(m/15);
Img = padarray(Img,[0,Len_BackWin],'both');   %补0
[m, n] = size(Img);
Len_Img2 = m;
Len_Img1 = n;
%% 设置目标窗和背景窗的大小（预设矩阵，窗大小）
G_Aim = zeros(Width_Win, Len_AimWin);
G_Back1 =  zeros(Width_Win, Len_BackWin);
G_Back2 = G_Back1;
G_Back = [G_Back1,G_Back2];
%设置输出
OUT1 = Img;
%设置循环步长
step = Width_Win;

%% ca-cfar
for i = Len_BackWin+floor(Len_AimWin/2)+1 : Len_Img1-(Len_BackWin+floor(Len_AimWin/2))%从左到右（列循环）
    for j = 1:step:Len_Img2-Width_Win+1 %从上到下（行循环）
        G_Aim = Img(j:j+Width_Win-1, i-floor(Len_AimWin/2):i+floor(Len_AimWin/2)); %目标窗(长、宽)
        G_Back1 = Img(j:j+Width_Win-1, i-floor(Len_AimWin/2)-Len_BackWin:i-floor(Len_AimWin/2)-1); %前背景窗(宽度和目标窗一样)
        G_Back2 = Img(j:j+Width_Win-1, i+floor(Len_AimWin/2)+1:i+floor(Len_AimWin/2)+Len_BackWin); %后背景窗(宽度和目标窗一样)
        G_Back = [G_Back1,G_Back2];
        Mean_Back = mean(mean(G_Back)); %背景窗均值
        [m1,n1] = size(G_Back);
        [m2,n2] = size(G_Aim);
        BackSize = m1*n1;
        AimSize = m2*n2;
        R = BackSize+AimSize;              %总点数
        TN = BackSize*(Pfa^(-1/R)-1).*Mean_Back;       %判决门限
        %判断是否为目标
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

%% 点迹凝聚
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
