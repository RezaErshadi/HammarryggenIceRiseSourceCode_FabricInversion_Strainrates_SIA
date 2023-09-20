close all
clear
clc
P = ThisIsAn_ApRES_Script(mfilename('fullpath'));
%% Dir
load('SeismicColorMap100.mat');
% Load the inversion results for all the pRES points along the profile
DataDir = fullfile(P.Project,'results','InvRes','HIR_Profile_Estimated.mat');
load(DataDir);

HIFA = HIR_Prof.HA;
VIFA = HIR_Prof.EigVal3 - HIR_Prof.EigVal2;
RP = HIR_Prof.RadarPos;
Z = HIR_Prof.Zinv;
s1 = HIR_Prof.EigVal1;
s2 = HIR_Prof.EigVal2;
s3 = HIR_Prof.EigVal3;
lns3s2 = log(s3./s2);
lns2s1 = log(s2./s1);
lns3s1 = log(s3./s1);
K = lns3s2 ./ lns2s1;
C = lns3s1;
v1 = HIR_Prof.v1;
v2 = v1+90;
avgV1 = mean(v1);
stdV1 = std(v1);
avgV2 = avgV1+90;
%%
[xSEBETK,pRESpoints,H] = func_getTopoParam(P,RP./1000);
SE = pRESpoints(3,:);
BE = pRESpoints(4,:);

xlbl(1,:) = [-4.5 -3 -2 -1 0 1 2 3 4.5];
xlbl(2,:) = round([-4.5 -3 -2 -1 0 1 2 3 4.5].*H,0)./1000;
xlbl(3,:) = xlbl(2,:) + pRESpoints(2,8);
%% SIA
% Load the SIA results
SIAdir = fullfile(P.Project,'data','gis','SIA','Profile_csv','pRESline_interpolated-1m_SIA_Clara_EditVer.csv');
SIA = readtable(SIAdir);
%% interpolate pRES line with 1m resolution
% llX = HIR_Prof.LatLong(3,:);
% llY = HIR_Prof.LatLong(4,:);
% vx = min(llX):1:max(llX);
% vy = interp1(llX,llY,vx,'linear');
% figure,
% plot(vx,vy,'.k')
% hold on
% plot(llX,llY,'or')
%%
% Interpolate the inverted eigenvalues and fabric orientation along the profile
X = repmat(RP,size(Z,1),1);
Y = repmat(Z,1,size(RP,2));
[Xq,Yq] = meshgrid([RP(1):10:RP(end)],[Z(1):1:Z(end)]);
ip_method = 'linear';
ip_s1 = interp2(X,Y,s1,Xq,Yq,ip_method);
ip_s2 = interp2(X,Y,s2,Xq,Yq,ip_method);
ip_s3 = interp2(X,Y,s3,Xq,Yq,ip_method);
ip_HIFA = ip_s2-ip_s1;
ip_v2 = interp2(X,Y,v2,Xq,Yq,ip_method);
%% 
% interpolating the surface and bed elevation
ip_SE = interp1(RP,SE,Xq(1,:));
ip_BE = interp1(RP,BE,Xq(1,:));
Elv = ip_SE - Yq;
% H = mean(ip_SE(:) - ip_BE(:));
NrmZ = (H-Yq(:,1))./H;
xlbl = round([-4.5 -3 -2 -1 0 1 2 3 4.5].*H,0);
% ylblHIFA = round((H-[100 300 350 400])./H,2);
% ylblVSR = round((H-[100 150 200 250 300])./H,2);
%%
VSRDir = fullfile(P.Project,'results','StrainRate','Falk_July2023','correctedStrainRates_late.mat');
VSR = load(VSRDir);
VSR = VSR.data;

vsr = VSR.vsr_fit;
vsre = VSR.vsre_fit;
ylblVSR = [0.8 0.7 0.6 0.5 0.4 0.3];

% [XXX,YYY] = meshgrid(RP,VSR.fine_Ygrid(:,1));
% vsr = VSR.fine_vsr(:,1:15);

% figure,
% pc = pcolor(XXX,YYY,vsr);
% set(pc, 'EdgeColor', 'none');
% set(gca,'YDir','reverse')

% mmVSR = VSR.fine_vsr(:,1:15);
% zVSR = VSR.fine_Ygrid(:,1);

% ylblVSR = [0.9 0.7 0.5 0.45];
% for i = 1:length(ylblVSR)
%     temp1(:,i) = [ylblVSR(i)+0.025 ; ylblVSR(i)-0.025];
% end
% VSRzAVG = H - (temp1 .* H);
% 
% for i = 1:15
%     for ii = 1:size(VSRzAVG,2)
%         [~,i1] = min(abs(zVSR-VSRzAVG(1,ii)));
%         [~,i2] = min(abs(zVSR-VSRzAVG(2,ii)));
%         mVSR(ii,i) = mean(mmVSR(i1:i2,i),'omitnan');
%     end
% end

% ylblVSR = [0.8 0.7 0.6 0.5 0.4 0.3];
% VSRzAVG = H - (ylblVSR .* H);
% mVSR = nan(length(VSRzAVG),15);
% for i = 1:15
%     i1 = 1;
%     for ii = 1:length(VSRzAVG)
%         [~,i2] = min(abs(zVSR-VSRzAVG(ii)));
%         mVSR(ii,i) = mean(mmVSR(i1:i2,i),'omitnan');
%         i1 = i2+1;
%     end
% end

% VSR = load(VSRDir);
% mmVSR = VSR.VSRmm;
% zVSR = VSR.Ymm_vsr;
% dVSR = zVSR(1,:)-zVSR;
% topGap = round(SE - VSR.Ymm_vsr(1,:));
% ylblVSR = [0.8 0.6 0.5 0.4 0.35 0.3];
% VSRzAVG = H - (ylblVSR .* H);
% ylblVSR = round((H-[100 150 200 250 300])./H,2);
% a = [100 150 200 250 300];
% mVSR = nan(length(VSRzAVG),15);
% for i = 1:15
%     i1 = 1;
%     for ii = 1:length(VSRzAVG)
%         [~,i2] = min(abs(dVSR(:,i)-VSRzAVG(ii)));
%         mVSR(ii,i) = mean(mmVSR(i1:i2,i),'omitnan');
%         i1 = i2+1;
%     end
% end
%%
% Extract horizontal profile at specific depths
ylblHIFA = [0.8 0.7 0.6 0.5 0.4 0.3];
HIFAzAVG = H - (ylblHIFA .* H);
% a = [100 300 350 400];
% mHIFA = nan(length(HIFAzAVG),15);
% for i = 1:15
%     i1 = 1;
%     for ii = 1:length(HIFAzAVG)
%         [~,i2] = min(abs(Z-HIFAzAVG(ii)));
%         mHIFA(ii,i) = mean(HIFA(i1:i2,i),'omitnan');
%         i1 = i2+1;
%     end
% end
% use interpolated data
mHIFA = nan(length(HIFAzAVG),size(ip_HIFA,2));
i1 = 1;
for ii = 1:length(HIFAzAVG)
    [~,i2] = min(abs(Yq(:,1)-HIFAzAVG(ii)));
    mHIFA(ii,:) = mean(ip_HIFA(i1:i2,:),'omitnan');
    i1 = i2+1;
end
%%
% Extract l3 profile at specific depths
ylbll3 = ylblHIFA;
l3zAVG = H - (ylbll3 .* H);
ml3 = nan(length(HIFAzAVG),size(ip_s3,2));
i1 = 1;
for ii = 1:length(HIFAzAVG)
    [~,i2] = min(abs(Yq(:,1)-HIFAzAVG(ii)));
    ml3(ii,:) = mean(ip_s3(i1:i2,:),'omitnan');
    i1 = i2+1;
end
%%
xc = HIR_Prof.LatLong(3,8);
yc = HIR_Prof.LatLong(4,8);
dXY = round(sqrt(((xc-SIA.X).^2)+((yc-SIA.Y).^2)),0);
[~,icent] = min(dXY);
d1 = Xq(1,:) - Xq(1,1);
d2 = round(sqrt(((SIA.X(1)-SIA.X).^2)+((SIA.Y(1)-SIA.Y).^2)),0);
for i = 1:length(d1)
    [~,a] = min(abs(d1(i) - d2));
    dd = dXY(a);
    if a < icent
        dd = - dXY(a);
    end
    flowinfo(i,:) = [SIA.X(a), SIA.Y(a), dd, SIA.FD_CCW_TN(a), SIA.SD_CCW_TN(a)];
    flowdirectionatpoint = flowinfo(i,end-1);
    straindirectionatpoint = flowinfo(i,end);
    SIA_(1,i) = dd;
    SIA_(2,i) = flowdirectionatpoint;
    SIA_(3,i) = straindirectionatpoint;
    V2vsFD(:,i) = abs(flowdirectionatpoint - ip_v2(:,i));
    V2vsSD(:,i) = abs(straindirectionatpoint - ip_v2(:,i));
end

ProfileLength = Xq(1,:);
v2_Average = mean(ip_v2);
FlowDirection_Average = SIA_(2,:);
StrainDirection_Average = SIA_(3,:);

v2_VS_FD = abs(v2_Average-FlowDirection_Average);
v2_VS_FD = abs(90 - abs(90 - v2_VS_FD));
% 
v2_VS_SD = abs(v2_Average-StrainDirection_Average);
v2_VS_SD = abs(90 - abs(90 - v2_VS_SD));
% 
figure,
subplot(2,1,1)
plot(ProfileLength,v2_Average);
hold on
plot(ProfileLength,FlowDirection_Average);
plot(ProfileLength,StrainDirection_Average);
ylim([0 360])

subplot(2,1,2)
plot(ProfileLength,v2_VS_FD);
hold on 
plot(ProfileLength,v2_VS_SD);
plot([Xq(1) Xq(end)] , [90 90])


% subplot(3,1,3)
% plot(ProfileLength,v2_VS_SD);
% hold on 
% plot([Xq(1) Xq(end)] , [90 90])
% ylim([0 360])


% V2vsFD(V2vsFD>90) = V2vsFD(V2vsFD>90) - 180;
% V2vsFD = abs(V2vsFD);
% 
% V2vsSD(V2vsSD>90) = V2vsSD(V2vsSD>90) - 180;
% V2vsSD = abs(V2vsSD);

mean(v2_VS_FD)
mean(v2_VS_FD(1:249))
mean(v2_VS_FD(251:end))

mean(v2_VS_SD)
mean(v2_VS_SD(1:249))
mean(v2_VS_SD(251:end))
%%
fig1 = CLASS_FixedPlot.SetFigureSize(0,0,0.4,0.95);
fnt = 18;
rw = 8;
cl = 1;
tmp = [ 1 ylblHIFA ];
for i = 1:length(tmp)-1
    lgndtxt1(i) = string(tmp(i))+"-"+string(tmp(i+1));
end
ii = 1;
ax{ii} = subplot(rw,cl,[1:2]);
LineStyle = [":";":";":";":";"-";"-"];
mycolors = [0 0.7 0; 1 0.5 0; 0 0 0.8;  0.8 0 0; 0.6 0 0.8;  0 0 0];
for i = 1:size(mHIFA,1)
    plot(ax{ii},Xq(1,:),mHIFA(i,:),'LineWidth',2,'LineStyle',LineStyle(i),'Color',mycolors(i,:))
    hold(ax{ii},'on')
end
lgnd1 = legend(ax{ii},lgndtxt1,'Location','NorthWest','NumColumns',2);
title(lgnd1,'z')
ylim(ax{ii},[0 0.1])
yticks(ax{ii},[0:0.025:0.1])
yticklabels(ax{ii},{'0','0.025','0.05','0.075','0.1'})
ylabel(ax{ii},'$\Delta\lambda_H$ [-]','Interpreter','latex')
xticks(ax{ii},xlbl)
% xticklabels(ax{ii},string(xlbl)+"m")
xlim(ax{ii},[-2600 2600])
xticklabels(ax{ii},[])
set(ax{ii},'FontSize',fnt)
title(ax{ii},'(a)')
%-----------------------------------
tmp = [ 1 ylblHIFA ];
for i = 1:length(tmp)-1
    lgndtxt2(i) = string(tmp(i))+"-"+string(tmp(i+1));
end
ii = 2;
ax{ii} = subplot(rw,cl,[3:4]);
for i = 1:size(mHIFA,1)
    plot(ax{ii},Xq(1,:),ml3(i,:),'LineWidth',2,'LineStyle',LineStyle(i),'Color',mycolors(i,:))
    hold(ax{ii},'on')
end
% plot(ax{ii},Xq(1,:),ml3,'-','LineWidth',2)
% legend(ax{ii},{'0-100     [m]','100-300 [m]','300-350 [m]','350-400 [m]'},'Location','NorthWest')
lgnd2 = legend(ax{ii},lgndtxt2,'Location','Northeast','NumColumns',2);
title(lgnd2,'z')
ylim(ax{ii},[0.33 0.75])
yticks(ax{ii},[0.33,0.5,0.75])
% yticklabels(ax{ii},{'0.25','0.5','0.75','1'})
% ylabel(ax{ii},'\lambda3 [-]')
ylabel(ax{ii},'$\lambda_3$ [-]','Interpreter','latex')
xticks(ax{ii},xlbl)
xlim(ax{ii},[-2600 2600])
% xticklabels(ax{ii},string(xlbl)+"m")
xticklabels(ax{ii},[])
set(ax{ii},'FontSize',fnt)
title(ax{ii},'(b)')
%-----------------------------------
SIA_(2,(SIA_(2,:)<=90)) = SIA_(2,(SIA_(2,:)<=90))+180;
% SIA_(3,(SIA_(3,:)>=270)) = SIA_(3,(SIA_(3,:)>=270))-90;
ii = 3;
ax{ii} = subplot(rw,cl,[5:6]);
plot(ax{ii},[RP(1) RP(end)],[180 180],'--k','LineWidth',2)
hold(ax{ii},'on')
plot(ax{ii},ProfileLength,v2_Average,'-b','LineWidth',2);
plot(ax{ii},ProfileLength,FlowDirection_Average,':r','LineWidth',2);
plot(ax{ii},ProfileLength,StrainDirection_Average,'-r','LineWidth',2);
ylim(ax{ii},[0 360])
yticks(ax{ii},[0 45 90 135 180 225 270 315 360])
yticklabels(ax{ii},{'N-S','-135','-90','-45','N-S','45','90','135','N-S'})
ylabel(ax{ii},'Orientation [deg]','Interpreter','latex')
xticks(ax{ii},xlbl)
xlim(ax{ii},[-2600 2600])
xticklabels(ax{ii},[])
legend(ax{ii},{'N-S','Mean $\vec{v}_2$ direction (pRES)','Surface flow direction (SIA)','Max horizontal strain rate direction (SIA)'},'Location','Northeast','NumColumns',2,'Interpreter','latex')
set(ax{ii},'FontSize',fnt)
title(ax{ii},'(c)')

% yyaxis(ax{ii},'left')
% plot(ax{ii},[RP(1) RP(end)],[180 180],'--k','LineWidth',2)
% hold(ax{ii},'on')
% plot(ax{ii},ProfileLength,v2_Average,'-','LineWidth',2);
% ylim(ax{ii},[45 315])
% yticks(ax{ii},[45 90 135 180 225 270 315])
% yticklabels(ax{ii},{'-135','-90','-45','N-S','45','90','135'})
% ylabel(ax{ii},'$\vec{v}_2$ [deg]','Interpreter','latex')
% yyaxis(ax{ii},'right')
% plot(ax{ii},ProfileLength,v2_VS_FD,':','LineWidth',2);
% plot(ax{ii},ProfileLength,v2_VS_SD,'-','LineWidth',2);
% ylim(ax{ii},[0 90])
% ylabel(ax{ii},'Deviation from $\vec{v}_2$ [deg]','Interpreter','latex')
% yticks(ax{ii},[0:15:90])
% % yticklabels(ax{ii},{'-135','15','30','45','60','75','90'})
% % xticks(ax{ii},RP([1,8,15]))
% % xticklabels(ax{ii},HIR_Prof.SiteName([1,8,15]))
% xticks(ax{ii},xlbl)
% xlim(ax{ii},[-2600 2600])
% xticklabels(ax{ii},[])
% % legend(ax{ii},{'N-S','mean v2 (pRES)','Flow direction (SIA)','max strain direction (SIA)'},'Location','SouthWest','NumColumns',2)
% set(ax{ii},'FontSize',fnt)
% title(ax{ii},'(c)')
%-----------------------------------
tmp = [ 1 ylblVSR ];
for i = 1:length(tmp)-1
    lgndtxt4(i) = string(tmp(i))+"-"+string(tmp(i+1));
end
ii = 4;
ax{ii} = subplot(rw,cl,[7:8]);
% plot(RP,mVSR,'.-','LineWidth',2,'MarkerSize',30)


for i = 1:length(ylblVSR)
    errorbar(RP,vsr(i,:),vsre(i,:),'.-','LineWidth',2,'MarkerSize',20,'LineStyle',LineStyle(i),'Color',mycolors(i,:))
hold on
end

lgnd4 = legend(ax{ii},lgndtxt4,'Location','SouthWest','NumColumns',2);
title(lgnd4,'z')
% legend(ax{ii},{'0-100     [m]','100-150 [m]','150-200 [m]','200-250 [m]','250-300 [m]'},'Location','SouthWest')


minVSR = min(vsr(:));
% maxVSR = max(mVSR(:));
ylim(ax{ii},[(minVSR+(0.25*minVSR)) 0])
ylabel(ax{ii},'$\mathbf{\dot{\varepsilon}}_{zz}$ [$a^{-1}$]','Interpreter','latex')
xticks(ax{ii},xlbl)
xlim(ax{ii},[-2600 2600])
% cellArray = num2cell([string(abs(round(xlbl/H, 1))) + newline + "(" + string(xlbl) + "m)"]);
% cellArray = num2cell([string(abs(round(xlbl/H,1)))+"\n("+string(xlbl)]+"m)");
% charArray = cell(size(cellArray));
% for i = 1:numel(cellArray)
%     charArray{i} = char(cellArray{i});
% end
% xticklabels(ax{ii},charArray)
xticklabels(ax{ii},string(abs(round(xlbl/H, 1))))
xtickangle(ax{ii},90)
xlabel(ax{ii},'x [-]')
set(ax{ii},'FontSize',fnt)
xtickangle(ax{ii},0)
title(ax{ii},'(d)')
%%
% v2_vs_NS = mean(ip_v2)-180;
% mean(v2_vs_NS)
% mean(v2_vs_NS(1:201))
% mean(v2_vs_NS(302:end))

% v2_vs_FD = mean(ip_v2) - SIA_(2,:);
% mean(v2_vs_FD)
% mean(v2_vs_FD(1:251))
% mean(v2_vs_FD(251:end))

v2_vs_Strain = mean(ip_v2) - SIA_(3,:);
mean(v2_vs_Strain)
mean(v2_vs_Strain(1:251))
mean(v2_vs_Strain(251:end))

%%
% fig0.InvertHardcopy = 'off';
% print(fig0,"HIR_Profile_Int2D.png",'-dpng','-r300')

fig1.InvertHardcopy = 'off';
% print(fig1,"HIR_Profile_DepthAverage_2.png",'-dpng','-r300')