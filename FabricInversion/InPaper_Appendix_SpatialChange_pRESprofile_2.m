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
BEfromRADAR = HIR_Prof.SBR_Elv(3,:); % Bed elevation at each site from pRES data

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
% Extract horizontal profile at specific depths
ylblHIFA = [0.8 0.6 0.5 0.4 0.35 0.3];
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

v2_VS_FD = abs(90 - abs(90 - V2vsFD));
v2_VS_SD = abs(90 - abs(90 - V2vsSD));

% ProfileLength = Xq(1,:);
% v2_Average = mean(ip_v2);
% FlowDirection_Average = SIA_(2,:);
% StrainDirection_Average = SIA_(3,:);
% 
% v2_VS_FD = abs(v2_Average-FlowDirection_Average);
% v2_VS_FD = abs(90 - abs(90 - v2_VS_FD));
% % 
% v2_VS_SD = abs(v2_Average-StrainDirection_Average);
% v2_VS_SD = abs(90 - abs(90 - v2_VS_SD));
% 
% 
% 
% V2vsFD(V2vsFD>90) = V2vsFD(V2vsFD>90) - 180;
% V2vsFD = abs(V2vsFD);
% 
% V2vsSD(V2vsSD>90) = V2vsSD(V2vsSD>90) - 180;
% V2vsSD = abs(V2vsSD);
%%
fig0 = CLASS_FixedPlot.SetFigureSize(0.05,0.05,0.55,0.8);
rw = 4;
cl = 4;
% CM1 = 'gnuplot2_r';
CM2 = 'cividis_r';
% --------------------------------------
ii = 1;
ax{ii} = subplot(rw,cl,[1,2,5,6]);
hold(ax{ii},'on')
plot(ax{ii},Xq(1,:),ip_SE,'-k')%,'linewidth',2)
plot(ax{ii},Xq(1,:),ip_BE,':k')%,'linewidth',2)

plot(ax{ii},[RP(1) RP(end)],[mean(ip_SE(:)) mean(ip_SE(:))],'-b','linewidth',2)
plot(ax{ii},[RP(1) RP(end)],[mean(ip_BE(:)) mean(ip_BE(:))],':b','linewidth',2)

plot(RP,SE,'or','MarkerFaceColor','r','MarkerSize',5)
plot(RP,BEfromRADAR,'xr','LineWidth',2)

pc = pcolor(ax{ii},Xq,Elv,ip_HIFA);
set(pc, 'EdgeColor', 'none');
ncm = getPyPlot_cMap('viridis',20);
colormap(ax{ii},flip(ncm))
cb = colorbar(ax{ii});
cb.Label.String = '\Delta\lambda_{H} [-]';
% cb.Location = "southoutside";
cb.Ruler.TickLabelRotation=90;
cb.Ticks = [0 :0.1:0.5];
% cb.TickLabels = {'0.33','0.50','0.75','1.00'};
caxis(ax{ii},[0 0.5])
xticks(ax{ii},xlbl)
% xticklabels(ax{ii},[string(abs(round(RP(1)/H,1))),'',string(abs(round(RP(3)/H,1))),'',string(abs(round(RP(5)/H,1))),'','',string(abs(round(RP(8)/H,1))),'','',string(abs(round(RP(11)/H,1))),'',string(abs(round(RP(13)/H,1))),'',string(abs(round(RP(end)/H,1)))])
% xtickangle(ax{ii},90)
% xlabel(ax{ii},'x')
xticklabels(ax{ii},[])
xlim(ax{ii},[min(RP)-50 max(RP)+50])
ylim(ax{ii},[min(BE)-25 max(SE)+25])
ytickangle(ax{ii},90)
newYticks = linspace(mean(ip_BE(:)),mean(ip_SE(:)),5);
newDepth = linspace(0,H,5)./H;
yticks(ax{ii},newYticks)
yticklabels(ax{ii},newDepth)
ylabel(ax{ii},'z [-]')
title(ax{ii},'(a)')
legend(gca,'Surface topo.','Bed topo.','Mean surface topo. (z=1)','Mean bed topo. (z=0)','pRES points','Detected bed','Location','southwest','NumColumns',3,'Position',[0.182346723044398 0.605778196322641 0.223572938689218 0.0452003727865795]);
set(ax{ii},'FontSize',20)
% --------------------------------------
ii = 2;
ax{ii} = subplot(rw,cl,[3:4,7,8]);
hold(ax{ii},'on')
pc = pcolor(ax{ii},Xq,Elv,ip_s3);
set(pc, 'EdgeColor', 'none');
plot(ax{ii},Xq(1,:),ip_SE,'-k')%,'linewidth',2)
plot(ax{ii},Xq(1,:),ip_BE,':k')%,'linewidth',2)
plot(RP,SE,'or','MarkerFaceColor','r','MarkerSize',5)
plot(RP,BEfromRADAR,'xr','LineWidth',2)
plot(ax{ii},[RP(1) RP(end)],[mean(ip_SE(:)) mean(ip_SE(:))],'-b','linewidth',2)
plot(ax{ii},[RP(1) RP(end)],[mean(ip_BE(:)) mean(ip_BE(:))],':b','linewidth',2)
ncm = getPyPlot_cMap('viridis',14);
colormap(ax{ii},flip(ncm))
cb = colorbar(ax{ii});
cb.Label.String = '\lambda3 [-]';
% cb.Location = "southoutside";
cb.Ticks = [0.3 0.5 0.75 1];
cb.TickLabels = {'0.33','0.50','0.75','1.00'};
cb.Ruler.TickLabelRotation=90;
caxis(ax{ii},[0.3 1])
xticks(ax{ii},xlbl)
% xticklabels(ax{ii},[string(abs(round(RP(1)/H,1))),'',string(abs(round(RP(3)/H,1))),'',string(abs(round(RP(5)/H,1))),'','',string(abs(round(RP(8)/H,1))),'','',string(abs(round(RP(11)/H,1))),'',string(abs(round(RP(13)/H,1))),'',string(abs(round(RP(end)/H,1)))])
% xtickangle(ax{ii},90)
xticklabels(ax{ii},[])
xlim(ax{ii},[min(RP)-50 max(RP)+50])
% xlabel(ax{ii},'X/H')
ylim(ax{ii},[min(BE)-25 max(SE)+25])
yticks(ax{ii},newYticks)
yticklabels(ax{ii},[])
title(ax{ii},'(b)')
set(ax{ii},'FontSize',20)
% % --------------------------------------
ii = 3;
ax{ii} = subplot(rw,cl,[9,10,13,14]);
hold(ax{ii},'on')
pc = pcolor(ax{ii},Xq,Elv,v2_VS_FD);
set(pc, 'EdgeColor', 'none');
plot(ax{ii},Xq(1,:),ip_SE,'-k')%,'linewidth',2)
plot(ax{ii},Xq(1,:),ip_BE,':k')%,'linewidth',2)
plot(RP,SE,'or','MarkerFaceColor','r','MarkerSize',5)
plot(RP,BEfromRADAR,'xr','LineWidth',2)
plot(ax{ii},[RP(1) RP(end)],[mean(ip_SE(:)) mean(ip_SE(:))],'-b','linewidth',2)
plot(ax{ii},[RP(1) RP(end)],[mean(ip_BE(:)) mean(ip_BE(:))],':b','linewidth',2)
ylim(ax{ii},[min(BE)-25 max(SE)+25])
yticks(ax{ii},newYticks)
yticklabels(ax{ii},[])
xticks(ax{ii},xlbl)
xticklabels(ax{ii},[string(abs(round(xlbl/H,1)))])
% xticklabels(ax{ii},[string(abs(round(RP(1)/H,1))),'',string(abs(round(RP(3)/H,1))),'',string(abs(round(RP(5)/H,1))),'','',string(abs(round(RP(8)/H,1))),'','',string(abs(round(RP(11)/H,1))),'',string(abs(round(RP(13)/H,1))),'',string(abs(round(RP(end)/H,1)))])
ytickangle(ax{ii},90)
xlabel(ax{ii},'x [-]')
xlim(ax{ii},[min(RP)-50 max(RP)+50])
yticks(ax{ii},newYticks)
yticklabels(ax{ii},newDepth)
% ytickangle(ax{ii},90)
ylabel(ax{ii},'z [-]')
set(ax{ii},'FontSize',20)
ncm = getPyPlot_cMap(CM2,9);
colormap(ax{ii},ncm)
caxis(ax{ii},[0 90])
cb = colorbar(ax{ii});
% cb.Label.String = 'v2 relative to flow direction [deg]';
cb.Label.String = '[deg]';
cb.Ticks = [0:30:90];
% cb.Location = "southoutside";
title(ax{ii},'(c)')
set(ax{ii},'FontSize',20)
% % --------------------------------------
ii = 4;
ax{ii} = subplot(rw,cl,[11,12,15,16]);
hold(ax{ii},'on')
pc = pcolor(ax{ii},Xq,Elv,v2_VS_SD);
set(pc, 'EdgeColor', 'none');
plot(ax{ii},Xq(1,:),ip_SE,'-k')%,'linewidth',2)
plot(ax{ii},Xq(1,:),ip_BE,':k')%,'linewidth',2)
plot(RP,SE,'or','MarkerFaceColor','r','MarkerSize',5)
plot(RP,BEfromRADAR,'xr','LineWidth',2)
plot(ax{ii},[RP(1) RP(end)],[mean(ip_SE(:)) mean(ip_SE(:))],'-b','linewidth',2)
plot(ax{ii},[RP(1) RP(end)],[mean(ip_BE(:)) mean(ip_BE(:))],':b','linewidth',2)
ylim(ax{ii},[min(BE)-25 max(SE)+25])
yticks(ax{ii},newYticks)
yticklabels(ax{ii},[])
xticks(ax{ii},RP)
xticklabels(ax{ii},[])
set(ax{ii},'FontSize',20)
xticks(ax{ii},xlbl)
xticklabels(ax{ii},[string(abs(round(xlbl/H,1)))])
% xtickangle(ax{ii},45)
xlabel(ax{ii},'x [-]')
xlim(ax{ii},[min(RP)-50 max(RP)+50])
ncm = getPyPlot_cMap(CM2,9);
colormap(ax{ii},ncm)
caxis(ax{ii},[0 90])
cb = colorbar(ax{ii});
% cb.Label.String = 'v2 relative to e_{max} direction [deg]';
cb.Label.String = '[deg]';
cb.Ticks = [0:30:90];
% cb.Location = "southoutside";
% cb.Ruler.TickLabelRotation=90;
title(ax{ii},'(d)')
set(ax{ii},'FontSize',20)

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
fig0.InvertHardcopy = 'off';
print(fig0,"HIR_Profile_Int2D.png",'-dpng','-r300')

% fig1.InvertHardcopy = 'off';
% print(fig1,"HIR_Profile_DepthAverage_2.png",'-dpng','-r300')