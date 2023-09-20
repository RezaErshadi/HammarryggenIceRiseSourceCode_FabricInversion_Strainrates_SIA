clear
close all
clc
P = ThisIsAn_ApRES_Script(mfilename('fullpath'));
%% Read Inverted Fabric from ApRES
load('SeismicColorMap100.mat');
InversionPath = fullfile(P.Project,'results','InvRes','HIR_Profile_Estimated.mat');
load(InversionPath);

HIFA = HIR_Prof.HA;
l1_R = HIR_Prof.EigVal1(:,8);
l2_R = HIR_Prof.EigVal2(:,8);
l3_R = HIR_Prof.EigVal3(:,8);
lns3s2_R = log(l3_R./l2_R);
lns2s1_R = log(l2_R./l1_R);
k_R = lns3s2_R ./ lns2s1_R;
C_R = log(l3_R./l1_R);

RP = HIR_Prof.RadarPos./1000;
BEfromRADAR = HIR_Prof.SBR_Elv(3,:); % Bed elevation at each site from pRES data

[xSEBETK,pRESpoints,H] = func_getTopoParam(P,RP);
xq = xSEBETK(:,1);
xp7 = pRESpoints(2,1);
xp14 = pRESpoints(2,end);
lpRES = abs(xp14 - xp7);

p0Thickness = H;

Z_R = HIR_Prof.Zinv;
NrmZ = (p0Thickness-Z_R)./p0Thickness;
%% Ice Core analysis
% IceCorePath = fullfile(P.Project,'data','icecore','TIR_fabricGS_VTsibulskaya_Sheet1_EdittedByRE.csv');
IceCorePath = fullfile(P.Project,'data','icecore','TIR_fabricGS_VTsibulskaya_AVGSheet1_EdittedByRE.csv');
IceCore = readtable(IceCorePath);
l1_IC = table2array(IceCore(:,4));
l2_IC = table2array(IceCore(:,3));
l3_IC = table2array(IceCore(:,2));
lns3s2_IC = log(l3_IC./l2_IC);
lns2s1_IC = log(l2_IC./l1_IC);
k_IC = lns3s2_IC ./ lns2s1_IC;
C_IC = log(l3_IC./l1_IC);

% k_IC(k_IC<0.2) = 0; % UAG
% k_IC(k_IC>=0.2 & k_IC<0.5) = -1; % G
% k_IC(k_IC>=0.5 & k_IC<0.98) = -2; % G
% k_IC(k_IC>=0.98 & k_IC<=1.02) = -3; % GCTZ
% k_IC(k_IC>1.02 & k_IC<=2) = -4; % C
% k_IC(k_IC>2 & k_IC<=5) = -5; % C
% k_IC(k_IC>5) = -6; % UAC
% k_IC = -1*k_IC;

Z_IC = table2array(IceCore(:,1));
%% plot
f = CLASS_FixedPlot.SetFigureSize(0,0,0.60,0.8);
fsz = 20;
rw = 1;
cl = 4;

% ===============================================================
ax{1} = subplot(rw,cl,[1]);
i = 1;
plot(ax{i},l1_IC,Z_IC,'xk','linewidth',1.5,'MarkerSize',10)
hold(ax{i},'on')
plot(ax{i},l2_IC,Z_IC,'dk','linewidth',1.5,'MarkerSize',10)
plot(ax{i},l3_IC,Z_IC,'*k','linewidth',1.5,'MarkerSize',10)
plot(ax{i},l1_R,Z_R,'.-r','linewidth',1.5,'MarkerSize',25)
plot(ax{i},l2_R,Z_R,'.-b','linewidth',1.5,'MarkerSize',25)
plot(ax{i},l3_R,Z_R,'.-g','linewidth',1.5,'MarkerSize',25)
set(ax{i},'YDir','reverse')
xlim(ax{i},[0 1])
xticks(ax{i},[0:0.25:1])
ylim(ax{i},[0 401])
yticks(ax{i},[0:50:p0Thickness round(p0Thickness,0)])
ylabel(ax{i},'Depth [m]')
xlabel(ax{i},'\lambda [-]')
title(ax{i},'Eigenvalues')
legend(ax{i},'Meas. \lambda_1','Meas. \lambda_2','Meas. \lambda_3','Est. \lambda_1','Est. \lambda_2','Est. \lambda_3','Location','northeast')
set(ax{i},'FontSize',fsz)
box(ax{i},'off')
% ===============================================================
ax{2} = subplot(rw,cl,[2]);
i = 2;
plot(ax{i},l2_IC-l1_IC,Z_IC,'ok','linewidth',1.5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10)
hold(ax{i},'on')
plot(ax{i},l2_R-l1_R,Z_R,'.-','color','b','linewidth',1.5,'MarkerSize',25)

plot(ax{i},l3_IC-l2_IC,Z_IC,'>k','linewidth',1.5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10)
plot(ax{i},l3_R-l2_R,Z_R,'>-','color','r','MarkerFaceColor','r','linewidth',1.5,'MarkerSize',8)

set(ax{i},'YDir','reverse')
yticklabels(ax{i},[])
xlim(ax{i},[0 1])
xticks(ax{i},[0:0.25:1])
ylim(ax{i},[0 401])
xlabel(ax{i},'\Delta\lambda [-]')
title(ax{i},'Anisotropy')
legend(ax{i},'Meas. \Delta\lambda_H','Est. \Delta\lambda_H','Meas. \Delta\lambda_V','Est. \Delta\lambda_V','Location','northeast')
set(ax{i},'FontSize',fsz)
box(ax{i},'off')
% ===============================================================
ax{4} = subplot(rw,cl,[3,4]);
i = 4;
% lgn1 = semilogx(ax{i},[1e0 1e0],[0 p0Thickness],'--k','LineWidth',2);
lgn1 = plot(ax{i},[1 1],[0 p0Thickness],'--k','LineWidth',2);
hold(ax{i},'on')
nrmlzdC_IC = C_IC.*1500./7;
nrmlzdC_R = C_R.*1500./7;
% lgn3 = semilogx(ax{i},k_R,Z_R,'.-k','linewidth',1.5,'MarkerSize',15);
lgn3 = plot(ax{i},k_R,Z_R,'.-k','linewidth',1.5,'MarkerSize',15);
lgn4 = scatter(ax{i},k_IC,Z_IC,150,C_IC,'filled','>','MarkerEdgeColor',[0 0 0]);
lgn5 = scatter(ax{i},k_R,Z_R,150,C_R,'filled','MarkerEdgeColor',[0 0 0]);
colormap(ax{i},flip(gray))
caxis(ax{i},[0 3.5])
cb = colorbar(ax{i});
cb.Label.String = 'C [-]';
cb.Position = [0.862645348837207 0.692799315452397 0.0104166666666666 0.208412240931568];
cb.Ticks = [0 1 2 3 3.5];
cb.TickLabels = {'0','1','2','3','\uparrow'};
set(ax{i},'YDir','reverse')
yticklabels(ax{i},[])
xlim(ax{i},[0 35])
ylim(ax{i},[0 401])
xlabel(ax{i},'K [-]')
title(ax{i},'Fabric type')
legend(ax{i},[lgn1,lgn4,lgn5],{'k=1','Meas. K','Est. K'},'Position',[0.80281007751938 0.835041938490213 0.0523255813953488 0.0661696178937557])
% xticks(ax{i},[1e-2 2e-1 1e0 5e0 1e1 2e1 1e2])
% xticks(ax{i},[1e-2 1e-1 1e0 1e1 100])
% xticklabels(ax{i},{'0.01','0.1','1','10','100'})
% xticklabels(ax{i},{'0\Leftarrow','0.1','0.2','0.5','1','5','10','20','\Rightarrow\infty'})
xtickangle(ax{i},0)
set(ax{i},'FontSize',fsz)
box(ax{i},'off')

% annotation(f,'textbox',...
%     [0.800387596899219 0.137931034482757 0.0629844961240305 0.0326188257222736],...
%     'String',{'Girdle zone'},...
%     'Rotation',0,...
%     'FontSize',20,...
%     'FontName','Helvetica Neue',...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
% 
% annotation(f,'textbox',...
%     [0.883236434108525 0.130475302889093 0.0629844961240305 0.0326188257222736],...
%     'String','Cluster zone',...
%     'Rotation',0,...
%     'FontSize',20,...
%     'FontName','Helvetica Neue',...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
% 
% annotation(f,'textbox',...
%     [0.71317829457364 0.0512581547064293 0.0707364341085303 0.0326188257222738],...
%     'String','Uniaxial girdle',...
%     'FontSize',20,...
%     'FontName','Helvetica Neue',...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
% 
% annotation(f,'textbox',...
%     [0.867248062015501 0.0512581547064293 0.0760658914728712 0.0326188257222738],...
%     'String','Uniaxial cluster',...
%     'FontSize',20,...
%     'FontName','Helvetica Neue',...
%     'FitBoxToText','off',...
%     'EdgeColor','none');

% ===============================================================
f.InvertHardcopy = 'off';
print(f,"HIR_IceCoreVsRadar.png",'-dpng','-r300');

%%
% L1radarAtIcecore = interp1(Z_R,l1_R,Z_IC);
% L1 = [l1_IC,L1radarAtIcecore];
% % L1(:,3) = ((L1(:,1) - L1(:,2)) ./ L1(:,1)) .* 100;
% L1(:,3) = L1(:,1) - L1(:,2);
% meanL = mean(L1(:,3));
%  
% L2radarAtIcecore = interp1(Z_R,l2_R,Z_IC);
% L2 = [l2_IC,L2radarAtIcecore];
% % L2(:,3) = ((L2(:,1) - L2(:,2)) ./ L2(:,1)) .* 100;
% L2(:,3) = L2(:,1) - L2(:,2);
% meanL(2) = mean(L2(:,3));
% 
% L3radarAtIcecore = interp1(Z_R,l3_R,Z_IC);
% L3 = [l3_IC,L3radarAtIcecore];
% % L3(:,3) = ((L3(:,1) - L3(:,2)) ./ L3(:,1)) .* 100;
% L3(:,3) = L3(:,1) - L3(:,2);
% meanL(3) = mean(L3(:,3));

dlhradarAtIcecore = interp1(Z_R,l2_R-l1_R,Z_IC);
dlhAtIcandRad = [l2_IC-l1_IC,dlhradarAtIcecore];
% L3(:,3) = ((L3(:,1) - L3(:,2)) ./ L3(:,1)) .* 100;
dlhAtIcandRad(:,3) = dlhAtIcandRad(:,1) - dlhAtIcandRad(:,2);

dlhAtIcandRad(:,4) = ((dlhAtIcandRad(:,1) - dlhAtIcandRad(:,2)) ./ dlhAtIcandRad(:,1)) .* 100;

meandlhAtIcandRad = mean(dlhAtIcandRad(:,3));
mean([dlhradarAtIcecore;l2_IC-l1_IC])

