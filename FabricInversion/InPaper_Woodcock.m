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

% k_R(k_R<0.2) = 0; % UAG
% k_R(k_R>=0.2 & k_R<0.5) = -1; % G
% k_R(k_R>=0.5 & k_R<0.98) = -2; % G
% k_R(k_R>=0.98 & k_R<=1.02) = -3; % GCTZ
% k_R(k_R>1.02 & k_R<=2) = -4; % C
% k_R(k_R>2 & k_R<=5) = -5; % C
% k_R(k_R>5) = -6; % UAC
% k_R = -1*k_R;

Z_R = HIR_Prof.Zinv;
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


%% Create WOODCOCK 1977
lnx = 0:0.01:6; % ln(lambda2 / lambda1)
lny = 0:0.01:6; % ln(lambda3 / lambda2)
[lnX,lnY] = meshgrid(lnx,lny);
K = lnY./lnX;
C = (K+1).*lnX; % ln(lambda3 / lambda1)
a = exp(C);
b = exp(lnY);
c = exp(lnX);
l1 = 1./(a+c+1);
l2 = c./(a+c+1);
l3 = a./(a+c+1);
l2l1 = l2-l1;
l3l1 = l3-l1;
l3l2 = l3-l2;

%% Plot WOODCOCK 1977 FIGURE1
f2 = CLASS_FixedPlot.SetFigureSize(0,0,0.5,0.65);
h = pcolor(lnX,lnY,l2l1);
fsz = 20;
set(h, 'EdgeColor', 'none');
clmp = getPyPlot_cMap('Greys',20);
colormap((clmp))
% set(gca,'ColorScale','log')
caxis([0 0.5])
cb = colorbar;
cb.Label.String = '\lambda_2-\lambda_1 [-]';
% cb.Ticks = [2e-1 1e0 5e0];
% cb.TickLabels = {'UAG','GCTZ','UAC'};
hold on
% contour(lnX,lnY,l2l1,[0:0.05:0.5],'--g','LineWidth',2);
contour(lnX,lnY,K,[0 0.2 0.5 1 2 5 inf],'--b','LineWidth',2);
% clabel(cc,hh,'FontSize',fsz,'Color','g','LabelSpacing',750)
contour(lnX,lnY,C,[0:2:10],'-r','LineWidth',2);
% clabel(cc,hh,'FontSize',fsz,'Color','g','LabelSpacing',750)
contour(lnX,lnY,l3l2,[0:0.1:0.9 0.95 0.975],':','LineColor','g','LineWidth',2);

plot(lns2s1_IC,lns3s2_IC,'ok','MarkerFaceColor','k','MarkerSize',10)
plot(lns2s1_R(5:26),lns3s2_R(5:26),'sk','MarkerFaceColor','g','MarkerSize',10)

% plot(rotated_ellipse_IceCore(1,:),rotated_ellipse_IceCore(2,:),'w','LineWidth',2);
% plot(rotated_ellipse_Radar(1,:),rotated_ellipse_Radar(2,:),'b','LineWidth',2);

pbaspect([1 1 1])
xlabel('ln(\lambda_2/\lambda_1)')
ylabel('ln(\lambda_3/\lambda_2)')
set(gca,'FontSize',fsz)
box on


legend('','K','C','\lambda_3-\lambda_2','Meas.','Est.','Location','northeast')

f2.InvertHardcopy = 'off';
% print(f2,"AllIn1CompleteWoodcock1977_EllipseData.png",'-dpng','-r300');


