close all
clear
clc
P = ThisIsAn_ApRES_Script(mfilename('fullpath'));
%% Dir
% InfoDir = fullfile(P.data,'QuadPolarimetric_PointsInfo.csv');
nm = 'p0';
DataDir = fullfile(P.data,'radar',nm);
DatFiles = dir(fullfile(DataDir,'*.dat'));
maxRange = 701;

whichfiles = [1,2,3,4];
for i = 1:length(whichfiles)
    filePath = fullfile(DataDir,DatFiles(whichfiles(i)).name);
    DAT(i) = FUNC_ReadMonsterFile(filePath,maxRange,1);
end
Z = DAT.ZT;
%%
HH = DAT(1).Signal;
HV = DAT(2).Signal;
VH = DAT(3).Signal;
VV = DAT(4).Signal;

pwrHH=20.*log10(abs(HH));
pwrHV=20.*log10(abs(HV));
pwrVH=20.*log10(abs(VH));
pwrVV=20.*log10(abs(VV));

%%
dA = 1;
ao = 0:dA:179; 
f = 3.0000e+08;
C_DepthWin = maxRange * 0.1;
C_ConvWin = maxRange * 0.1;
DenoisingFlag.PA = [  "1", "MovingAverage"  , string(maxRange*0.01) ;
                      "0", "Conv1D"         , string(maxRange*0.1) ;
                      "2", "Conv2D"         , string(maxRange*0.01) ;
                      "0", "DenoisePCA"     , string(1)];
DenoisingFlag.PD = [  "1", "MovingAverage"  , string(maxRange*0.01) ;
                      "0", "Conv1D"         , string(maxRange*0.01) ;
                      "0", "Conv2D"         , string(maxRange*0.01) ;
                      "0", "DenoisePCA"     , string(1)];
                  
%
[HH,VH,HV,VV] = QuadpoleSynthesizer(HH,VH,HV,VV,ao,0);
ObsDta = CLASS_S2P.Signal2Param(HH,VH,HV,VV,Z,ao,f,C_DepthWin,C_ConvWin,DenoisingFlag,"radar");
%%

% rw = 1;
% col = 2;
% 
% fig = SetFigureSize(0,0,0.30,0.5);
% 
% subplot(rw,col,1)
% 
% plot(mean([pwrHH pwrHV pwrVH],2),Z,'-k','LineWidth',1.2)
% ylim([0 max(Z)])
% set(gca,'YDir','reverse')
% xlabel('Power [dB]')
% ylabel(gca,'Depth [m]')
% set(gca,"FontSize",16)
% 
% subplot(rw,col,2)
% plot(mean(ObsDta{13},2),Z,'-k','LineWidth',1.2)
% hold on
% plot([0.4 0.4],[Z(1) Z(end)],'--r')
% set(gca,'YDir','reverse')
% xlabel('Coherence magnitude [-]')
% ylim([0 max(Z)])
% xlim([0 1])
% set(gca,"FontSize",16)

%%
fig = SetFigureSize(0,0,0.3,0.3);
yyaxis left
plot(Z,mean([pwrHH pwrHV pwrVH],2),'LineWidth',1.2)
xlim([0 max(Z)])
ylabel('Power [dB]')
xlabel(gca,'Depth [m]')
set(gca,"FontSize",16)

yyaxis right
plot(Z,mean(ObsDta{13},2),'LineWidth',1.2)
hold on
plot([Z(1) Z(end)],[0.4 0.4],'--k')

rectangle('Position',[Z(1) 0.4 Z(end) 1],...
    'FaceColor',[0 1 0 0.1],'EdgeColor','g','LineWidth',0.1)

rectangle('Position',[Z(1) 0 Z(end) 0.4],...
    'FaceColor',[1 0 0 0.1],'EdgeColor','r','LineWidth',0.1)

xlim([0 max(Z)])
ylabel('Coherence magnitude [-]')
ylim([0 1])
set(gca,"FontSize",16)

%%
fig.InvertHardcopy = 'off';
% print(fig,"HIR_Power&MeanCohMag_"+nm+".png",'-dpng','-r300');
%%
function f = SetFigureSize(ss1,ss2,w,h)
    f = figure;
    set(f,'Color',[1 1 1]);
    set(f, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]); % full screen figure
    set(f, 'Units', 'centimeters');
    scrn = get(f, 'OuterPosition'); % get the size of the screen in CM

    wdt = scrn(3) * w;
    hgt = scrn(4) * h;

    s1 = scrn(3) * ss1;
    s2 = scrn(4) * ss2;

    set(f, 'OuterPosition', [s1, s2, wdt, hgt]); % change the figure size to the new size
    set(f, 'Units', 'Normalized');
end