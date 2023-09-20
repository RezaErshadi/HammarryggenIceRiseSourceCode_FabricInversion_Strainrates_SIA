clear;
close all;
% clc
P = ThisIsAn_ApRES_Script(mfilename('fullpath'));
%%
DataDir = fullfile(P.Project,'results','InvRes','HIR_Profile_Estimated.mat');
load(DataDir);
RP = HIR_Prof.RadarPos./1000;
BEfromRADAR = HIR_Prof.SBR_Elv(3,:); % Bed elevation at each site from pRES data

[xSEBETK,pRESpoints,H,Smassbalance,SvelocitySIA,SvelocitySAT] = func_getTopoParam(P,RP);
xq = xSEBETK(:,1);
xp7 = pRESpoints(2,1);
xp14 = pRESpoints(2,end);
lpRES = abs(xp14 - xp7);

[~,im1] = min(abs(xSEBETK(:,1)-xp7));
[~,im2] = min(abs(xSEBETK(:,1)-xp14));
x_pRES = xSEBETK(im1:im2,1);
se_pRES = xSEBETK(im1:im2,2);
be_pRES = xSEBETK(im1:im2,3);
tk_pRES = xSEBETK(im1:im2,4);

meanSE = round(mean(se_pRES),0);
meanBE = round(mean(be_pRES),0);

xlbl(1,:) = [-4.5 -3 -2 -1 0 1 2 3 4.5];
xlbl(2,:) = round([-4.5 -3 -2 -1 0 1 2 3 4.5].*H,0)./1000;
xlbl(3,:) = xlbl(2,:) + pRESpoints(2,8);

SVsia = [table2array(SvelocitySIA(:,1))./1000 table2array(SvelocitySIA(:,end))];
SVsat = [table2array(SvelocitySAT(:,1))./1000 table2array(SvelocitySAT(:,end))];
SMB = [table2array(Smassbalance(:,1))./1000 table2array(Smassbalance(:,end))];
aSV = lm(SVsia);
aSMB = lm(SMB);
%%
fig1 = CLASS_FixedPlot.SetFigureSize(0,0,0.2,0.9);
redTrans = 0.1; 
rw = 4;
clm = 1;
% -------------------------
subplot(rw,clm,[1:2])
plot(pRESpoints(2,:),pRESpoints(3,:),'or','MarkerFaceColor','r')
hold on
plot(pRESpoints(2,:),BEfromRADAR,'xr')
% plot(xSEBETK(:,1),xSEBETK(:,2),'-k','LineWidth',6)
% hold on
% plot(xSEBETK(:,1),xSEBETK(:,3),'-k','LineWidth',6)
c = jet(length(xq));
nTK = rescale(xSEBETK(:,4),1,length(xq));
for i = 1:length(xq)
    plot([xq(i) xq(i)],[xSEBETK(i,2) xSEBETK(i,3)],'Color',c(round(nTK(i),0),:),'LineWidth',3)
end
colormap(gca,c);
% colorbar(gca,'southoutside')
clb = colorbar(gca,"southoutside");
clb.Position = [0.167363534226055 0.599013968775677 0.261415535541386 0.0144487869992782];
caxis(gca,[min(xSEBETK(:,4)) max(xSEBETK(:,4))])

clb.Ticks = round(linspace(min(xSEBETK(:,4)),max(xSEBETK(:,4)),5),0) ; %Create 8 ticks from zero to 1
% clb.TickLabels = num2cell(1:8) ;

plot(pRESpoints(2,:),pRESpoints(3,:),'or','MarkerFaceColor','r')
plot(pRESpoints(2,:),BEfromRADAR,'xr')

legend(gca,'pRES points','Detected bed','Location','northwest')

xlim([min(xq) max(xq)])
set(gca,'xticklabel',{[]})
ylabel('[m.a.s.l.]')
set(gca,'FontSize',18)
rectangle('Position',[xp7 -600 lpRES 1000],...
    'FaceColor',[1 0 0 redTrans],'EdgeColor','r','LineWidth',0.1)
% -------------------------
subplot(rw,clm,[3])
plot(SVsia(:,1),SVsia(:,2),'.-b','linewidth',2)
hold on
plot(SVsat(:,1),SVsat(:,2),'.r','linewidth',2,'MarkerSize',10)
rectangle('Position',[xp7 aSV(1) lpRES diff(aSV)],...
    'FaceColor',[1 0 0 redTrans],'EdgeColor','r','LineWidth',0.1)
ylabel('[m/a]')
xlim([min(xq) max(xq)])
ylim(aSV)
set(gca,'xticklabel',{[]})
% title('Surface velocity')
legend(gca,'SIA','Satellite','Location','northwest')
set(gca,'FontSize',18)
% -------------------------
subplot(rw,clm,[4])
plot(SMB(:,1),SMB(:,2),'.-b','linewidth',2)
hold on
plot([xp7 xp14],[400 400],'--k','linewidth',2)
rectangle('Position',[xp7 aSMB(1) lpRES diff(aSMB)],...
    'FaceColor',[1 0 0 redTrans],'EdgeColor','r','LineWidth',0.1)
ylabel('[mm w.e. a^{-1}]')
% xlabel('Length [km]')
% xticks(linspace(SurveElevation(1,1),SurveElevation(end,1),5))
xlim([min(xq) max(xq)])
ylim(aSMB)
% title('Surface mass balance')
legend(gca,'Lenaerts and others, 2014','Cavitte and others., 2022','Location','west')
set(gca,'FontSize',18)
% -------------------------
% subplot(rw,clm,[5])
% plot(x_pRES,se_pRES,'-k','LineWidth',3)
% hold on
% plot(x_pRES,be_pRES,'-k','LineWidth',3)
% c = c(im1:im2,:);
% nTK = rescale(tk_pRES,1,length(x_pRES));
% % for i = 1:length(x_pRES)
% %     plot([x_pRES(i) x_pRES(i)],[se_pRES(i) be_pRES(i)],'Color',c(round(nTK(i),0),:),'LineWidth',7)
% % end
% % plot(pRESpoints(2,:),BEfromRADAR,'xr')
% % colormap(gca,c);
% plot([min(x_pRES)-0.3 max(x_pRES)+0.3],[meanSE meanSE],'.-b','LineWidth',2)
% plot([min(x_pRES)-0.3 max(x_pRES)+0.3],[meanBE meanBE],'.-b','LineWidth',2)
% plot(pRESpoints(2,:),pRESpoints(3,:),'or','MarkerFaceColor','r')
% plot(pRESpoints(2,:),BEfromRADAR,'or','MarkerFaceColor','r')
% xlim([min(x_pRES)-0.3 max(x_pRES)+0.3])
% xticks(gca,xlbl(3,:))
% xticklabels(gca,xlbl(1,:))
% ylabel('[m.a.s.l.]')
% set(gca,'FontSize',18)
%%
% fig1.InvertHardcopy = 'off';
% print(fig1,"TIR_Profile_QGIS.png",'-dpng','-r300');
%%
function a = lm(aa)
l = max(aa(:,2)) - min(aa(:,2));
a = [min(aa(:,2))-(l*0.05) max(aa(:,2))+(l*0.05)];
end