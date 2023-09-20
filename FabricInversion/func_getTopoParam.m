function [xSEBETK,pRESpoints,H,SMB,SVsia,SVsat] = func_getTopoParam(P,RP)
dtapath = fullfile(P.Project,'data');
opts = spreadsheetImportOptions("NumVariables", 4);
opts.VariableNames = ["Distance", "Lat", "Lon", "Par"];
opts.VariableTypes = ["double", "double", "double", "double"];

opts.Sheet = "SurfaceElevation";
REMA8m = readtable(fullfile(dtapath,'HIR_ExtendedProfile.xlsx'), opts, "UseExcel", false);
REMA8m(1,:) = [];

opts.Sheet = "BedElevation_BedMachineV3";
BedElv = readtable(fullfile(dtapath,'HIR_ExtendedProfile.xlsx'), opts, "UseExcel", false);
BedElv(1,:) = [];
%%
xRP = 11.978;
for i = 1:length(RP)-1
    xRP(i+1) = xRP(i) + abs(RP(i+1) - (RP(i)));
end
xp7 = xRP(1);
lpRES = xRP(end) - xRP(1);
xp14 = xp7+lpRES;

% Surface elevation directly from REMA
SEpp = [table2array(REMA8m(:,1))./1000 table2array(REMA8m(:,end))];

% bed elevation directly from BMV3
BEpp = [table2array(BedElv(:,1))./1000 table2array(BedElv(:,end))];
BEpp = BEpp(1:end-1,:);

xq = SEpp(:,1);

% x values from REMA
xSEBETK(:,1) = xq;

% surface elevation from REMA
xSEBETK(:,2) = SEpp(:,2);

% Interpolated Bed Elevation from BMV3
xSEBETK(:,3) = interp1(BEpp(:,1),BEpp(:,2),xq);

% Calculate the Thickness from the difference between REMA and Interpolated bed elevation from BMV3 thickness
xSEBETK(:,4) = abs( xSEBETK(:,2) - xSEBETK(:,3) );

pRESpoints(1,:) = RP;
pRESpoints(2,:) = xRP;
for i = 1:length(xRP)
    [~,ii] = min(abs(xSEBETK(:,1)-xRP(i)));
    pRESpoints(3,i) = xSEBETK(ii,2);
    pRESpoints(4,i) = xSEBETK(ii,3);
    pRESpoints(5,i) = xSEBETK(ii,4);
end

[~,im1] = min(abs(xSEBETK(:,1)-xp7));
[~,im2] = min(abs(xSEBETK(:,1)-xp14));
x_pRES = xSEBETK(im1:im2,1);
se_pRES = xSEBETK(im1:im2,2);
be_pRES = xSEBETK(im1:im2,3);
tk_pRES = xSEBETK(im1:im2,4);

% meanSE = round(mean(se_pRES),0);
% meanBE = round(mean(be_pRES),0);
meanTK = round(mean(tk_pRES),0);

% Thickness at ice divide which by chance laso is equal to the mean
% thickness along the pRES profile (549 m)
H = round(meanTK,0);

pRESpoints(6,:) = pRESpoints(1,:).*1000./H; % normlized distance

% SE = pRESpoints(3,:);
% BE = pRESpoints(4,:);
%%
opts.Sheet = "SMB";
SMB = readtable(fullfile(dtapath,'HIR_ExtendedProfile.xlsx'), opts, "UseExcel", false);
SMB(1,:) = [];

opts.Sheet = "SurfaceVelocity(SIA)";
SVsia = readtable(fullfile(dtapath,'HIR_ExtendedProfile.xlsx'), opts, "UseExcel", false);
SVsia(1,:) = [];

opts.Sheet = "SurfaceVelocity(SAT)";
SVsat = readtable(fullfile(dtapath,'HIR_ExtendedProfile.xlsx'), opts, "UseExcel", false);
SVsat(1,:) = [];