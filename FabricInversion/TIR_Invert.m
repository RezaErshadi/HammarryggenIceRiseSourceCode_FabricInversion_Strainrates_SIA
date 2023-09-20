clear
close all
clc
[PP,SavedPath,Prjcts,ps] = FUNC_ApRES_PathFix;
tm = string(datetime(now,'ConvertFrom','datenum'));
tic; % Start Timer
%% -------------------------------- Input Parameters
    %----------------- Reading data input parameters -----------------***** 
ProjectName = "TIR_Brussels"; % Project's folder name
SiteName = []; % leave empty for batch mode
Save_Inversion = 1;
EigenValuesFile = [];
maxRange = 401;                         % |*****CHANGE THIS*****|
BedRange = [];                   % |*****CHANGE THIS*****|
% Show the plot (1) or not (0)?
FigVis = 1;
% Azimuthal spacing (synthesizing data)
dA = 1;
% Azimuthal Vector
ao = 0:dA:179; 
% Depth window which smooths coherence
%(averaging)
C_DepthWin = maxRange * 0.1;           % |*****YOU MIGHT NEED TO CHANGE THIS*****|
%(1D convolution)
C_ConvWin = maxRange * 0.1;            % |*****YOU MIGHT NEED TO CHANGE THIS*****|
% Smoothing PowerAnomaly (PA) and PhaseDifference (PD)
% First column digit shows which method applies first
% The last column values are the factor of smoothing for that method
% For PCA, the value is the number of Pricipal components
DenoisingFlag.PA = [  "1", "MovingAverage" , string(maxRange*0.05) ;
                      "0", "Conv1D" , string(maxRange*0.1) ;
                      "2", "Conv2D" , string(maxRange*0.05) ;
                      "0", "DenoisePCA" , string(1)];
DenoisingFlag.PD = [  "1", "MovingAverage" , string(maxRange*0.05) ;
                      "0", "Conv1D" , string(maxRange*0.01) ;
                      "0", "Conv2D" , string(maxRange*0.01) ;
                      "0", "DenoisePCA" , string(1)];
    %*****----------------- Inversion input parameters -----------------*****              
ResZmdl = 10;
resInv = ResZmdl; % Inversion resolution can not be finer than model resolution
FabricOrientationApproach = "Manual";
switch FabricOrientationApproach
    case "Auto"
        MHB = [];
    case "SemiAuto"
        MHB = [ 100 ; 300 ; maxRange ];
    case "Manual"
        MHB = [ [maxRange] , [90] ];
end
options1 = optimoptions('fmincon','Display','iter-detailed');
options1.UseParallel = true;
options1.OptimalityTolerance = 1e-6;
options1.FunctionTolerance = 1e-6;
options1.ConstraintTolerance = 1e12;
InvOrder = "r&v1"; % r&v1 , v1&r
% v1 parameters
v1CF = 5; % v1 cost function
v1ST = 1e-6; % v1 step tolerance
v1AM = "Pointwise"; % v1 Approach Method: Legendre, Pointwise
if v1AM == "Legendre" 
    v1NLC= 30; % v1 Number of Legendre Coefficients
else
    v1NLC = [];
end
% r parameters
rCF = 5; % r cost function
rST = 1e-6; % r step tolerance
rAM = "Pointwise"; % r Approach Method: Legendre, Pointwise
if rAM == "Legendre"
    rNLC = 10; % r Number of Legendre Coefficients
else
    rNLC = [];
end
%% Pack the input parameters
if isempty(SiteName)
    SiteName = Func_GetAllTheSites(PP,ps,ProjectName);
end
prjcDir = strcat(PP,ps,ProjectName);
InfoDir = strcat(prjcDir,ps,'_Info');
inputsName = ["PP","SavedPath","Prjcts","ps","tm"...
"ProjectName","InfoDir","EigenValuesFile","FigVis",...
"maxRange","BedRange","dA","ao",...
"C_DepthWin","C_ConvWin","DenoisingFlag",...
"ResZmdl","resInv","FabricOrientationApproach","MHB","options1","InvOrder",...
"v1ST","v1AM","v1NLC","v1CF"...
"rST","rAM","rNLC","rCF"];
for i = 1:length(inputsName)
    eval(strcat('InputParameters.',inputsName(i),' = ',inputsName(i),';'));
end
%%
for i = 1:size(SiteName,1)
    close all
    [InversionOutput,SiteNumber,fg,ax] = FUNC_RunInversion(InputParameters,SiteName(i));
    if Save_Inversion ~= 0
        SavedPathPrjc = strcat(SavedPath,ps,ProjectName);
        if ~exist(SavedPathPrjc, 'dir')
            mkdir(SavedPathPrjc)
        end
        LabelName = strcat(SavedPathPrjc,ps,string(SiteNumber),'_',ProjectName,'_',SiteName(i),'_',tm);
        save(LabelName+"_InversionResults.mat",'InversionOutput');
        print(fg,LabelName+"_InversionFigure.png",'-dpng','-r300');
    end
end





























