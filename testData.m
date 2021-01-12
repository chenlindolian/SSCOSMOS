clc,clear all
addpath('~/Desktop/SSCOSMOS/sharecode/SSCOSMOS/SimulationData/');
%% load in QSMParams
OriNum = 3;
FileStrAppArray = cell(OriNum, 1);
for OriInd = 1:OriNum
    FileStrAppArray{OriInd} = ['_Ori', num2str(OriInd)];
end

QSMParams.OriNum = OriNum;

QSMParams.filenamePhase = cell(QSMParams.OriNum, 1);
for OriInd = 1:OriNum
    QSMParams.filenamePhase{OriInd} = ['EveQSM_SimForward', FileStrAppArray{OriInd}, '_v2.1.mat']; 
end

QSMParams.dataType = 0;         % 0 for Simulation, 1 for GRE phase data
QSMParams.fov = [216, 180, 180];
QSMParams.sizeVol = [216, 180, 180];
QSMParams.sizeRecAll = [216, 180, 180];
QSMParams.voxSize = QSMParams.fov./QSMParams.sizeVol;

kernelchoice = 2;           %  1 for Lap SSCOSMOS kernel, 2 for SMV SSCOSMOS kernel
switch kernelchoice       
    case 1         
        QSMParams.kernel = 'Lap'; 
        QSMParams.radiusArray = 1; 
        alpha = 0;                      
    case 2
        QSMParams.kernel = 'SMV';      
        QSMParams.radiusArray = 1:3;    % for SMV only, if using radius of 1 should be equivalent to Lap
end

maxit = 300;
tol = 1e-6;     

filenameSave = [pwd, 'SSCOSMOS/SSCOSMOS_Ori_', num2str(QSMParams.OriNum), '_', QSMParams.kernel,'_',num2str(length(QSMParams.radiusArray))];
        

switch QSMParams.kernel
    case 'Lap'
         [chi_mo, itn, mask_eval] = SScosmos_Lap_lsmr(QSMParams, maxit, tol);
         
    case 'SMV'
         [chi_mo, itn, mask_eval] = SScosmos_SMV_lsmr(QSMParams, maxit, tol);

end


% Estimate performance
S = load(QSMParams.filenamePhase{1});
if ~exist('mask_eval', 'var')
    maskErode = imerode3(S.BrainMask, 1);   
else
    maskErode = mask_eval;
end

% Target
N = S.Params.sizeVol;
chi_tissue_ref = zeros(N);   % chi_tissue referenced
chi_tissue_ref(maskErode>0) = S.chi_tissue(maskErode>0) - mean((S.chi_tissue(maskErode>0)));

chi_mo_mean = zeros(N);
chi_mo_mean(maskErode==1) = chi_mo(maskErode==1) - mean(chi_mo(maskErode==1));
    
rmse_Lcosmos = 100 * norm((chi_mo_mean(maskErode>0) - chi_tissue_ref(maskErode>0))) / norm(chi_tissue_ref(maskErode>0));
disp(['RMSE cosmos = ', num2str(rmse_Lcosmos)])

ssim_Lcosmos = compute_ssim(chi_mo_mean, chi_tissue_ref);
disp(['SSIM cosmos = ', num2str(ssim_Lcosmos)])

hfen_Lcosmos = compute_hfen(chi_mo_mean, chi_tissue_ref);
disp(['HFEN cosmos = ', num2str(hfen_Lcosmos)])

% save results
   
save([filenameSave, '.mat'], 'chi_mo',  'itn',  'maskErode');

saveNII(chi_mo_mean, filenameSave, QSMParams, 1);   

dispFlag = 1;
if dispFlag == 1
     plot3plane(permute(chi_mo_mean, [2,1,3]), 113, 90, 72, -0.15, 0.15, 0, 1);                       
     plot3plane(permute(chi_mo_mean-chi_tissue_ref, [2,1,3]), 113, 90, 72, -0.15, 0.15, 0, 1);         
end

