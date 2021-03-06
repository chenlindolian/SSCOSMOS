function  [chi_mo, itn, mask_eval] = SScosmos_Lap_lsmr(QSMParams, maxit, tol)
% Quantitative Susceptibility Mapping (QSM) using single-step multiple orientations

%% Author: Lin Chen, Xu Li
% Affiliation: Radiology @ JHU
% Email address: lchen151@jhmi.edu
%                xuli@mri.jhu.edu
% 
% % Quantitative Susceptibility Mapping (QSM) using Laplacian kernel based SSCOSMOS
% Input:  
%       QSMParams includes fields:
%       OriNum: Orientation Number      
%       filenamePhase:  cell array containing the filenames of Coregistered
%       Raw phase (coregistered by real and imaginary part) for human brain
%       filenameMag: cell array containing the filenames of Coregistered magnitude
%       filenameCoregParams: cell array containing the filenames of scan
%       parameters (fov, sizeVol, VoxSize, TAng: field of view, size of Volume, voxel size, 
%       coregisterd angles with orientation angles)
%       filenameMask: cell array containing the filenames of Coregistered brain mask           
%       maxit: maximum number of iteration
%       tol: convergence tolerence for LSMR
%    
% Output:
%     chi_mo: isotropic susceptibility map
%     itn: total iteration number
%     mask eval: Brain mask 

OriNum = QSMParams.OriNum;
ParamsArray = cell(OriNum, 1);

N = QSMParams.sizeVol;
voxel_size = QSMParams.voxSize;

fittingDataArray = zeros([N, OriNum]);

kernel2 = generate_LapKernel(N);
kernel2_conj = conj(kernel2);

mask_eval = ones(N);
for OriInd = 1:OriNum
   switch QSMParams.dataType
        case {0} % for simulation
             S = load(QSMParams.filenamePhase{OriInd});
             GREPhase = S.GREPhase;
           
             GREMag = S.GREMag;
             BrainMask = imerode3(S.BrainMask, 1);      % boundary inconsistancy   
             mask_unrelyPhase = create_mask_unrelyPhase(S.GREPhase, pi/3);   % mask threshold may need to change
 
             BrainMask = BrainMask.*(BrainMask - mask_unrelyPhase);
             BrainMask = imfill3(BrainMask);
             Mask = BrainMask; % reliable mask of phase
                
        case {1} % for huamn brain
             S = load(QSMParams.filenameCoregParams{OriInd});
             
             nii = load_untouch_nii(QSMParams.filenamePhase{OriInd});
             GREPhase = nii.img;
             GREPhase = permute(GREPhase,[2,1,3,length(size(GREPhase))]); 
            
             nii = load_untouch_nii(QSMParams.filenameMask{OriInd});
             Mask = nii.img;
             Mask = permute(Mask,[2,1,3]);
             
             nii = load_untouch_nii(QSMParams.filenameMag{OriInd});
             GREMag = nii.img;                       
             GREMag = permute(GREMag, [2, 1, 3:length(size(GREMag))]);  
              
   end

 mask_eval = mask_eval & Mask;

 Params = S.Params;
    
 Params.D = ifftshift(conv_kernel_rot_c0(Params,Params.TAng));
    
 Params.C = 1./(2*pi*Params.gamma*Params.B0)*1e6;
  
       
   if Params.nEchoes > 1
        temp = zeros(N);
        weight = AverageEchoWeight(GREMag,Params.TEs);
         for echoii = 1:Params.nEchoes 
            selectedEcho = Params.echoNums(echoii);
            temp = temp + LapPhase(GREPhase(:,:,:,selectedEcho), kernel2)/(Params.TEs(selectedEcho)).*weight(:,:,:,selectedEcho);
         end
        temp = temp./(sum(weight,4)+0.00001); % averaging
        fittingDataArray(:,:,:,OriInd) = temp.*Params.C; 
   else
        fittingDataArray(:,:,:,OriInd) = LapPhase(GREPhase, kernel2)/Params.TEs.*Params.C;                        
   end
   
 switch QSMParams.dataType
    case {0}
        Params.Weight = ones(N);
    case {1}    
        QSMParams.WeightEcho = 1:Params.nEchoes;
    
        Params.Weight = sqrt(sum(single(GREMag(:,:,:,QSMParams.WeightEcho)).^2, 4));

        Params.Weight = Params.Weight.*Mask;
        Params.Weight = Params.Weight./mean(tovec(Params.Weight(Mask > 0)));
    
 end  
    ParamsArray{OriInd} = Params;
end

%% using lsmr method 

VoxNum = prod(N);
b = zeros(OriNum*VoxNum, 1, 'single');

outSMV = mycreate_SMVkernel(mask_eval, 1, 1, 1, N, voxel_size);
mask_eval = sum(outSMV.SMV_mask, 4);

for OriInd = 1:OriNum     
    Params = ParamsArray{OriInd}; 
    temp = mask_eval.*Params.Weight.*fittingDataArray(:,:,:,OriInd); 
    b(((OriInd - 1)*VoxNum+1) : OriInd*VoxNum) = temp(:);
end

disp('solving QSM sscosmos inverse problem using lsmr ...');
tic
lambda = 0;
atol = tol;
btol = tol;
conlim = 1e+8;
itnlim = maxit;
[chi_mo, istop, itn, normr, normAr, normA, condA, normx]= lsmr(@afun,b,lambda,atol,btol,conlim,itnlim);
toc

chi_mo = real(reshape(chi_mo, N));

function y = afun(x,transp_flag)
   if strcmp(transp_flag,'transp')      % y = A'*x
       y = zeros(VoxNum, 1, 'single');
       x = single(x); 
       for orient_i = 1:OriNum          
            Params = ParamsArray{orient_i};
            temp = x(((orient_i - 1)*VoxNum+1) : orient_i*VoxNum);            
            
            temp = reshape(temp, N).*mask_eval; 
            
            temp = ifftn(fftn(temp).*kernel2_conj.*Params.D);
            y = y + temp(:);             
       end
              
   elseif strcmp(transp_flag,'notransp') % y = A*x     
       y = zeros(OriNum*VoxNum, 1);
       x = single(x);                               % change to single format
                         
       x = reshape(x, N); 
        
       for orient_i = 1:OriNum
            Params = ParamsArray{orient_i}; 
            
            temp = mask_eval.*ifftn(fftn(x).*kernel2.*Params.D);                       
            
            y(((orient_i - 1)*VoxNum+1) : orient_i*VoxNum) = temp(:);           
       end
       
   end
end

function output = LapPhase(input, L)
    % laplacian of phase using sine/cosine, input is phase with wraps
    PhaseSin = sin(input);
    PhaseCos = cos(input);    
    output  = PhaseCos.*ifftn(L.*fftn(PhaseSin)) - PhaseSin.*ifftn(L.*fftn(PhaseCos));      
end

end