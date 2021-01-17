function  [chi_mo, itn, mask_eval] = SScosmos_SMV_lsmr(QSMParams, maxit, tol)
%       or   function  [chi_mo, itn, mask_eval, freqshift] = SScosmos_SMV_lsmr(QSMParams, maxit, tol)   
%% Author: Lin Chen, Xu Li
% Affiliation: Radiology @ JHU
% Email address: lchen151@jhmi.edu
%                xuli@mri.jhu.edu
% 
% % Quantitative Susceptibility Mapping (QSM) using SMV based SSCOSMOS
%
% Input:  
%       QSMParams includes fields:
%       OriNum: Orientation Number
%       radiusArray: Radius of SMV kernels         
%       filenamePhase:  cell array containing the filenames of Coregistered Unwrapped phase 
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
%     mask eval: Brain mask after SMV kernel

OriNum = QSMParams.OriNum;
ParamsArray = cell(OriNum, 1);
outsmvArray = cell(OriNum, 1);

N = QSMParams.sizeVol;
voxel_size = QSMParams.voxSize;

fittingDataArray = zeros([N, OriNum]);
Maskarray = zeros([N, OriNum]);

for OriInd = 1:OriNum
    switch QSMParams.dataType
        case {0} % for simulation
             S = load(QSMParams.filenamePhase{OriInd});
             GREPhase = S.PhaseUnwrap; % or use Phase unwrapped with different phase unwrapping methods
             BrainMask = imerode3(S.BrainMask, 1);      % boundary inconsistancy   
             mask_unrelyPhase = create_mask_unrelyPhase(S.GREPhase, pi/3);   % mask threshold may need to change

             BrainMask = BrainMask.*(BrainMask - mask_unrelyPhase);
             BrainMask = imfill3(BrainMask);
             Mask = BrainMask; % reliable mask of phase             Mask = S.BrainMask;
             GREMag = S.GREMag;
             
        case {1} % for huamn brain
             nii = load_untouch_nii(QSMParams.filenamePhase{OriInd});
             GREPhase = nii.img;
             GREPhase = permute(GREPhase,[2,1,3,length(size(GREPhase))]); 
            
             nii = load_untouch_nii(QSMParams.filenameMask{OriInd});
             Mask = nii.img;
             Mask = permute(Mask,[2,1,3]);
             
             nii = load_untouch_nii(QSMParams.filenameMag{OriInd});
             GREMag = nii.img;                       
             GREMag = permute(GREMag, [2, 1, 3:length(size(GREMag))]);  
             
             S = load(QSMParams.filenameCoregParams{OriInd});
    end
  
    Params = S.Params;
    Params.D = ifftshift(conv_kernel_rot_c0(Params,Params.TAng));
    Params.C = 1./(2*pi*Params.gamma*Params.B0)*1e6;
    
    min_radius = min(QSMParams.radiusArray);
    max_radius = max(QSMParams.radiusArray);
    step_size_radius = QSMParams.radiusArray(2) - QSMParams.radiusArray(1);                                                
    outSMV = mycreate_SMVkernel(Mask, min_radius, max_radius, step_size_radius, N, voxel_size); % creat SMV kernels

    numKernel = size(outSMV.SMV_kernel,4);
    Maskarray(:,:,:,OriInd) = sum(outSMV.SMV_mask,4);
    outSMV.SMV_kernel_conj=conj(outSMV.SMV_kernel);
       
    if Params.nEchoes > 1
        temp = zeros(N);
        weight = AverageEchoWeight(GREMag,Params.TEs);
        for echoii = 1:Params.nEchoes 
            GREphase = GREPhase(:,:,:,echoii);
            Fphi = fftn(GREphase);
            phasetemp = zeros(N);
            for kk = 1:numKernel
                phasetemp = phasetemp + outSMV.SMV_mask(:,:,:,kk).*ifftn(outSMV.SMV_kernel(:,:,:,kk).*Fphi);
            end
            temp = temp + phasetemp/(Params.TEs(echoii)).*weight(:,:,:,echoii);
        end
        temp = temp./(sum(weight,4)+0.00001); % averaging
        fittingDataArray(:,:,:,OriInd) = temp.*Params.C; 
    else
        phasetemp = zeros(N);     
        Fphi = fftn(GREPhase);
        for kk = 1:numKernel
          phasetemp = phasetemp +  outSMV.SMV_mask(:,:,:,kk).* ifftn(outSMV.SMV_kernel(:,:,:,kk).*Fphi);  
        end
        fittingDataArray(:,:,:,OriInd) = phasetemp/Params.TEs.*Params.C;                         
    end
   
    % weight
    switch QSMParams.dataType
      case {0}
         Params.Weight = ones(N);
      case {1}
         QSMParams.WeightEcho = 1:Params.nEchoes;

         Params.Weight = sqrt(sum(single(GREMag(:,:,:,QSMParams.WeightEcho)).^2, 4));

         Params.Weight = Params.Weight.*Maskarray(:,:,:,OriInd);
         Params.Weight = Params.Weight./mean(tovec(Params.Weight(Maskarray(:,:,:,OriInd)> 0)));

    end
    
    ParamsArray{OriInd} = Params;
    outsmvArray{OriInd} = outSMV;
end

mask_eval = sum(Maskarray, 4) >= 3;

%% using lsmr method 
VoxNum = prod(N);
b = zeros(OriNum*VoxNum, 1, 'single');

for OriInd = 1:OriNum     
    Params = ParamsArray{OriInd}; 
    temp = Params.Weight.*fittingDataArray(:,:,:,OriInd); 
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
      x = single(x);                                     % change to single format
      for orient_i = 1:OriNum          
        Params = ParamsArray{orient_i};
        outSMV = outsmvArray{orient_i};
        tempin = x(((orient_i - 1)*VoxNum+1) : orient_i*VoxNum);            
        tempin = reshape(tempin, N);            
        tempin = Params.Weight.*tempin;
        
        output = zeros(size(tempin));
        for ii = 1:numKernel
            output = output + ifftn(Params.D.*outSMV.SMV_kernel_conj(:,:,:,ii).*fftn(outSMV.SMV_mask(:,:,:,ii).*tempin));
        end

        y = y + output(:); 
      end
              
   elseif strcmp(transp_flag,'notransp') % y = A*x     
       y = zeros(OriNum*VoxNum, 1, 'single');
       x = single(x);                               % change to single format
                             
       x = reshape(x, N); 
        
       for orient_i = 1:OriNum
            Params = ParamsArray{orient_i}; 
            outSMV = outsmvArray{orient_i};
            DFx = Params.D.*fftn(x);
            tempin = zeros(size(x));
            for ii = 1:numKernel
                tempin = tempin + outSMV.SMV_mask(:,:,:,ii).*ifftn(outSMV.SMV_kernel(:,:,:,ii).*DFx);    
            end
            tempin = Params.Weight.*tempin;
            y(((orient_i - 1)*VoxNum+1) : orient_i*VoxNum) = tempin(:);           
       end
       
   end
end

%% To add frequency shift term 
% freqshift = chi_mo(end-OriNum+1:end);   % global frequency shift (critical for convergence test)
% chi_mo = chi_mo(1:end-OriNum);
% 
% chi_mo = real(reshape(chi_mo, N));

% function y = afun(x,transp_flag)
%    if strcmp(transp_flag,'transp')      % y = A'*x
%        y = zeros(VoxNum+OriNum, 1, 'single');
%        x = single(x);                                     % change to single format
%        for orient_i = 1:OriNum          
%             Params = ParamsArray{orient_i};
%             outSMV = outsmvArray{orient_i};
%             tempin = x(((orient_i - 1)*VoxNum+1) : orient_i*VoxNum);            
%             tempin = reshape(tempin, N);  
%             tempin = Params.Weight.*tempin;
%             y(VoxNum + orient_i) = sum(tempin(:));
%             output = zeros(size(tempin));
%             for ii = 1:numKernel
%                 output = output + ifftn(Params.D.*outSMV.SMV_kernel_conj(:,:,:,ii).*fftn(outSMV.SMV_mask(:,:,:,ii).*tempin)); 
%             end
% 
%             y(1:VoxNum) = y(1:VoxNum) + output(:);             
%        end
%               
%    elseif strcmp(transp_flag,'notransp') % y = A*x     
%        y = zeros(OriNum*VoxNum, 1, 'single');
%        x = single(x);                               % change to single format
%        fs = x(end-OriNum+1:end);                    % modeled frequency shift       
%        x = reshape(x(1:end-OriNum), N); 
%         
%         
%        for orient_i = 1:OriNum
%             Params = ParamsArray{orient_i}; 
%             outSMV = outsmvArray{orient_i};
%             DFx = Params.D.*fftn(x);
%             tempin = zeros(size(x));
%             for ii = 1:numKernel
%                 tempin = tempin + outSMV.SMV_mask(:,:,:,ii).*ifftn(outSMV.SMV_kernel(:,:,:,ii).*DFx); 
%             end
%             tempin = Params.Weight.*(tempin+ fs(orient_i));
%             y(((orient_i - 1)*VoxNum+1) : orient_i*VoxNum) = tempin(:);           
%        end


end
