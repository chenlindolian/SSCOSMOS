function mask_unrelyPhase = create_mask_unrelyPhase(GREPhase, GradMag_thresh)
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
%
% This function creates unreliable phase mask for QSM
% Reason:
% Due to phase wraps and phase noise in MRI, some phase measurements (and their laplacian) 
% are not reliable. Even if they can be unwrapped, they give errors that are not seperable
% from local phase, i.e. they are non-harmonic as compared to the harmonic
% background phase. Here we threshold the median-filtered gradient magnitude of
% raw phase to get the unreliable phase mask

% matlab 'medfilt3' function is only availabe in R2016b and after
useinternal = 0;        % use matlab internal precompiled medfilt3
versionstring = version;
if strfind(versionstring, 'R2017') % introduced in R2016b
    useinternal = 1;
    disp('User matlab internal mex version.')
end

if nargin < 2
    GradMag_thresh = pi/2;  % practically determined
end
Nkernel = [5, 5, 5];    % may need to change according to image resolution

TV = TVOP;

num_echoes = size(GREPhase, 4);
mask_unrelyPhase = zeros(size(GREPhase(:,:,:,1)));  % 3D mask

for echoii = 1:num_echoes
    phaseGrad = TV*GREPhase(:,:,:,echoii);    
    
    if useinternal
       phaseGradMag = medfilt3(sqrt(sum(phaseGrad.^2, 4)./3), 'replicate', Nkernel); % matlab internal one
    else
       phaseGradMag = mymedfilt3(sqrt(sum(phaseGrad.^2, 4)./3), Nkernel);
    end
    mask_unrelyPhase = mask_unrelyPhase | (phaseGradMag > GradMag_thresh);
end

