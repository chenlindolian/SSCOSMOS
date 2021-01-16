function out = mycreate_SMVkernel(chi_mask, min_radius, max_radius, step_size_radius, N, voxel_size, phase_total)
% out = mycreate_SMVkernel(chi_mask, min_radius, max_radius,
% step_size_radius, N, voxel_size, phase_total)
% adapted from B.B.'s SSQSM code
% added switch for only creatnig the kernels without phase data
% X.L., 2017-08-16

if nargin < 7
    data_flag = 0;
else
    data_flag = 1;
end

% smv_radii from max to min, with min voxel_size
smv_radii = (max_radius:-step_size_radius:min_radius)*min(voxel_size);
disp(['SMV kernel size = ', num2str(smv_radii)])
num_kernel = length(smv_radii);

% meshgrid
[Y,X,Z] = meshgrid(-N(2)/2:(N(2)/2-1),-N(1)/2:(N(1)/2-1),-N(3)/2:(N(3)/2-1));

X = X * voxel_size(1);
Y = Y * voxel_size(2);
Z = Z * voxel_size(3);

unrely_tol = 1e-3;      % tol for unreliable boundary voxels 
SMV_kernel = zeros([N, num_kernel]);    % SMV kernel with different radius
mask_Sharp = zeros([N, num_kernel]);    % corresponidng mask and rims
mask_prev = zeros(N);

if data_flag == 1
    phase_Sharp = zeros(N);   % if there is phase data
end

count = 1;
for k = 1:num_kernel

    SMV = gen_SMVkernel_voxel_scaled( X, Y, Z, smv_radii(k));
    mask_rely = gen_SMVMask_new( SMV, chi_mask, unrely_tol);   
    
    if sum(mask_rely(:)) == 0
        continue
    end
    
    SMV_kernel(:,:,:,count) = SMV;
    mask_Sharp(:,:,:,count) = (mask_rely-mask_prev);
    mask_prev = mask_rely;
    
    % if there is phase data
    if data_flag == 1
        phase_Sharp = phase_Sharp +  mask_Sharp(:,:,:,count).* ifftn(SMV_kernel(:,:,:,count) .* fftn(phase_total));    
    end
    
    if count == 1
        Del_Sharp_Inv = SMV;
    end
    count = count + 1;
end
    
out.SMV_kernel = SMV_kernel;
out.SMV_mask = mask_Sharp;

if data_flag == 1
    out.mask_eval = phase_Sharp~=0;
    out.SMV_phase = phase_Sharp;
    out.SMV_inv_kernel = Del_Sharp_Inv; % For inverse Sharp
end

end

%% subfunctions
function SMV = gen_SMVkernel_voxel_scaled( X, Y, Z, smv_rad)
  
smv = (X.^2 + Y.^2 + Z.^2) <= smv_rad^2;    %nonnegtive,radially symmetric
smv = smv / sum(smv(:));                    % normalization

smv_kernel = zeros(size(X));
smv_kernel(1+end/2,1+end/2,1+end/2) = 1;    % delta
smv_kernel = smv_kernel - smv;              % s', delta - s

SMV = fftn(fftshift(smv_kernel));

end

function mask_rely = gen_SMVMask_new(SMV, chi_mask, unrely_tol)
  
mask_unrely = ifftn(SMV .* fftn(chi_mask));             % chi_mask is the support region (harmoic too)
mask_unrely = abs(mask_unrely) > unrely_tol;            % mask of unreliable phase estimates

mask_rely = chi_mask .* (chi_mask - mask_unrely);       % mask of reliable phase estimates
                                                        % seems to like detecting the boundary region (Laplician check)                                                       
end