function lap_opK = generate_LapKernel(N)
% generate numerical laplacian kernel of size N
ksize = [3, 3, 3];               
khsize = (ksize-1)/2;

kernelcore = [];
kernelcore(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
kernelcore(:,:,2) = [0 1 0; 1 -6 1; 0 1 0];
kernelcore(:,:,3) = [0 0 0; 0 1 0; 0 0 0];

Kernel = zeros(N);
Kernel( 1+N(1)/2 - khsize(1) : 1+N(1)/2 + khsize(1), 1+N(2)/2 - khsize(2) : 1+N(2)/2 + khsize(2), ...
    1+N(3)/2 - khsize(3) : 1+N(3)/2 + khsize(3) ) = kernelcore;

% laplacian kernel in k space
lap_opK = fftn(fftshift(Kernel));   % make real delta function, only do fftn

end

