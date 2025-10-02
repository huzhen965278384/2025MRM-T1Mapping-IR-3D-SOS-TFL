%% Example Reconstruction code
% This code uses the k-space data from a single TI image to demonstrate the reconstruction process 
% described in our paper: Rapid 3D Whole-Brain High-Resolution T1 Quantification: 
% Accelerating Standard Inversion Recovery with Stack-of-Spirals Turbo FLASH Acquisition. 
% Please cite our paper using DOI: 10.1002/mrm.70112.
% Please contact us if you have any questions:
% - - - - - - - - - Dan Zhu (dzhu6@uw.edu), Zhen Hu (zhu22@jhmi.edu)  - - - - - - - - - 

clear all; close all; clc;
%% Preparation
addpath(genpath('ReconFunctions'));
load("Example_recon_data.mat")
% Variables inside the Example_recon_data.mat:
% 1. k-space data (k_raw)
% 2. Trajectory (Traj), Density weighting vector (wi)
% 3. Sensitivity map (S), inverse sensitivity map (Sinv)
% 4. B0 map (B0Map)
% 5. Image Size (Img_Size), image recon size (Img_Size_recon)
% 6. Dwell time (Dwell_time) 

%% Through-plane SENSE reconstruction
fprintf('Start the through-plane SENSE reconstruction.\n')
OS_factor = 1.25; % in-plane oversampling factor
Img_Size_OS = round(Img_Size(1:2).*OS_factor); % in-plane oversampling for spiral

FT = NUFFT(Traj,wi,[0,0],Img_Size_OS); % non-uniform Fourier transform
b = zeros(size(k_raw).*[1,1,2,1]); % wFSI = wk = b
b(:,:,1:2:end,:) = k_raw.*sqrt(wi); % zero-padding for through-plane SENSE
mask = false([1,1,size(b,3)]);mask(:,:,1:2:end) = true; % SENSE undersampling pattern

%% In-plane Compressed Sensing reconstruction
fprintf('Start the in-plane Compressed Sensing Reconstruction.\n')
% Define the reconstruction parameters: 
par_recon.CSIter = 15;  % Iteration number = 15
par_recon.TV_weight = 0.003; % Weight for the compressed sensing regularization. 
par_recon.Tol = 0; % Stop criteria 

imfft = @(I)FT*fftshift(fft(fftshift(I,3),[],3),3);
imifft = @(k)ifftshift(ifft(ifftshift(FT'*k,3),[],3),3);
A = @(I)imfft(I.*S).*mask; 
At = @(k)sum(imifft(k.*mask).*Sinv,4); % A*I = k; At*A*I = At*k; 
I_recon = At(b);

weight_scales = max(abs(I_recon(:)));
regs = reg_fun('TV',struct('dims',1:3,'weight',par_recon.TV_weight*weight_scales)); % Use total variation as sparsity transform.
I_recon = CG_Nonlinear_wNorms(A,At,b,regs,I_recon,par_recon.CSIter,par_recon.Tol*weight_scales); % Nonlinear conjugate gradient
I_recon = CG_Nonlinear_wNorms(A,At,b,regs,I_recon,par_recon.CSIter,par_recon.Tol*weight_scales);
I_recon = CG_Nonlinear_wNorms(A,At,b,regs,I_recon,par_recon.CSIter,par_recon.Tol*weight_scales);

%% Spiral deblurring
par_recon.deblur_radius = 8; 

I_recon = zeropad(I_recon,[Img_Size_OS,Img_Size(3)],'size'); 
I_recon = Spiral_Deblur(I_recon,real(B0Map),Traj, Dwell_time, par_recon.deblur_radius); % deblur
I_recon = zeropad(I_recon,Img_Size,'size');
I_recon = imresize3(I_recon,Img_Size_recon);

imshow(abs(I_recon(:,:,61)), [0,18]); colormap gray; colorbar; % show the center slice of the 3D recon image
save("Example_recon_fit_result.mat","I_recon");
