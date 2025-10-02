%% Fitting Example Code
% This code demonstrate the 3-parameter fitting model descirbed in our paper:
% Rapid 3D Whole-Brain High-Resolution T1 Quantification: Accelerating Standard Inversion Recovery with Stack-of-Spirals Turbo FLASH Acquisition. 
% Please cite the paper using DOI: 10.1002/mrm.70112.
% Please feel free to contact me with any questions: 
% - - - - - - - - - Zhen Hu (zhu22@jhmi.edu) - - - - - - - - - 

%% The 8 inversion time (TI) used in the paper: 
Nt = 7; 
t_begin = 100; t_end = 3000; %ms
TI = round(logspace(log(t_begin)/log(10), log(t_end)/log(10), Nt)); % 7 TI are logarithmically spaced from 100 ms to 3000 ms
TI = [TI,6000];  % A following 8th TI of 6000 ms
Nt = size(TI,2);

%% Load 8 TI images:
% The example fitting dataset consists of eight inversion time (TI) images acquired using our 3D SOS-TFL readout with a 2000 ms saturation time (Tsat)
load("Example_fitting_data.mat"); data = ImgsTI;

%% Fitting model and algorithm
% A three-parameter model was used to fit T1, equilibrium magnetization M0, and inversion degree α. 
% M_z(TI)=M_0∙(1-(1-(1-e^(-T_sat/T_1 ))∙cos⁡(α)) ∙e^(-TI/T_1 ) )    

T_sat = 2000; % ms
fit_model = @(x,TI)(x(1).*(1-(1-(1-exp(-T_sat./x(2))).*cos(deg2rad(x(3)))).*exp(-TI./x(2)))); % x(1)/x(2)/x(3):M0/T1/α
fit_model_abs = @(x,TI)abs(fit_model(x, TI));
%% Fitting calculation
Maps = three_params_fitting(data, TI, fit_model_abs); 
imshow(Maps(:,:,61,2),[0,2500]);colormap jet; colorbar; 
%The central slice of the 3D T1 map, corresponding to the fourth column in Figure 3.

%%
function Fitted_maps = three_params_fitting(Img, TI, fit_model)
% Pixel-wise 3-parameter fitting function.
% Input: 
%       Img: A 4D matrix with all TI images concatenated in the fourth
%           dimension. Dimension: Nx*Ny*N_slice*N_TI.
%       TI: Inversion time array. Dimension: 1*N_TI
%       fit_model: 3-parameter fitting model
% Output: 
%       Fitted_Maps: a 4D matrix containing the proton density map, T1
%          map and the inversion degree map. Dimension: Nx*Ny*N_slice*3.

    % Matlab nonlinear curve-fitting problems solver using Levenberg–Marquardt algorithm. 
    options = optimset('lsqcurvefit'); 
    options = optimset(options,'TolX',1e-16,'TolFun',1e-16,'MaxIter',1000,'MaxFunEvals',1000);  
    options = optimset(options,'Display','off');
    options = optimset(options,'Algorithm','levenberg-marquardt');

    N_slice = size(Img,3);
    Fitted_maps = zeros([size(Img,1:3),3]);
    mask = double(abs(squeeze(Img(:,:,:,end)))>= abs(mean(squeeze(Img(:,:,:,end)),'all'))*0.1); % using the last TI image to draw a human brain mask.
    
    for x_pixel = 1:size(Img,1)
        for y_pixel = 1:size(Img,2) % Note: if you can use parfor in your computer, then use parfor to accelerate.
            for i_slice = 1:N_slice
    
                if mask(x_pixel, y_pixel, i_slice) == 1
    
                    xdata = TI;
                    ydata = squeeze(Img(x_pixel,y_pixel,i_slice,:))';

                    x0 = [abs(max(ydata)), 1000, 180]; % initial estimate of [M0, T1 (ms), nominated FA (degree)].
                    lb = [abs(max(ydata))/2, 100, 90];  % lower bond
                    ub = [abs(max(ydata))*2, 6000, 180]; % upper bond
                    [x_fit,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(fit_model, x0, xdata, ydata,lb,ub,options); %resnorm: squared 2-norm of the residual at x: sum((fun(x,xdata)-ydata).^2).
                    Fitted_maps(x_pixel, y_pixel, i_slice,:) = x_fit;
                     
                end % end if 
    
            end % i_slice
   
        end % y_pixel
        fprintf("---------- x_pixel = %d finished ---------.\n", x_pixel);
    end % x_pixel
end