function [Img_deblur] = Spiral_Deblur(Img,B0,Traj,dwell_time,r)
%SPIRAL_DEBLUR for spiral deblurring with k-space trajectory
% [idx] = kSpaceLocs(N,Profile_Order)
% Input:
%    -Img: images needs to be deblured
%    -B0: B0 map in Hz
%    -Traj: k-space trajectory (complex)
%    -dwell_time: acquisition dwelling time
%    -r: deblur radius
% Ouput:
%    -Img_deblur: deblurred images
% ------------------------------------------------------------------------
%     by Dan Zhu (dzhu6@uw.edu) Mar. 2025, Version 1.0
%     Please let me know if there's any bug to be fixed 
% ------------------------------------------------------------------------
%   Detailed explanation goes here
    Img_deblur=zeros(size(Img));Img=Img(:,:,:,:);
    [N1,N2,N3,N4]=size(Img);
    dwell_time=double(dwell_time);

    % Calculate t(k)
    tk=[0,dwell_time*(0:size(Traj,1)-1)*1e-6];
    radius=[0;abs(Traj(:,1))];
        
    % calculate rotation matrix
    kx=kSpaceLocs(2*r+1);
    rot_radius=min(abs(kx.'+1j*kx),max(radius));
    tk_rot=interp1(radius,tk,rot_radius);
    tk_rot=nan20(tk_rot);
    tk_rot(isnan(tk_rot))=0;
    rot_matrix=exp(-1j*2*pi*tk_rot);
    Img=zeropad(Img,[N1+r*2,N2+r*2,N3,N4],'size');
    
    % Low-pass preconditioning L1-L5 10.1002/mrm.29928
    lpfilt=lpfilt_sphere([2*r+1,2*r+1],0)*0.7+0.3; % L4 (min=0.3) best SSIM
   
    % deblur
    for sl=1:N3
        B0_S=B0(:,:,sl);
        rot_matrix1=rot_matrix.^permute(B0_S(:),[2,3,1]);
        rot_matrix1=imfft(rot_matrix1,[1;2])/(2*r+1).*lpfilt;
        rot_matrix1=reshape(rot_matrix1,[],N1*N2);
        for idx=1:N4
            Img_col=im2col(Img(:,:,sl,idx),[2*r+1 2*r+1],'sliding');
            Img_deblur(:,:,sl,idx)=reshape(sum(Img_col.*rot_matrix1,1),[N1,N2]);
        end
    end
end

function [idx] = kSpaceLocs(N,Profile_Order)
%KSPACE_LOCS Claculate the k-space locations covered by N pixels
% [idx] = kSpaceLocs(N,Profile_Order)
% Input:
%    -N: number of encoding steps (k-space matrix size in desired dim)
%    -Profile_Order: 'low-high' for centric order, 'linear' for linear
%                    order
% Ouput:
%    -idx: k-space locations nomalized from -0.5 to 0.5
% ------------------------------------------------------------------------
%     by Dan Zhu (dzhu6@uw.edu) Dec. 2020, Version 1.0
%     Please let me know if there's any bug to be fixed 
% ------------------------------------------------------------------------
    if nargin<2
        Profile_Order='linear';
    end
    idx=linspace(-0.5,0.5,floor(N/2)*2+1);idx=idx(1:N);
    switch lower(Profile_Order)
        case 'linear'
        case 'low-high'
            [~,od]=sort(abs(idx));
            idx=idx(od);
    end
end
