function [ Ipad ] = zeropad( I , pad_size, pad_type, pad_value)
%ZEROPADDING Pad image with zeros
%   ZEROPAD(I) increase image size by 2 for each dimension by adding a
%   border outsize the original image with intensities of zero. I could be
%   a vector or a multi-dimension matrix.
%   ZEROPAD(I,N) or ZEROPAD(I,N,'factor') increase image size by N 
%   for each dimension if N is a positive scaler. If N<1, it will truncate 
%   the image with a factor of N. If N is a vector, the size of each 
%   dimension is processed seperately. The lenth of N should meet with the 
%   dimension of image if N is not a scalar.
%   ZEROPAD(I,pad_size,'size') increase image size to pad_size. If
%   pad_size is a scaler, it will pad all dimension to pad_size. Otherwise
%   the length of pad_size should meet with the dimension of image.
%   ZEROPAD By Dan Zhu (dzhu6@uw.edu) 2018

% Checking validity parameters
    % if no input, give an error
    if nargin==0
        error('No input images');
    end
    
    % if no pad_size, use twice
    if nargin==1
        pad_size=2;
        pad_type='factor';
        warning(['No pad size specified. Pad the image to twice its ', ...
            'original size by default']);
    end
    
    % if no pad_size, use twice
    if nargin==2
        pad_type='factor';
        warning(['No pad type specified. New image size is calculated by ', ...
            'multipicating the giving factor by default']);
    end
    
    % if pad_type not recognized, use twice
    if ~ischar(pad_type) ||(~strcmp(pad_type,'factor') && ~strcmp(pad_type,'size'))
        pad_type='factor';
        warning(['Unrecognized pad type. New image size is calculated by ', ...
            'multipicating the giving factor by default']);
    end
    
    % if pad_size is a cell get only first item
    if iscell(pad_size)
        if ~isempty(pad_size)
            pad_size=pad_size{1};
        else
            pad_size=2;
            pad_type='factor';
        end
        warning(['The parameter for dimensions should not be a cell. '...
            'Using only the first element.']);
    end
    
    % pad_size should be numeric
    if ~isnumeric(pad_size)  % pad_size not numeric
        pad_size=2;
        pad_type='factor';
        warning(['Unrecognized pad size. Pad the image to twice its ', ...
            'original size by default']);
    end
    
    % pad_size should be a scalar or vector
    if ~isvector(pad_size)  % pad_size not a vector
        pad_size=pad_size(:).';
        warning('Pad size vectorized.');
    end
    
    % if pad_size is a scalar, transfer to vector
    if isscalar(pad_size)
        dim_Image=true_ndims(I);
        if dim_Image>1
            pad_size=repmat(pad_size,[1,dim_Image]);
        end
    end
    
    % if no pad_value input, pad_value=0;
    if ~exist('pad_value')
        pad_value=0;
    end
    
    % if pad_size is a vector, it's length should meet the dimension of I
    dim_Image=true_ndims(I);
    if length(pad_size)>dim_Image
        if prod(pad_size(dim_Image+1:end))>1
            warning('pad_size is too long, truncate pad_size to image dimension.');
        end
        pad_size=pad_size(1:dim_Image);
    elseif ~isscalar(pad_size) && length(pad_size)<dim_Image
        pad_size=pad_size(:).';
        if strcmp(pad_type,'size')
            pad_size1=size(I);
        else
            pad_size1=ones(1,dim_Image);
        end
        pad_size1(1:length(pad_size))=pad_size;
        pad_size=pad_size1;
        warning('pad_size is too short, ignore truncation in other dims.');
    end
    

    switch pad_type
        case 'size' % if pad_type is 'size' pad_size must be positive integer
            if sum(round(pad_size)~=pad_size)
                pad_size=fix(pad_size);
                warning(['pad_size must be integer! Fractional part of all ',...
                  'decimal numbers are ignored']);
            end
            if sum(pad_size<=0)
                pad_size(pad_size<=0)=1;
                warning(['pad_size must be positive! All non-positive terms ',...
                  'are changed to 1']);
            end
        case 'factor'
            pad_size=pad_size.*size(I,1:dim_Image);
            pad_type='size';
            if sum(round(pad_size)~=pad_size)
                pad_size=round(pad_size);
            end
            if sum(pad_size<=0)
                pad_size(pad_size<=0)=1;
                warning(['pad_size must be positive! All non-positive terms ',...
                  'are changed to 1']);
            end
    end
    
% Calculate the index for each dimension
    idx_ori={};idx_pad={};
    for dim=1:true_ndims(I)
        l=size(I,dim);  % length in this dimension, used for scaling
        l_pad=pad_size(dim);
        
        l_center=round(l/2); % one more pixel before center if odd
        l_center_pad=round(l_pad/2);
        
        if l<=l_pad
            idx_ori{dim}=1:l;
            l_left=l_center-1;
            l_right=l-l_center;
            idx_pad{dim}=(l_center_pad-l_left):(l_center_pad+l_right);
        else
            l_left=l_center_pad-1;
            l_right=l_pad-l_center_pad;
            idx_ori{dim}=(l_center-l_left):(l_center+l_right);
            idx_pad{dim}=1:l_pad;
        end  
    end
    
% Zeropadding
    if isscalar(pad_size)
        pad_size=[pad_size,1];
    end
    Ipad=ones(pad_size(:).')*pad_value;
    Ipad(idx_pad{:})=I(idx_ori{:});
end

