function regs = reg_fun(fun_name,opts)
%REG_FUN generates regularization functions
% regs = reg_fun(fun_name,opts)
% Input:
%    -fun_name: name of regularization norm
%    -opts: options
% Ouput:
%    -regs: regularization function handle
% ------------------------------------------------------------------------
%     by Dan Zhu (dzhu6@uw.edu) Feb 2024, Version 1.0
%     Please let me know if there's any bug to be fixed 
% ------------------------------------------------------------------------
    if ~iscell(fun_name),fun_name={fun_name};end;fun_name=fun_name(:);
    if ~iscell(opts),opts={opts};end;opts=opts(:);
    regs=repmat({},[size(fun_name,1),4]);
    for idx=1:size(fun_name,1)
        if strcmpi(fun_name{idx},'TV')
            regs{idx,1}=1;opt=opts{idx};dims=1;regs{idx,2}=1;
            if isfield(opt,'dims') && ~isempty(opt.dims),dims=opt.dims;end
            if isfield(opt,'weight') && ~isempty(opt.weight),regs{idx,2}=opt.weight;end
            dI=@(I,d)-diff(take(I,d,[1:size(I,d),size(I,d)]),[],d);
            dIp1=@(I,d)cat(d,take(I,d,1),diff(take(I,d,1:size(I,d)-1),[],d),-take(I,d,size(I,d)-1));
            regs{idx,3}=@(I)cell2mat(arrayfun(@(d)dI(I,d),permute(dims(:),[2:ndims(I)+1,1]),'UniformOutput',false));
            if length(dims)>1
                regs{idx,4}=@(I)sum(cell2mat(arrayfun(@(idx)dIp1(take(I,ndims(I),idx),dims(idx)),permute((1:length(dims))',[2:ndims(I),1]),'UniformOutput',false)),ndims(I));
            else
                regs{idx,4}=@(I)dIp1(I,dims);
            end
        end
    end
end

function Data_index = take(Data,dim,index)
    for idx=1:ndims(Data)
        inds{idx}=1:size(Data,idx);
    end
    if length(dim)>1
        for idx=1:length(dim)
            inds{dim(idx)}=index{idx};
        end
    else
        inds{dim}=index;
    end
    Data_index=Data(inds{:});
end