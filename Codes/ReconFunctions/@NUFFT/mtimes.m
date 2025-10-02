function ress = mtimes(a,bb)
% performs the normal nufft

% Dan modification for multi dimensions & faster
N=size(bb);
if a.adjoint; ress=zeros([a.imSize,N(3:end)]);
else;  ress=zeros([a.dataSize,N(3:end)]); end
% Dan modification for multi dimensions & faster

for m=1:prod(N(3:end))
    if a.adjoint
        ress(:,:,m) = reshape(nufft_adj(reshape(bb(:,:,m),N(1)*N(2),1).*a.w(:),a.st)/sqrt(prod(a.imSize)),a.imSize);
    else
        ress(:,:,m)  = reshape(nufft(bb(:,:,m),a.st)/sqrt(prod(a.imSize)).*a.w(:),a.dataSize);
    end
end
