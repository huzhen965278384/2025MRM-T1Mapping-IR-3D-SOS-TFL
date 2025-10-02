function [x,obj] = CG_Nonlinear_wNorms(A,At,b,regs,x,nIter,TOL)
%CG_NONLINEAR_WNorms Solve min(0.5*(Ax-b)^2+weight*lp(x)^p) problem with
%Conjugate Gradient. l1(x) is the l1 norm (such as TV, wavelet) and l2(x)
%is the l2 norm (such as Tikhonov). There can be multiple lp(x) terms.
% Input
%   -A: Ecoding Matrix in function handle
%   -At: Transpose of Ecoding Matrix in function handle
%   -b: Acquisition
%   -regs: a cell with {p,weight,reg,regT) p: norm; weight: weight of this 
%        norm; reg:regularization in function handle; regT: transpose of reg
%        could be in multiple columns for multip.
%   -x: Initial (REQUIRED)
%   -nIter: Number of Iterations (Optional, **Default: size(x0))
%   -Tol: Tolerance (Optional, Default: norm(b)*1e-3)
% Output:
%   -x: Fitted results
%   -obj: object function per iteration
%
% ** Conjugate Gradient is thoretically proven to be converged in R
% iterarions, for R is the rank of A.
%
% Object Function: f(x)=0.5||Ax-b||^2
% Gradient: G=AT(Ax-b)
%   ----------------------------------------------------------------------
%     by Dan Zhu (dzhu6@uw.edu) Feb 2024, Version 1.1
%        add wavelet denoising
%     by Dan Zhu (dzhu6@uw.edu) March 2019, Version 1.0
%     Please let me know if there's any bug to be fixed 
%   ----------------------------------------------------------------------
% input check
    if ~exist('x','var') || isempty(x),x=At(b);end
    if ~exist('nIter','var') || isempty(nIter),nIter=300;end
    if ~exist('TOL','var') || isempty(TOL),TOL=1e-5*norm(b(:))^2;end
% function handles
    % regx=reg(x), use regx indstead of reg(x) to avoid recomputing (more efficient)
    lp=@(regx,p)sum((abs(regx).^2+eps).^(p/2),'all');                                       % lp-norm
    glp=@(regx,regT,p)regT(p*regx.*(abs(regx).^2+eps).^(p/2-1));                            % gradient of lp-norm
    regx_fun=@(x)arrayfun(@(idx)regs{idx,3}(x),1:size(regs,1),'UniformOutput',false);       % reg(x) for all
    f=@(regx)sum(arrayfun(@(idx)regs{idx,2}*lp(regx{idx},regs{idx,1}),1:size(regs,1)));     % object function 4 all
    g=@(regx)sum(cell2mat(permute(arrayfun(@(idx)regs{idx,2}*glp(regx{idx},regs{idx,4},regs{idx,1}),...
        (1:size(regs,1))','UniformOutput',false),[2:ndims(regx{1})+1,1])),ndims(regx{1})+1);% gradient function 4 all
% Line search parameters
    ls_alpha=1e-2;    % tolerance weight minimum obj descent
    ls_t0=1;          % initial line search starting step size
    ls_beta=0.618;    % step size shrink ratio(Golden Ratio)
    maxlsiter=150;    % max line search iteration
% Initialization
    regx=regx_fun(x);res=A(x)-b;
    obj(1)=0.5*sum(abs(res).^2,'all')+f(regx);
    g0=At(res)+g(regx);
    d=-g0; % d: descent direction
% Iteration
    fprintf('iter:');
    for iter=1:nIter
        fprintf('%d\t',iter);
        % line search iteration
        t=ls_t0;                                                % initialize step size
        regx=regx_fun(x+t*d);                                   % initialize new regx
        Ad=A(d);                                                % for efficiency: A(x+td)-b=Ax-b+tAd=res+tAd
        obj(iter+1)=0.5*sum(abs(res+t*Ad).^2,'all')+f(regx);    % initialize object function
        lsiter=0;                                               % initialize line search iteration
        min_descent=ls_alpha*abs(g0(:)'*d(:));                  % alpha*t*abs(gTd)
        while obj(iter+1)>obj(iter)-t*min_descent && lsiter<maxlsiter
            t=t*ls_beta;                                        % update step size
            regx=regx_fun(x+t*d);                               % update new regx
            obj(iter+1)=0.5*sum(abs(res+t*Ad).^2,'all')+f(regx);% update object function
            lsiter=lsiter+1;                                    % update line search iteration
        end                                                     % end of line search iteration
        if lsiter == maxlsiter                  % failing line search
            error('Aborted for max line search iteration Reached. Might have a bug in operators.');
        end
        if lsiter>2,ls_t0=ls_t0*ls_beta; end    % starting step is large, decrease 
        if lsiter<1,ls_t0=ls_t0/ls_beta; end    % starting step is small, increase

        x=x+t*d;                                % update x
        res=res+t*Ad;                           % update res
        normG=norm(g0(:));                      % old l2 norm of g0
        g0=At(res)+g(regx);                     % update g
        beta=(norm(g0(:))/normG)^2;             % descent direction update weight=norm(new_g)/norm(old_g)
        d=-g0+beta*d;                           % update descent direction
        if abs(obj(iter)-obj(iter+1))<TOL,break;end  % Relative Error=f(x_(i-1))-f(x_i)
    end
fprintf('\n');
end