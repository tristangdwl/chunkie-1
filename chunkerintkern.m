function fints = chunkerintkern(chnkr,kern,dens,targs,opts)
%CHUNKERINTKERN compute the convolution of the integral kernel with
% the density defined on the chunk geometry. 
%
% input:
%   chnkr - chunks description of curve
%   kern - integral kernel taking inputs kern(s,t,sn,tn) where 
%          s is a source, sn is the normal at the source, t is a target
%          and tn is the normal at the target
%   opdims - input and output dimensions of the kernel (opdims(1) dimension
%           output, opdims(2) dimension of input)
%   dens - density on boundary, should have size opdims(2) x k x nch
%          where k = chnkr.k, nch = chnkr.nch
%   targs - targ(1:2,i) gives the coords of the ith target
%   opts - structure for setting various parameters
%       opts.targstau - if provided, the normals at the targets
%       opts.usesmooth - if = 1, then just use the smooth integration
%          rule for each chunk. if = 0, adaptive integration 
%          (quadgk) is used. if = 2, a hybrid method is used 
%          where the smooth rule is used for points where the 
%          smooth rule is accurate to opts.gausseps digits for 
%          gauss's id and uses adaptive at other points (default 0)
%       opts.quadgkparams - if non-empty this is a cell structure
%       containing string,value pairs to be sent to quadgk (default {})
%       opts.gausseps - if the hybrid method is used, the smooth 
%       rule is applied at points for which Gauss' ID is accurate 
%       with the smooth rule to absolute error opts.gausseps (default 1e-8)
%
% output:
%   fints - opdims(1) x nt array of integral values
%
% see also QUADGK

% determine operator dimensions using first two points

rs = chnkr.r(:,1:2);
ds = chnkr.d(:,1:2); dsn = sqrt(sum(ds.^2,1)); 
ds = bsxfun(@rdivide,ds,dsn);

ftemp = kern(rs(:,1),rs(:,2),ds(:,1),ds(:,2));
opdims = size(ftemp);

if nargin < 6
    opts = [];
end

[dim,nt] = size(targs);
assert(dim==2,'only dimension two tested');

if ~isfield(opts,'usesmooth'); opts.usesmooth = false; end
if ~isfield(opts,'quadgkparams'); opts.quadgkparams = {}; end
if ~isfield(opts,'gausseps'); opts.gausseps = 1e-8; end
if ~isfield(opts,'targstau'); opts.targstau = zeros(dim,nt); end
if ~isfield(opts,'verb'); opts.verb= false; end

targstau = opts.targstau;

[dim,nt2] = size(targstau);

assert(dim==2 && nt2==nt,...
    'opts.targstau should have same dimensions as targs');

if opts.usesmooth == 1
    fints = chunkerintkern_smooth(chnkr,kern,opdims,dens, ...
        targs,targstau,opts);
elseif opts.usesmooth == 0
    fints = chunkerintkern_adap(chnkr,kern,opdims,dens, ...
        targs,targstau,opts);
elseif opts.usesmooth == 2
    fints = zeros(opdims(1),nt);
    optssw = []; optssw.gausseps = opts.gausseps; 
    optssw.justsmoothworks = true;
    sw = chunkerin(chnkr,targs,optssw);
    fints(:,sw) = reshape(...
        chunkerintkern_smooth(chnkr,kern,opdims,dens, ...
        targs(:,sw),targstau(:,sw),opts),opdims(1),nnz(sw));
    fints(:,~sw) = reshape(...
        chunkerintkern_adap(chnkr,kern,opdims,dens, ...
        targs(:,~sw),targstau(:,~sw),opts),opdims(1),nnz(~sw));
    fints = fints(:);
end

end



function fints = chunkerintkern_smooth(chnkr,kern,opdims,dens, ...
    targs,targstau,opts)

k = chnkr.k;
nch = chnkr.nch;

assert(numel(dens) == opdims(2)*k*nch,'dens not of appropriate size')
dens = reshape(dens,opdims(2),k,nch);

[~,w] = lege.exps(k);
[~,nt] = size(targs);

tau = taus(chnkr);

fints = zeros(opdims(1)*nt,1);

% assume smooth weights are good enough
for i = 1:nch
    densvals = dens(:,:,i); densvals = densvals(:);
    dsdtdt = sqrt(sum(abs(chnkr.d(:,:,i)).^2,1));
    dsdtdt = dsdtdt(:).*w(:)*chnkr.h(i);
    dsdtdt = repmat( (dsdtdt(:)).',opdims(2),1);
    densvals = densvals.*(dsdtdt(:));
    kernmat = kern(chnkr.r(:,:,i),targs, ...
        tau(:,:,i),targstau);

    fints = fints + kernmat*densvals;
end

end

function fints = chunkerintkern_adap(chnkr,kern,opdims,dens, ...
    targs,targstau,opts)

k = chnkr.k;
nch = chnkr.nch;

assert(numel(dens) == opdims(2)*k*nch,'dens not of appropriate size')
dens = reshape(dens,opdims(2),k,nch);

[~,~,u] = lege.exps(k);
[~,nt] = size(targs);

fints = zeros(opdims(1)*nt,1);

% using adaptive quadrature
[rc,dc] = exps(chnkr);
for i = 1:nch
    if opts.verb; fprintf('chunk %d integral\n',i); end
    rci = rc(:,:,i);
    dci = dc(:,:,i);
    densvals = dens(:,:,i); densvals = densvals.';
    densc = u*densvals; % each column is set of coefficients
                    % for one dimension of density on chunk
    for j = 1:nt
        indj = (j-1)*opdims(1);
        for l = 1:opdims(1)
            ind = indj+l;
            temp = chunkerintchunk_kernfcoefs(kern,opdims,l,...
                densc,rci,dci,targs(:,j),targstau(:,j));

            fints(ind) = fints(ind) + temp*chnkr.h(i);
        end
    end
end

end
