function sys = build_blob_sys(nover,d,ctr,zk1,zk2)
% intervals at y=+- d/2, going from -L/2, L/2


cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 0;
cparams.maxchunklen = 4.0/zk2*2^(-nover); % setting a chunk length helps when the
                              % frequency is known
                              
pref = []; 
pref.k = 16;
start = tic; chnkr_temp = chunkerfunc(@(t) ellcfun(t,d,d),cparams,pref); 
t1 = toc(start);
chnkrs = [];
for i = 1:size(ctr,2)
    chnkri = move(chnkr_temp,[],ctr(:,i),[],[]);
    chnkrs = [chnkrs,chnkri];
end
chnkr = merge(chnkrs);
chnkr.n = normals(chnkr);



fprintf('%5.2e s : time to build geo\n',t1)

% plot(chnkr)
% daspect([1,1,1])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   .  .  .  derivative single layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fkern = @(s,t) chnk.helm2d.kern(zk1,s,t,'sprime',1)-chnk.helm2d.kern(zk2,s,t,'sprime',1);
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; sysmatds = chunkermat(chnkr,fkern,opts);
t1 = toc(start)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   .  .  .  single layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fkern = @(s,t) chnk.helm2d.kern(zk1,s,t,'s',1)-chnk.helm2d.kern(zk2,s,t,'s',1);
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; sysmats = chunkermat(chnkr,fkern,opts);
t1 = toc(start)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   .  .  .  double layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fkern = @(s,t) chnk.helm2d.kern(zk1,s,t,'d',1)-chnk.helm2d.kern(zk2,s,t,'d',1);
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; sysmatd = chunkermat(chnkr,fkern,opts);
t1 = toc(start)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   .  .  .  derivative double layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fkern = @(s,t) chnk.helm2d.kern(zk1,s,t,'dprime',1)-...
    chnk.helm2d.kern(zk2,s,t,'dprime',1);
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; sysmatdd = chunkermat(chnkr,fkern,opts);
t1 = toc(start)


sysmat = [sysmatd,sysmats;-sysmatdd,-sysmatds];
sysmat = sysmat+eye(size(sysmat,1));

% if strcmp(optscomp.type,'flam')
    % FLAM is slower to build but faster to solve
    occ = 1024;
    rank_or_tol = 1e-13;
    
    % xxflam = repmat(real(chnkr.r(1,:)),2,1);
    % xyflam = repmat(real(chnkr.r(2,:)),2,1);
    xxflam = repmat(real(chnkr.r(1,:)),1,2);
    xyflam = repmat(real(chnkr.r(2,:)),1,2);
    xflam = [xxflam(:).';xyflam(:).'];
    tic;
    F = rskelf(sysmat,xflam,occ,rank_or_tol); 
    tflam = toc
    sysinv = @(rhs) rskelf_sv(F,rhs);
    sys.F = F;
% else
%     tic;
%     [L, U] = lu(sysmat);
%     tlu = toc
%     sysinv = @(rhs) U\(L\rhs);
%     sys.L = L;
%     sys.U = U;
% end

sys.sysmat = sysmat;
sys.sysinv = sysinv;
sys.chnkr = chnkr;
sys.d = d;
sys.zk1 = zk1;
sys.zk2 = zk2;
sys.ctr = ctr;
sys.shape = 'ell';


end






function [r,d,d2] = ellcfun(t,a,b) 
% parameterization of an ellipse
r = [a*cos(t(:).');b*sin(t(:).')];
d = [-a*sin(t(:).');b*cos(t(:).')];
d2 = [-a*cos(t(:).');-b*sin(t(:).')];

end