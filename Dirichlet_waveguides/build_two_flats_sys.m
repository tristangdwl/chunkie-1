function sys = build_two_flats_sys(nover,d,zk)
% intervals at y=+- d/2, going from -L/2, L/2


nch = 2.^(nover+2);
cparams = []; cparams.ta = 0; cparams.tb = pi;
chnkr_tmp = chunkerfuncuni(@flatinterface,nch,cparams);
chnkr_tmp = sort(chnkr_tmp);

chnkr1 = move(chnkr_tmp,[0; d/2]);
chnkr2 = move(chnkr_tmp,[0;-d/2]);
chnkr = merge([chnkr1,chnkr2]);

opts = []; opts.eps = 1e-8;
fkern_self = @(s,t) 0.25*1i*besselh(0,zk*sqrt((cos(s.r(1,:))-cos(t.r(1,:).')).^2+(s.r(2,:)-t.r(2,:).').^2));
sysmat = chunkermat(chnkr,fkern_self,opts);

sys = [];
sys.d = d;
sys.chnkr = chnkr;
sys.zk = zk;

occ = 1024; rank_or_tol = opts.eps;
tic;
F = rskelf(sysmat,real(chnkr.r(:,:)),occ,rank_or_tol); 
tflam = toc
% 
sys.sysinv = @(rhs) rskelf_sv(F,rhs);
sys.F = F;
end