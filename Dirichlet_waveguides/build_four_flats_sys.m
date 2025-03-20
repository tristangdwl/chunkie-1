function sys = build_four_flats_sys(nover,ds,zk)


dl = ds(1);
dr = ds(2);

nch = 2.^(nover+2);
cparams = []; cparams.ta = 0; cparams.tb = pi;
chnkr_tmp = chunkerfuncuni(@flatinterface,nch,cparams);
chnkr_tmp = sort(chnkr_tmp);

chnkr1 = move(chnkr_tmp,[0;dl/2]);
chnkr2 = move(chnkr_tmp,[0;-dl/2]);
chnkrl = merge([chnkr1,chnkr2]);

% nch = 2.^(nover+2);
% cparams = []; cparams.ta = 0; cparams.tb = pi;
% chnkr_tmp = chunkerfuncuni(@flatinterface,nch,cparams);
% chnkr_tmp = sort(chnkr_tmp);

chnkr1 = move(chnkr_tmp,[0;dr/2]);
chnkr2 = move(chnkr_tmp,[0;-dr/2]);
chnkrr = merge([chnkr1,chnkr2]);

chnkr = merge([chnkrl,chnkrr]);

opts = []; opts.eps = 1e-8;
fkern_self = @(s,t) 0.25*1i*besselh(0,zk*sqrt((cos(s.r(1,:))-cos(t.r(1,:).')).^2+(s.r(2,:)-t.r(2,:).').^2));
sysmatll = chunkermat(chnkrl,fkern_self,opts);
sysmatrr = chunkermat(chnkrr,fkern_self,opts);

% fkern_lr = @(s,t) 0.25*1i*besselh(0,zk*sqrt((2+cos(s.r(1,:))-cos(t.r(1,:).')).^2+(s.r(2,:)-t.r(2,:).').^2));
% fkern_rl = @(s,t) 0.25*1i*besselh(0,zk*sqrt((-2+cos(s.r(1,:))-cos(t.r(1,:).')).^2+(s.r(2,:)-t.r(2,:).').^2));
% sysmat = chunkermat(chnkrr,fkern_self,opts);
wl = chnkrl.weights; wr = chnkrr.weights;


pts_l = [cos(chnkrl.r(1,:))-1; chnkrl.r(2,:)];
pts_r = [cos(chnkrr.r(1,:))+1; chnkrr.r(2,:)];
sysmatlr = chnk.helm2d.green(zk,pts_l,pts_r).*wl(:).';
sysmatrl = chnk.helm2d.green(zk,pts_r,pts_l).*wr(:).';

sysmat = [sysmatll,sysmatrl;sysmatlr,sysmatrr];



sys = [];
sys.ds = ds;
sys.chnkrl = chnkrl;
sys.chnkrr = chnkrr;
sys.chnkr = chnkr;
sys.zk = zk;

sys.pts = struct("r",[pts_l,pts_r]);


occ = 1024; rank_or_tol = opts.eps;
tic;
F = rskelf(sysmat,real(sys.pts.r),occ,rank_or_tol); 
tflam = toc
% 
sys.sysinv = @(rhs) rskelf_sv(F,rhs);
% sys.sysinv = @(rhs) sysmat\rhs;
sys.F = F;



end