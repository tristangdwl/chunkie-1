function sys = build_dirichlet_sys(nover,d,L,zk)

chnkr = leaky.complex_interface(nover,d,L);

sys = [];
sys.d = d;
sys.zk = zk;
sys.chnkr = chnkr;

fkern = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');
opts = []; opts.eps = 1e-8;
sysmat = chunkermat(chnkr,fkern,opts)-eye(chnkr.npt)/2;
occ = 1024; rank_or_tol = opts.eps;
tic;
F = rskelf(sysmat,real(chnkr.r(:,:)),occ,rank_or_tol); 
tflam = toc
% 
sys.sysinv = @(rhs) rskelf_sv(F,rhs);
sys.F = F;
end