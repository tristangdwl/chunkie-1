zk = 1;
ctr = [0,30,11;0,4,-4];
% ctr = ctr(:,2);
d = 1;
sys = build_blob_sys(2,d,ctr,zk,2*zk);
L_comp= 40;

tol = 1e-10;
sys.comp_data = comp_far(sys,L_comp,tol);
%%
idp = 0;
igrad = 0;
src = []; src.r = [-4;0]; src.n = [1;0]; 
chnkr = vert_chnkr(1,L_comp);
chnkr = move(chnkr,[],[-4;0]);
targ = chnkr;
% targ = []; targ.r = [-2;0]; targ.n = [1;0]; 

tic;
val = comp_obs(src, targ, sys, 's');
% val=compress.apply_comp_far(sys.comp_data.f2f,sys,src,targ,idp,igrad);
toc;
tic;
val_true = wave_kern.s(src,targ,sys);
toc;
[val,val_true,val-val_true];
max(abs(val-val_true),[],"all")
