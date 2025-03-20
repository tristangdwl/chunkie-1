addpath(genpath('~/FastAlgs_backup/chunkie-1'))
addpaths_loc();
load("two_tube_data.mat")

zk = sys_L.zk;
zk2 = 4*zk;


d_R = 0.5;
sys = build_blob_sys(2,d_R,[0;0],zk,2*zk);
% ctr_x = [6,8,10,20];
% ctrs = [ctr_x;0*ctr_x];
ctr = [10;0];

L_skel= 40;
tol = 1e-10;

Lplot = 26;
nplt = 160;

xts = linspace(-Lplot,Lplot,nplt);
yts = linspace(-Lplot,Lplot,nplt);
[X,Y] = meshgrid(xts,yts);
xtargs = X(:).'; ytargs = Y(:).';

targ = [];
targ.r = [xtargs;ytargs];
targ.n = [1;0]+0*targ.r;

il = (targ.r(1,:)<chnkr_L.r(1,1));
ir = (targ.r(1,:)>chnkr_R.r(1,1));
im = ~il & ~ir;

targl = [];
targl.r = targ.r(:,il);
targl.n = targ.n(:,il);

targm = [];
targm.r = targ.r(:,im);
targm.n = targ.n(:,im);

targr = [];
targr.r = targ.r(:,ir);
targr.n = targ.n(:,ir);

ntarg = size(targ.r,2);


if true

istart = tic;
sys_R = build_blob_sys(2,d_R,ctr,zk,zk2);
sys_R.comp_data = comp_far(sys_R,L_skel,tol);

wts_L = chnkr_L.weights; wts_L = wts_L(:);
wts_R = chnkr_R.weights; wts_R = wts_R(:);
wts = chnkr_L.weights; wts = wts(:);

fkern_SR = @(s,t) - comp_obs(s,t,sys_R,'s');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
sysmats_R_obs = chunkermat(chnkr_R,fkern_SR,opts);
% sysmats_R = fkern_SR(chnkr_R,chnkr_R).*wts_R(:).';
t2 = toc(start)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   .  .  .  derivative double layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fkern_Dp_R = @(s,t) - comp_obs(s,t,sys_R,'dp');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
sysmatdd_R_obs = chunkermat(chnkr_R,fkern_Dp_R,opts);
% sysmatdd_R = fkern_Dp_R(chnkr_R,chnkr_R).*wts_R(:).';
t2 = toc(start)
% %%
fkern_D_R = @(s,t) - comp_obs(s,t,sys_R,'d');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
sysmatd_R_obs = chunkermat(chnkr_R,fkern_D_R,opts);
% sysmatd_R = fkern_Dp_R(chnkr_R,chnkr_R).*wts_R(:).';
t2 = toc(start)


fkern_Sp_R = @(s,t) - comp_obs(s,t,sys_R,'sp');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
sysmatsd_R_obs = chunkermat(chnkr_R,fkern_Sp_R,opts);
% sysmatsd_R = fkern_Sp_R(chnkr_R,chnkr_R).*wts_R(:).';
t2 = toc(start)


sysmat_R_obs = [sysmatd_R_obs,sysmats_R_obs;-sysmatdd_R_obs,-sysmatsd_R_obs];
%%




sysmat = [sysmat_L,sysmat_LR;sysmat_RL,sysmat_R+sysmat_R_obs];
sysmat = sysmat-eye(size(sysmat,1));


%%
src = [];
src.r = [-31;0];
src.n = [1;0];
tic;
rhs1 = [d_guide_kern(src,chnkr_L,sys_L,'s');-d_guide_kern(src,chnkr_L,sys_L,'sp')];
rhs2 = zeros(2*chnkr_R.npt,1);
rhs = [rhs1;rhs2];
toc
tic;
soln = sysmat\rhs;
ddens = soln([1:chnkr_L.npt,2*chnkr_L.npt+(1:chnkr_R.npt)]);
sdens = soln([chnkr_L.npt+(1:chnkr_R.npt),2*chnkr_L.npt+chnkr_R.npt+(1:chnkr_R.npt)]);
toc

%%
uin = zeros(ntarg,1);
uscat = zeros(ntarg,1);

uin(il) = d_guide_kern(src,targl,sys_L,'s');
%%

w = chnkr.weights; w = w(:);
ssrc = [];
ssrc.r = chnkr.r(:,:);
ssrc.n = chnkr.n(:,:);
ssrc.charge = sdens.*w;

dsrc = [];
dsrc.r = chnkr.r(:,:);
dsrc.n = chnkr.n(:,:);
dsrc.charge = ddens.*w;

id_L = 1:chnkr_L.npt;
ssrc_L = [];
ssrc_L.r = ssrc.r(:,id_L);
ssrc_L.n = ssrc.n(:,id_L);
ssrc_L.charge = ssrc.charge(id_L);

dsrc_L = [];
dsrc_L.r = dsrc.r(:,id_L);
dsrc_L.n = dsrc.n(:,id_L);
dsrc_L.charge = dsrc.charge(id_L);

id_R = chnkr_L.npt+(1:chnkr_R.npt);
ssrc_R = [];
ssrc_R.r = ssrc.r(:,id_R);
ssrc_R.n = ssrc.n(:,id_R);
ssrc_R.chrg = ssrc.charge(id_R);

dsrc_R = [];
dsrc_R.r = dsrc.r(:,id_R);
dsrc_R.n = dsrc.n(:,id_R);
dsrc_R.chrg = dsrc.charge(id_R);

tic
uscat(il) = d_guide_kern_lump(ssrc_L,dsrc_L,targl,sys_L);
uscat(im) = d_two_flat_kern_lump(ssrc,dsrc,targm,sys_M);
uscat(ir) = wave_kern.lump_source(targr,ssrc_R,dsrc_R,sys_R);
toc

%%

figure(5); clf
% subplot(1,2,1)
tools.pplot(X,Y,reshape(imag(uin-uscat),nplt,nplt),'$\Im(u_{\rm{tot}})$')
chnkr_L_p = real(sys_L.chnkr.r(:,:)); chnkr_L_p(:,chnkr_L_p(1,:)>chnkr_L.r(1,1)) = NaN;
% chnkr_M_p = real(sys_M.pts.r); chnkr_M_p(:,abs(chnkr_M_p(1,:))>1) = NaN;
chnkr_M_p = real([cos(sys_M.chnkr.r(1,:));sys_M.chnkr.r(2,:)]); chnkr_M_p(:,chnkr_M_p(1,:)<0) = NaN;
chnkr_R_p = real(sys_R.chnkr.r(:,:)); chnkr_R_p(:,chnkr_R_p(1,:)<chnkr_R.r(1,1)) = NaN;
    hold on
    plot(real(chnkr_L_p(1,:)),real(chnkr_L_p(2,:)),'k.','markersize',1)
    plot(real(chnkr_M_p(1,:)),real(chnkr_M_p(2,:)),'k.','markersize',1)
    plot(real(chnkr_R_p(1,:)),real(chnkr_R_p(2,:)),'k.','markersize',1)
    plot(real(chnkr.r(1,:)),real(chnkr.r(2,:)),'r.','markersize',1e-3)
    hold off

ti = toc(istart)

uin_blob = uin;
uscat_blob = uscat;
utot_blob = uin - uscat;

end


if true


istart = tic;
sys_R = build_empty_sys(zk);
sys_R.comp_data = comp_far(sys_R,L_skel,tol);

wts_L = chnkr_L.weights; wts_L = wts_L(:);
wts_R = chnkr_R.weights; wts_R = wts_R(:);
wts = chnkr_L.weights; wts = wts(:);



fkern_SR = @(s,t) - comp_obs(s,t,sys_R,'s');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
sysmats_R_obs = chunkermat(chnkr_R,fkern_SR,opts);
% sysmats_R = fkern_SR(chnkr_R,chnkr_R).*wts_R(:).';
t2 = toc(start)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   .  .  .  derivative double layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fkern_Dp_R = @(s,t) - comp_obs(s,t,sys_R,'dp');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
sysmatdd_R_obs = chunkermat(chnkr_R,fkern_Dp_R,opts);
% sysmatdd_R = fkern_Dp_R(chnkr_R,chnkr_R).*wts_R(:).';
t2 = toc(start)
% %%
fkern_D_R = @(s,t) - comp_obs(s,t,sys_R,'d');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
sysmatd_R_obs = chunkermat(chnkr_R,fkern_D_R,opts);
% sysmatd_R = fkern_Dp_R(chnkr_R,chnkr_R).*wts_R(:).';
t2 = toc(start)


fkern_Sp_R = @(s,t) - comp_obs(s,t,sys_R,'sp');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
sysmatsd_R_obs = chunkermat(chnkr_R,fkern_Sp_R,opts);
% sysmatsd_R = fkern_Sp_R(chnkr_R,chnkr_R).*wts_R(:).';
t2 = toc(start)


sysmat_R_obs = [sysmatd_R_obs,sysmats_R_obs;-sysmatdd_R_obs,-sysmatsd_R_obs];
%%

sysmat = [sysmat_L,sysmat_LR;sysmat_RL,sysmat_R+sysmat_R_obs];
sysmat = sysmat-eye(size(sysmat,1));


%%
src = [];
src.r = [-31;0];
src.n = [1;0];
tic;
rhs1 = [d_guide_kern(src,chnkr_L,sys_L,'s');-d_guide_kern(src,chnkr_L,sys_L,'sp')];
rhs2 = zeros(2*chnkr_R.npt,1);
rhs = [rhs1;rhs2];
toc
tic;
soln = sysmat\rhs;
ddens = soln([1:chnkr_L.npt,2*chnkr_L.npt+(1:chnkr_R.npt)]);
sdens = soln([chnkr_L.npt+(1:chnkr_R.npt),2*chnkr_L.npt+chnkr_R.npt+(1:chnkr_R.npt)]);
toc

%%
uin = zeros(ntarg,1);
uscat = zeros(ntarg,1);

uin(il) = d_guide_kern(src,targl,sys_L,'s');
%%

w = chnkr.weights; w = w(:);
ssrc = [];
ssrc.r = chnkr.r(:,:);
ssrc.n = chnkr.n(:,:);
ssrc.charge = sdens.*w;

dsrc = [];
dsrc.r = chnkr.r(:,:);
dsrc.n = chnkr.n(:,:);
dsrc.charge = ddens.*w;

id_L = 1:chnkr_L.npt;
ssrc_L = [];
ssrc_L.r = ssrc.r(:,id_L);
ssrc_L.n = ssrc.n(:,id_L);
ssrc_L.charge = ssrc.charge(id_L);

dsrc_L = [];
dsrc_L.r = dsrc.r(:,id_L);
dsrc_L.n = dsrc.n(:,id_L);
dsrc_L.charge = dsrc.charge(id_L);

id_R = chnkr_L.npt+(1:chnkr_R.npt);
ssrc_R = [];
ssrc_R.r = ssrc.r(:,id_R);
ssrc_R.n = ssrc.n(:,id_R);
ssrc_R.chrg = ssrc.charge(id_R);

dsrc_R = [];
dsrc_R.r = dsrc.r(:,id_R);
dsrc_R.n = dsrc.n(:,id_R);
dsrc_R.chrg = dsrc.charge(id_R);

tic
uscat(il) = d_guide_kern_lump(ssrc_L,dsrc_L,targl,sys_L);
uscat(im) = d_two_flat_kern_lump(ssrc,dsrc,targm,sys_M);
uscat(ir) = wave_kern.lump_source(targr,ssrc_R,dsrc_R,sys_R);
toc

%%

figure(7); clf
% subplot(1,2,1)
tools.pplot(X,Y,reshape(imag(uin-uscat),nplt,nplt),'$\Im(u_{\rm{tot}})$')
chnkr_L_p = real(sys_L.chnkr.r(:,:)); chnkr_L_p(:,chnkr_L_p(1,:)>chnkr_L.r(1,1)) = NaN;
% chnkr_M_p = real(sys_M.pts.r); chnkr_M_p(:,abs(chnkr_M_p(1,:))>1) = NaN;
chnkr_M_p = real([cos(sys_M.chnkr.r(1,:));sys_M.chnkr.r(2,:)]); chnkr_M_p(:,chnkr_M_p(1,:)<0) = NaN;
chnkr_R_p = real(sys_R.chnkr.r(:,:)); chnkr_R_p(:,chnkr_R_p(1,:)<chnkr_R.r(1,1)) = NaN;
    hold on
    plot(real(chnkr_L_p(1,:)),real(chnkr_L_p(2,:)),'k.','markersize',1)
    plot(real(chnkr_M_p(1,:)),real(chnkr_M_p(2,:)),'k.','markersize',1)
    plot(real(chnkr_R_p(1,:)),real(chnkr_R_p(2,:)),'k.','markersize',1)
    plot(real(chnkr.r(1,:)),real(chnkr.r(2,:)),'r.','markersize',1e-3)
    hold off

ti = toc(istart)

uin_naught = uin;
uscat_naught = uscat;
utot_naught = uin - uscat;


uin_naught_ctr = 0;
uscat_naught_ctr = wave_kern.lump_source(struct("r",ctr),ssrc_R,dsrc_R,sys_R);
utot_naught_ctr = uin_naught_ctr - uscat_naught_ctr;

%%
src = [];
src.r = [-31;0];
src.n = [1;0];

tic;
rhs1 = zeros(2*chnkr_L.npt,1);
rhs2 = -[chnk.helm2d.kern(sys.zk1,struct("r",ctr),chnkr_R,'s',1);-chnk.helm2d.kern(sys.zk1,struct("r",ctr),chnkr_R,'sprime',1)];
rhs = [rhs1;rhs2];
toc
tic;
soln = sysmat\rhs;
ddens = soln([1:chnkr_L.npt,2*chnkr_L.npt+(1:chnkr_R.npt)]);
sdens = soln([chnkr_L.npt+(1:chnkr_R.npt),2*chnkr_L.npt+chnkr_R.npt+(1:chnkr_R.npt)]);
toc


%%
uin = zeros(ntarg,1);
uscat = zeros(ntarg,1);

uin(ir) = chnk.helm2d.kern(sys.zk1,struct("r",ctr),targr,'s',1);
%%

w = chnkr.weights; w = w(:);
ssrc = [];
ssrc.r = chnkr.r(:,:);
ssrc.n = chnkr.n(:,:);
ssrc.charge = sdens.*w;

dsrc = [];
dsrc.r = chnkr.r(:,:);
dsrc.n = chnkr.n(:,:);
dsrc.charge = ddens.*w;

id_L = 1:chnkr_L.npt;
ssrc_L = [];
ssrc_L.r = ssrc.r(:,id_L);
ssrc_L.n = ssrc.n(:,id_L);
ssrc_L.charge = ssrc.charge(id_L);

dsrc_L = [];
dsrc_L.r = dsrc.r(:,id_L);
dsrc_L.n = dsrc.n(:,id_L);
dsrc_L.charge = dsrc.charge(id_L);

id_R = chnkr_L.npt+(1:chnkr_R.npt);
ssrc_R = [];
ssrc_R.r = ssrc.r(:,id_R);
ssrc_R.n = ssrc.n(:,id_R);
ssrc_R.chrg = ssrc.charge(id_R);

dsrc_R = [];
dsrc_R.r = dsrc.r(:,id_R);
dsrc_R.n = dsrc.n(:,id_R);
dsrc_R.chrg = dsrc.charge(id_R);

tic
uscat(il) = d_guide_kern_lump(ssrc_L,dsrc_L,targl,sys_L);
uscat(im) = d_two_flat_kern_lump(ssrc,dsrc,targm,sys_M);
uscat(ir) = wave_kern.lump_source(targr,ssrc_R,dsrc_R,sys_R);
toc

%%

figure(6); clf
% subplot(1,2,1)
tools.pplot(X,Y,reshape(imag(uin-uscat),nplt,nplt),'$\Im(u_{\rm{tot}})$')
chnkr_L_p = real(sys_L.chnkr.r(:,:)); chnkr_L_p(:,chnkr_L_p(1,:)>chnkr_L.r(1,1)) = NaN;
% chnkr_M_p = real(sys_M.pts.r); chnkr_M_p(:,abs(chnkr_M_p(1,:))>1) = NaN;
chnkr_M_p = real([cos(sys_M.chnkr.r(1,:));sys_M.chnkr.r(2,:)]); chnkr_M_p(:,chnkr_M_p(1,:)<0) = NaN;
chnkr_R_p = real(sys_R.chnkr.r(:,:)); chnkr_R_p(:,chnkr_R_p(1,:)<chnkr_R.r(1,1)) = NaN;
    hold on
    plot(real(chnkr_L_p(1,:)),real(chnkr_L_p(2,:)),'k.','markersize',1)
    plot(real(chnkr_M_p(1,:)),real(chnkr_M_p(2,:)),'k.','markersize',1)
    plot(real(chnkr_R_p(1,:)),real(chnkr_R_p(2,:)),'k.','markersize',1)
    plot(real(chnkr.r(1,:)),real(chnkr.r(2,:)),'r.','markersize',1e-3)
    hold off

ti = toc(istart);

uin_pt = uin;
uscat_pt = uscat;
utot_pt = uin - uscat;
end


%%
figure(8); clf
subplot(1,3,1)
tools.pplot(X,Y,reshape(imag(utot_blob-utot_naught),nplt,nplt),'$\Im(u_{\rm{tot}})$')
chnkr_L_p = real(sys_L.chnkr.r(:,:)); chnkr_L_p(:,chnkr_L_p(1,:)>chnkr_L.r(1,1)) = NaN;
% chnkr_M_p = real(sys_M.pts.r); chnkr_M_p(:,abs(chnkr_M_p(1,:))>1) = NaN;
chnkr_M_p = real([cos(sys_M.chnkr.r(1,:));sys_M.chnkr.r(2,:)]); chnkr_M_p(:,chnkr_M_p(1,:)<0) = NaN;
chnkr_R_p = real(sys_R.chnkr.r(:,:)); chnkr_R_p(:,chnkr_R_p(1,:)<chnkr_R.r(1,1)) = NaN;
hold on
plot(real(chnkr_L_p(1,:)),real(chnkr_L_p(2,:)),'k.','markersize',1)
plot(real(chnkr_M_p(1,:)),real(chnkr_M_p(2,:)),'k.','markersize',1)
plot(real(chnkr_R_p(1,:)),real(chnkr_R_p(2,:)),'k.','markersize',1)
plot(real(chnkr.r(1,:)),real(chnkr.r(2,:)),'r.','markersize',1e-3)
hold off
clim([-0.1,0.1])


% A = [besselj(0,zk2*d_R), -besselh(0,zk*d_R); -zk2*besselj(1,zk2*d_R), zk*besselh(1,zk*d_R)];
% cs = A\[besselj(0,zk*d_R);-zk*besselj(1,zk*d_R)];
% beta = cs(2);

beta = -(zk2*besselj(1,zk2*d_R)*besselj(0,zk*d_R)-zk*besselj(0,zk2*d_R)*besselj(1,zk*d_R)) / ...
    (zk2*besselh(0,zk*d_R)*besselj(1,zk2*d_R)-zk*besselh(1,zk*d_R)*besselj(0,zk2*d_R));

foldy_rat = beta*utot_naught_ctr/besselj(0,zk*d_R)*4i;

subplot(1,3,2)
tools.pplot(X,Y,reshape(imag(foldy_rat*utot_pt),nplt,nplt),'$\Im(u_{\rm{tot}})$')
chnkr_L_p = real(sys_L.chnkr.r(:,:)); chnkr_L_p(:,chnkr_L_p(1,:)>chnkr_L.r(1,1)) = NaN;
% chnkr_M_p = real(sys_M.pts.r); chnkr_M_p(:,abs(chnkr_M_p(1,:))>1) = NaN;
chnkr_M_p = real([cos(sys_M.chnkr.r(1,:));sys_M.chnkr.r(2,:)]); chnkr_M_p(:,chnkr_M_p(1,:)<0) = NaN;
chnkr_R_p = real(sys_R.chnkr.r(:,:)); chnkr_R_p(:,chnkr_R_p(1,:)<chnkr_R.r(1,1)) = NaN;
hold on
plot(real(chnkr_L_p(1,:)),real(chnkr_L_p(2,:)),'k.','markersize',1)
plot(real(chnkr_M_p(1,:)),real(chnkr_M_p(2,:)),'k.','markersize',1)
plot(real(chnkr_R_p(1,:)),real(chnkr_R_p(2,:)),'k.','markersize',1)
plot(real(chnkr.r(1,:)),real(chnkr.r(2,:)),'r.','markersize',1e-3)
hold off
clim([-0.1,0.1])

subplot(1,3,3)
tools.pplot(X,Y,reshape(imag(foldy_rat*utot_pt+utot_blob-utot_naught),nplt,nplt),'$\Im(u_{\rm{tot}})$')
chnkr_L_p = real(sys_L.chnkr.r(:,:)); chnkr_L_p(:,chnkr_L_p(1,:)>chnkr_L.r(1,1)) = NaN;
% chnkr_M_p = real(sys_M.pts.r); chnkr_M_p(:,abs(chnkr_M_p(1,:))>1) = NaN;
chnkr_M_p = real([cos(sys_M.chnkr.r(1,:));sys_M.chnkr.r(2,:)]); chnkr_M_p(:,chnkr_M_p(1,:)<0) = NaN;
chnkr_R_p = real(sys_R.chnkr.r(:,:)); chnkr_R_p(:,chnkr_R_p(1,:)<chnkr_R.r(1,1)) = NaN;
hold on
plot(real(chnkr_L_p(1,:)),real(chnkr_L_p(2,:)),'k.','markersize',1)
plot(real(chnkr_M_p(1,:)),real(chnkr_M_p(2,:)),'k.','markersize',1)
plot(real(chnkr_R_p(1,:)),real(chnkr_R_p(2,:)),'k.','markersize',1)
plot(real(chnkr.r(1,:)),real(chnkr.r(2,:)),'r.','markersize',1e-3)
hold off
clim([-0.1,0.1])

%%
figure(9)
u1 = utot_blob-utot_naught;
u2 = utot_pt;
rat = norm(u1)/norm(u2);

tools.pplot(X,Y,reshape(abs(u1)-rat*abs(u2),nplt,nplt),'$\Im(u_{\rm{tot}})$')
chnkr_L_p = real(sys_L.chnkr.r(:,:)); chnkr_L_p(:,chnkr_L_p(1,:)>chnkr_L.r(1,1)) = NaN;
% chnkr_M_p = real(sys_M.pts.r); chnkr_M_p(:,abs(chnkr_M_p(1,:))>1) = NaN;
chnkr_M_p = real([cos(sys_M.chnkr.r(1,:));sys_M.chnkr.r(2,:)]); chnkr_M_p(:,chnkr_M_p(1,:)<0) = NaN;
chnkr_R_p = real(sys_R.chnkr.r(:,:)); chnkr_R_p(:,chnkr_R_p(1,:)<chnkr_R.r(1,1)) = NaN;
hold on
plot(real(chnkr_L_p(1,:)),real(chnkr_L_p(2,:)),'k.','markersize',1)
plot(real(chnkr_M_p(1,:)),real(chnkr_M_p(2,:)),'k.','markersize',1)
plot(real(chnkr_R_p(1,:)),real(chnkr_R_p(2,:)),'k.','markersize',1)
plot(real(chnkr.r(1,:)),real(chnkr.r(2,:)),'r.','markersize',1e-3)
hold off



return
%%


% planewave direction
phi = 0;
kvec = zk*[cos(phi);sin(phi)];


scale = 10;

zk = 1;
zk2 = 2*zk*scale;

d_R = 0.5/scale;

sys_R = build_blob_sys(2,d_R,ctr,zk,zk2);

% define boundary data
rhs_val = -planewave(kvec,sys_R.chnkr.r(:,:));
rhs_grad = -1i*sum(kvec.*sys_R.chnkr.n(:,:),1).*rhs_val;


rhs = [rhs_val(:); rhs_grad(:)];

dens = sys_R.sysmat\rhs;



Lplot = 10;
nplt = 160;

xts = linspace(-Lplot+ctr(1),Lplot+ctr(1),nplt);
yts = linspace(-Lplot+ctr(2),Lplot+ctr(2),nplt);
[X,Y] = meshgrid(xts,yts);
xtargs = X(:).'; ytargs = Y(:).';

targ = [];
targ.r = [xtargs;ytargs];
targ.n = [1;0]+0*targ.r;

ntarg = size(targ.r,2);
iiin = chunkerinterior(sys_R.chnkr,targ.r);
iout = ~iiin;

targin = [];
targin.r = targ.r(:,iiin);
targin.n = targ.n(:,iiin);

targout = [];
targout.r = targ.r(:,iout);
targout.n = targ.n(:,iout);

uin = zeros(1,ntarg);
uin(iout) = planewave(kvec(:),targout.r);

uscat = zeros(1,ntarg);
skern = @(s,t) chnk.helm2d.kern(zk,s,t,"s");
dkern = @(s,t) chnk.helm2d.kern(zk,s,t,"d");
uscat(iout) = chunkerkerneval(sys_R.chnkr,dkern,dens(1:sys_R.chnkr.npt),targout.r) + ...
     chunkerkerneval(sys_R.chnkr,skern,dens(sys_R.chnkr.npt+(1:sys_R.chnkr.npt)),targout.r);

skern = @(s,t) chnk.helm2d.kern(zk2,s,t,"s");
dkern = @(s,t) chnk.helm2d.kern(zk2,s,t,"d");
uscat(iiin) = chunkerkerneval(sys_R.chnkr,dkern,dens(1:sys_R.chnkr.npt),targin.r) + ...
     chunkerkerneval(sys_R.chnkr,skern,dens(sys_R.chnkr.npt+(1:sys_R.chnkr.npt)),targin.r);

% %%
utot = uin + uscat;
figure(9);clf
subplot(1,3,1)
tools.pplot(X,Y,reshape(real(uin),nplt,nplt),'$\Im(u_{\rm{in}})$')

subplot(1,3,2)
tools.pplot(X,Y,reshape(real(uscat),nplt,nplt),'$\Im(u_{\rm{scat}})$')

subplot(1,3,3)
tools.pplot(X,Y,reshape(real(utot),nplt,nplt),'$\Im(u_{\rm{tot}})$')

% %%
A = [besselj(0,zk2*d_R), -besselh(0,zk*d_R); -zk2*besselj(1,zk2*d_R), zk*besselh(1,zk*d_R)];
cs = A\[besselj(0,zk*d_R);-zk*besselj(1,zk*d_R)];
beta = cs(2);

foldy_rat = beta*planewave(kvec(:),ctr)/besselj(0,zk*d_R)*4i;
u_pt = chnk.helm2d.kern(zk,struct("r",ctr,"n",[1;0]),targ,"s")*foldy_rat;

M = max(abs(u_pt));

figure(10);clf
subplot(1,3,1)
tools.pplot(X,Y,reshape(real(uscat),nplt,nplt),'$\Im(u_{\rm{scat}})$')
clim([-M,M])

subplot(1,3,2)
tools.pplot(X,Y,reshape(real(u_pt),nplt,nplt),'$\Im(u_{\rm{pt}})$')
clim([-M,M])

uscat(iiin) =0;
u_pt(iiin) = 0;
subplot(1,3,3)
tools.pplot(X,Y,reshape(real(uscat(:)+u_pt(:)),nplt,nplt),'$\Im(u_{\rm{diff}})$')
% clim([-M,M])
