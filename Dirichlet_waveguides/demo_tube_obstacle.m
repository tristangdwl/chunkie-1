addpath('../')
addpaths_loc();


zk = 1;
test_acc = 1;

which_test = 0;

d_L = 8;
d_R = 1;
ctr = [6;0];
L = 40;

nover = 2;

sys_L = build_dirichlet_sys(2,d_L,L,zk);
sys_M = build_two_flats_sys(3,d_L,zk);
sys_R = build_blob_sys(2,d_R,ctr,zk,2*zk);

L_skel= L;
tol = 1e-10;
sys_R.comp_data = comp_far(sys_R,L_skel,tol);

chnkr = vert_chnkr(1,L);
chnkr_L = move(chnkr,[0;0]);
chnkr_L.n = normals(chnkr_L);
wts_L = chnkr_L.weights; wts_L = wts_L(:);

chnkr_R = move(chnkr,[],[4;0]);
chnkr_R.n = normals(chnkr_R);
wts_R = chnkr_R.weights; wts_R = wts_R(:);

chnkr = merge([chnkr_L,chnkr_R]);
chnkr.n = normals(chnkr);
wts = chnkr_L.weights; wts = wts(:);
chnkr.npt


irebuild = 0;
if irebuild
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   .  .  .  single double layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fkern_SL = @(s,t) d_guide_kern(s,t,sys_L,'s')-d_two_flat(s,t,sys_M,'s');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
sysmats_L = chunkermat(chnkr_L,fkern_SL,opts);
% sysmats_L = fkern_SL(chnkr_L,chnkr_L).*wts_L(:).';
t1 = toc(start)

fkern_SR = @(s,t) d_two_flat(s,t,sys_M,'s');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; sysmats_R = chunkermat(chnkr_R,fkern_SR,opts);
% sysmats_R = fkern_SR(chnkr_R,chnkr_R).*wts_R(:).';
t2 = toc(start)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   .  .  .  derivative double layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fkern_Dp_L = @(s,t) d_guide_kern(s,t,sys_L,'dp')-d_two_flat(s,t,sys_M,'dp');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
sysmatdd_L = chunkermat(chnkr_L,fkern_Dp_L,opts);
% sysmatdd_L = fkern_Dp_L(chnkr_L,chnkr_L).*wts_L(:).';
t1 = toc(start)

fkern_Dp_R = @(s,t) d_two_flat(s,t,sys_M,'dp');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
sysmatdd_R = chunkermat(chnkr_R,fkern_Dp_R,opts);
% sysmatdd_R = fkern_Dp_R(chnkr_R,chnkr_R).*wts_R(:).';
t2 = toc(start)
%%
start = tic; 
fkern_Sp_M = @(s,t) -d_two_flat(s,t,sys_M,'sp');
% sysmatsd_L = chunkermat(chnkr_L,fkern_Sp_M,opts);
sysmatsd_R = -chunkermat(chnkr_R,fkern_Sp_M,opts);
fkern_D_M = @(s,t) -d_two_flat(s,t,sys_M,'d');
% sysmatd_L = chunkermat(chnkr_L,fkern_D_M,opts);
sysmatd_R = -chunkermat(chnkr_R,fkern_D_M,opts);
t3 = toc(start)
sysmatsd_L = zeros(chnkr_L.npt);
sysmatd_L = zeros(chnkr_L.npt);
sysmat_L = [sysmatd_L,sysmats_L;-sysmatdd_L,-sysmatsd_L];
sysmat_R = [sysmatd_R,sysmats_R;-sysmatdd_R,-sysmatsd_R];

%%
w_L = chnkr_L.weights; w_L = w_L(:).';
w_R = chnkr_R.weights; w_R = w_R(:).';
fkern_S_M = @(s,t) -d_two_flat(s,t,sys_M,'s');
fkern_Dp_M = @(s,t) -d_two_flat(s,t,sys_M,'dp');

pts_L = [];
pts_L.r = chnkr_L.r(:,:);
pts_L.n = chnkr_L.n(:,:);

pts_R = [];
pts_R.r = chnkr_R.r(:,:);
pts_R.n = chnkr_R.n(:,:);

start = tic; 
sysmats_RL = fkern_S_M(pts_L,pts_R).*w_L;
sysmatsd_RL = fkern_Sp_M(pts_L,pts_R).*w_L;
t4 = toc(start)
start = tic; 
sysmatd_RL = fkern_D_M(pts_L,pts_R).*w_L;
sysmatdd_RL = fkern_Dp_M(pts_L,pts_R).*w_L;
t5 = toc(start)

start = tic; 
sysmats_LR = fkern_S_M(pts_R,pts_L).*w_R;
sysmatsd_LR = fkern_Sp_M(pts_R,pts_L).*w_R;
t4 = toc(start)
start = tic; 
sysmatd_LR = fkern_D_M(pts_R,pts_L).*w_R;
sysmatdd_LR = fkern_Dp_M(pts_R,pts_L).*w_R;
t5 = toc(start)

sysmat_RL = -[sysmatd_RL,sysmats_RL;-sysmatdd_RL,-sysmatsd_RL];
sysmat_LR = [sysmatd_LR,sysmats_LR;-sysmatdd_LR,-sysmatsd_LR];

%%


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

occ = 1024;
rank_or_tol = 1e-8;

xxflam = repmat(real(chnkr.r(1,:)),1,2);
xyflam = repmat(real(chnkr.r(2,:)),1,2);
xflam = [xxflam(:).';xyflam(:).'];
tic;
F = rskelf(sysmat,xflam,occ,rank_or_tol); 
tflam = toc
sysinv = @(rhs) rskelf_sv(F,rhs);
end
%%
src = [];
% src.r = [-5;5];
src.r = [-31;0];
src.n = [1;0];

% src_R = [];
% src_R.r = [21;0];
% src_R.n = [1;0];
tic;
rhs1 = [d_guide_kern(src,chnkr_L,sys_L,'s');-d_guide_kern(src,chnkr_L,sys_L,'sp')];
rhs2 = zeros(2*chnkr_R.npt,1);
rhs = [rhs1;rhs2];
toc
tic;
% soln = sysinv(rhs);
soln = sysmat\rhs;
ddens = soln([1:chnkr_L.npt,2*chnkr_L.npt+(1:chnkr_R.npt)]);
sdens = soln([chnkr_L.npt+(1:chnkr_R.npt),2*chnkr_L.npt+chnkr_R.npt+(1:chnkr_R.npt)]);
toc

%%
figure(2);clf
plot(real(chnkr.r(2,:)),imag(ddens),'.')
hold on
plot(real(chnkr.r(2,:)),imag(sdens),'.')
hold off
% 





%%
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
uin = zeros(ntarg,1);
uscat = zeros(ntarg,1);

uin(il) = d_guide_kern(src,targl,sys_L,'s');
% uin(ir) = d_guide_kern(src_R,targr,sys_R,'s');
% uin(il) = d_guide_kern(src,targl,sys_L,'sp');
% uin(ir) = d_guide_kern(src_R,targr,sys_R,'sp');
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

% ssrc.charge = [ssrc_L.charge;0*ssrc_R.charge];
% dsrc.charge = [dsrc_L.charge;0*dsrc_R.charge];
tic
uscat(il) = d_guide_kern_lump(ssrc_L,dsrc_L,targl,sys_L);

uscat(im) = d_two_flat_kern_lump(ssrc,dsrc,targm,sys_M);
uscat(ir) = wave_kern.lump_source(targr,ssrc_R,dsrc_R,sys_R);
toc

%%
% figure(1); clf
% subplot(1,3,1)
% tools.pplot(X,Y,reshape(imag(uin),nplt,nplt),'$u_{\rm{in}}$')
% 
% subplot(1,3,2)
% tools.pplot(X,Y,reshape(imag(uscat),nplt,nplt),'$u_{\rm{scat}}$')
% 
% subplot(1,3,3)
% tools.pplot(X,Y,reshape(imag(uin-uscat),nplt,nplt),'imag($u_{\rm{tot}})$')

figure(5); 
subplot(1,2,1)
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

subplot(1,2,2)
tools.pplot(X,Y,reshape((abs((uin-uscat))),nplt,nplt),'$|u_{\rm{tot}}|$')
hold on
    plot(real(chnkr_L_p(1,:)),real(chnkr_L_p(2,:)),'k.','markersize',1)
    plot(real(chnkr_M_p(1,:)),real(chnkr_M_p(2,:)),'k.','markersize',1)
    plot(real(chnkr_R_p(1,:)),real(chnkr_R_p(2,:)),'k.','markersize',1)
    % plot(real(chnkr.r(1,:)),real(chnkr.r(2,:)),'r.','markersize',1e-3)
    hold off

font_rate=10/12;f_width = 3;f_height = f_width*0.75;
set(gcf,'Position',[100   200   2.5*round(f_width*font_rate*144)   1.4*round(f_height*font_rate*144)])

    % exportgraphics(gcf,'Two_Dirichlet_tubes.pdf')

    %%
%      figure(3)
% h = pcolor(log10(abs(sysmat))); set(h,'edgecolor','none');
% colorbar
% title('$\log_{10}|K_{ij}|$','Interpreter','latex')
% set(gca,'fontsize',16)
% daspect([1,1,1])
% set(gca,'TickLabelInterpreter','latex')
% exportgraphics(gcf,'d_waveguide_sysmat.pdf')


if test_acc
src = [];
% src.r = [-5;5];
src.r = [8;0];
src.n = [1;0];
tic;
rhs1 = [d_guide_kern(src,chnkr_L,sys_L,'s');-d_guide_kern(src,chnkr_L,sys_L,'sp')];
% rhs2 = -0*[d_guide_kern(src_R,chnkr_R,sys_R,'s');-d_guide_kern(src_R,chnkr_R,sys_R,'sp')];
rhs2 = zeros(2*chnkr_R.npt,1);
rhs = [rhs1;rhs2];
toc
tic;
% soln = sysinv(rhs);
soln = sysmat\rhs;
ddens = soln([1:chnkr_L.npt,2*chnkr_L.npt+(1:chnkr_R.npt)]);
sdens = soln([chnkr_L.npt+(1:chnkr_R.npt),2*chnkr_L.npt+chnkr_R.npt+(1:chnkr_R.npt)]);
toc

%%
uin2 = zeros(ntarg,1);
uscat2 = zeros(ntarg,1);
uin2(il) = d_guide_kern(src,targl,sys_L,'s');

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

% ssrc.charge = [ssrc_L.charge;0*ssrc_R.charge];
% dsrc.charge = [dsrc_L.charge;0*dsrc_R.charge];
tic
uscat2(il) = d_guide_kern_lump(ssrc_L,dsrc_L,targl,sys_L);

uscat2(im) = d_two_flat_kern_lump(ssrc,dsrc,targm,sys_M);
uscat2(ir) = wave_kern.lump_source(targr,ssrc_R,dsrc_R,sys_R);
toc

%%
figure(6); clf

tools.pplot(X,Y,reshape(log10(abs(uin2-uscat2)),nplt,nplt),'$|u_{\rm{tot}}|$')
hold on
    plot(real(chnkr_L_p(1,:)),real(chnkr_L_p(2,:)),'k.','markersize',1)
    plot(real(chnkr_M_p(1,:)),real(chnkr_M_p(2,:)),'k.','markersize',1)
    plot(real(chnkr_R_p(1,:)),real(chnkr_R_p(2,:)),'k.','markersize',1)
    % plot(real(chnkr.r(1,:)),real(chnkr.r(2,:)),'r.','markersize',1e-3)
    hold off
font_rate=10/12;f_width = 3;f_height = f_width*0.75;
set(gcf,'Position',[100   200   1.4*round(f_width*font_rate*144)   1.4*round(f_height*font_rate*144)])
% exportgraphics(gcf,'two_dirichlet_error.pdf')

end