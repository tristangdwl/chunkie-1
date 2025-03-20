addpath(genpath('~/FastAlgs_backup/chunkie-1'))
addpaths_loc();
load("two_tube_data.mat")

zk = sys_L.zk;


d_R = 1;
sys = build_blob_sys(2,d_R,[0;0],zk,2*zk);
% ctr_x = [6,8,10,20];
% ctrs = [ctr_x;0*ctr_x];
ctrs = {[6;0],[8,8;4,-4],[6,12;0,0],[6,6,10;4,-4,0]};

L_skel= 40;
tol = 1e-10;

nctr = size(ctrs,2);
nctr = length(ctrs);
figure(5); clf;
t = tiledlayout("flow");t.TileSpacing = 'tight';

% sys.comp_data = comp_far(sys,L_skel,tol);

for i = 1:nctr
    istart = tic;
sys_R = build_blob_sys(2,d_R,ctrs{i},zk,2*zk);
% sys_R = build_empty_sys(zk);
% sys_R.chnkr = move(sys_R.chnkr,[],ctrs(:,i));
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
% 
% occ = 1024;
% rank_or_tol = 1e-8;
% 
% xxflam = repmat(real(chnkr.r(1,:)),1,2);
% xyflam = repmat(real(chnkr.r(2,:)),1,2);
% xflam = [xxflam(:).';xyflam(:).'];
% tic;
% F = rskelf(sysmat,xflam,occ,rank_or_tol); 
% tflam = toc
% sysinv = @(rhs) rskelf_sv(F,rhs);

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
% figure(2);clf
% plot(real(chnkr.r(2,:)),imag(ddens),'.')
% hold on
% plot(real(chnkr.r(2,:)),imag(sdens),'.')
% hold off
% % 





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
nexttile();
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

% subplot(1,2,2)
% tools.pplot(X,Y,reshape((abs((uin-uscat))),nplt,nplt),'$|u_{\rm{tot}}|$')
% hold on
%     plot(real(chnkr_L_p(1,:)),real(chnkr_L_p(2,:)),'k.','markersize',1)
%     plot(real(chnkr_M_p(1,:)),real(chnkr_M_p(2,:)),'k.','markersize',1)
%     plot(real(chnkr_R_p(1,:)),real(chnkr_R_p(2,:)),'k.','markersize',1)
%     % plot(real(chnkr.r(1,:)),real(chnkr.r(2,:)),'r.','markersize',1e-3)
%     hold off

% font_rate=10/12;f_width = 3;f_height = f_width*0.75;
% set(gcf,'Position',[100   200   2.5*round(f_width*font_rate*144)   1.4*round(f_height*font_rate*144)])
ti = toc(istart)
end
%%
% exportgraphics(gcf,'various_obs_locations.pdf','resolution',300)
 exportgraphics(gcf,'various_obs_nums.pdf','resolution',300)