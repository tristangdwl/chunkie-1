addpath(genpath('~/FastAlgs_backup/chunkie-1'))
addpaths_loc();

zk = 1;
zk2 = 20*zk;

irebuild = 1; % rebuild the systam matrix?
d_L = 8;
L_chnkr = 40;

nover = 2;

sys_L = build_dirichlet_sys(2,d_L,L_chnkr,zk);
sys_R = build_two_flats_sys(3,d_L,zk);


d_R = 0.1;
ctr = [6;0];
% sys_blob = build_blob_sys(2,d_R,[0;0],zk,2*zk);


chnkr = vert_chnkr(1,L_chnkr);
chnkr = sort(chnkr);
% wts = chnkr.weights; wts = wts(:);
% chnk_struc = [];
% chnk_struc.r = chnkr.r(:,:); chnk_struc.n = chnkr.n(:,:);

if irebuild
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   .  .  .  single double layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
fkern_S = @(s,t) d_guide_kern(s,t,sys_L,'s')-d_two_flat(s,t,sys_R,'s');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
sysmats = chunkermat(chnkr,fkern_S,opts);
% sysmats = fkern_S(chnk_struc,chnk_struc).*wts(:).';
t1 = toc(start)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   .  .  .  derivative double layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fkern_Sp = @(s,t) d_guide_kern(s,t,sys_L,'sp')-d_two_flat(s,t,sys_R,'sp');

fkern_Dp = @(s,t) d_guide_kern(s,t,sys_L,'dp')-d_two_flat(s,t,sys_R,'dp');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
sysmatdd = chunkermat(chnkr,fkern_Dp,opts);
% sysmatdd = fkern_Dp(chnk_struc,chnk_struc).*wts(:).';
t1 = toc(start)
sysmatd = zeros(chnkr.npt);
sysmatsd = zeros(chnkr.npt); 


sysmat = [sysmatd,sysmats;-sysmatdd,-sysmatsd];
sysmat = sysmat-eye(size(sysmat,1));

end


%% set up targs

Lplot = 26;
nplt = 160;

xts = linspace(-Lplot,Lplot,nplt);
yts = linspace(-Lplot,Lplot,nplt);
[X,Y] = meshgrid(xts,yts);
xtargs = X(:).'; ytargs = Y(:).';

targ = [];
targ.r = [xtargs;ytargs];
targ.n = [1;0]+0*targ.r;

il = (targ.r(1,:)<chnkr.r(1,1));
ir = ~il;

targl = [];
targl.r = targ.r(:,il);
targl.n = targ.n(:,il);

targr = [];
targr.r = targ.r(:,ir);
targr.n = targ.n(:,ir);

ntarg = size(targ.r,2);


%% unperturbed solution
src = [];
src.r = [-31;0];
src.n = [1;0];
tic;
rhs = [d_guide_kern(src,chnkr,sys_L,'s');-d_guide_kern(src,chnkr,sys_L,'sp')];
toc
tic;
soln = sysmat\rhs;
ddens = soln(1:chnkr.npt);
sdens = soln(chnkr.npt+(1:chnkr_L.npt));
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

tic
uscat(il) = d_guide_kern_lump(ssrc,dsrc,targl,sys_L);
uscat(ir) = d_two_flat_kern_lump(ssrc,dsrc,targr,sys_R);
toc

%%

figure(7); clf
% subplot(1,2,1)
tools.pplot(X,Y,reshape(imag(uin-uscat),nplt,nplt),'$\Im(u_{\rm{tot}})$')
chnkr_L_p = real(sys_L.chnkr.r(:,:)); chnkr_L_p(:,chnkr_L_p(1,:)>chnkr.r(1,1)) = NaN;
chnkr_R_p = real(sys_R.chnkr.r(:,:)); chnkr_R_p(:,chnkr_R_p(1,:)<chnkr.r(1,1)) = NaN;
    hold on
    plot(real(chnkr_L_p(1,:)),real(chnkr_L_p(2,:)),'k.','markersize',1)
    plot(real(chnkr_R_p(1,:)),real(chnkr_R_p(2,:)),'k.','markersize',1)
    plot(real(chnkr.r(1,:)),real(chnkr.r(2,:)),'r.','markersize',1e-3)
    hold off


uin_naught = uin;
uscat_naught = uscat;
utot_naught = uin - uscat;


uin_naught_ctr = 0;
uscat_naught_ctr = d_two_flat_kern_lump(ssrc,dsrc,struct("r",ctr,"n",[1;0]),sys_R);
utot_naught_ctr = uin_naught_ctr - uscat_naught_ctr;

%% Foldy-Lax solution
src = [];
src.r = [-31;0];
src.n = [1;0];

tic;
rhs = -[d_two_flat(struct("r",ctr),chnkr,sys_R,'s');-d_two_flat(struct("r",ctr),targr,sys_R,'sp')];
toc
tic;
soln = sysmat\rhs;
ddens = soln(1:chnkr.npt);
sdens = soln(chnkr.npt+(1:chnkr.npt));
toc


%%
uin = zeros(ntarg,1);
uscat = zeros(ntarg,1);

% uin(ir) = chnk.helm2d.kern(sys.zk1,struct("r",ctr),targr,'s',1);
uin(ir) = d_two_flat(struct("r",ctr),targr,sys_R,'s');
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

tic
uscat(il) = d_guide_kern_lump(ssrc,dsrc,targl,sys_L);
uscat(ir) = d_two_flat_kern_lump(ssrc,dsrc,targr,sys_R);
toc

%%

figure(6); clf
% subplot(1,2,1)
tools.pplot(X,Y,reshape(imag(uin-uscat),nplt,nplt),'$\Im(u_{\rm{tot}})$')
chnkr_L_p = real(sys_L.chnkr.r(:,:)); chnkr_L_p(:,chnkr_L_p(1,:)>chnkr.r(1,1)) = NaN;
chnkr_R_p = real(sys_R.chnkr.r(:,:)); chnkr_R_p(:,chnkr_R_p(1,:)<chnkr.r(1,1)) = NaN;
    hold on
    plot(real(chnkr_L_p(1,:)),real(chnkr_L_p(2,:)),'k.','markersize',1)
    plot(real(chnkr_R_p(1,:)),real(chnkr_R_p(2,:)),'k.','markersize',1)
    plot(real(chnkr.r(1,:)),real(chnkr.r(2,:)),'r.','markersize',1e-3)
    hold off

uin_pt = uin;
uscat_pt = uscat;
utot_pt = uin - uscat;


%%
figure(8); clf
A = [besselj(0,zk2*d_R), -besselh(0,zk*d_R); -zk2*besselj(1,zk2*d_R), zk*besselh(1,zk*d_R)];
cs = A\[besselj(0,zk*d_R);-zk*besselj(1,zk*d_R)];
beta = cs(2);

rats = beta*utot_naught_ctr/besselj(0,zk*d_R);

tools.pplot(X,Y,reshape(abs(utot_naught + rats*utot_pt),nplt,nplt),'$\Im(u_{\rm{tot}})$')
chnkr_L_p = real(sys_L.chnkr.r(:,:)); chnkr_L_p(:,chnkr_L_p(1,:)>chnkr.r(1,1)) = NaN;
chnkr_R_p = real(sys_R.chnkr.r(:,:)); chnkr_R_p(:,chnkr_R_p(1,:)<chnkr.r(1,1)) = NaN;
hold on
plot(real(chnkr_L_p(1,:)),real(chnkr_L_p(2,:)),'k.','markersize',1)
plot(real(chnkr_R_p(1,:)),real(chnkr_R_p(2,:)),'k.','markersize',1)
plot(real(chnkr.r(1,:)),real(chnkr.r(2,:)),'r.','markersize',1e-3)
hold off


