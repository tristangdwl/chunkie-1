zk = 1;
test_acc = 1;
d_L = 4;
% d_R = 8;
L = 40;

nover = 2;

sys_L = build_dirichlet_sys(2,d_L,L,zk);
sys_R = build_two_flats_sys(3,d_L,zk);

chnkr = vert_chnkr(1,L);
% chnkr = vert_chnkr2(2,d_L,L);
% cparams = []; cparams.ta = -d_L/2;cparams.tb = d_L/2;
% chnkr = chunkerfuncuni(@flatinterface,6,cparams);
% chnkr = move(chnkr,[],[],pi/2);

chnkr = sort(chnkr);
chnkr.npt
% nref = 15;
% for i = 1:nref
% 
%     % flag_u = min(abs(chnkr.r(2,1,:) -1));
%     opts_flag.fac = .8;
%     bdrys = [[0;d_L/2],[0;-d_L/2]];
%     flag = flagnear(chnkr,bdrys,opts_flag);
%     flag = sum(flag,1)>0; 
%     ns = 1:chnkr.nch;
%
%     opts_r = [];
%     opts_r.splitchunks = ns(flag);
% 
%     chnkr = refine(chnkr,opts_r);
% 
% end
% chnkr.n = normals(chnkr);
chnkr.npt

irebuild = 0;
if irebuild
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   .  .  .  single double layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fkern_S = @(s,t) d_guide_kern(s,t,sys_L,'s')-d_two_flat(s,t,sys_R,'s');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; sysmats = chunkermat(chnkr,fkern_S,opts);
t1 = toc(start)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   .  .  .  derivative double layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fkern_Sp = @(s,t) d_guide_kern(s,t,sys_L,'sp')-d_two_flat(s,t,sys_R,'sp');

fkern_Dp = @(s,t) d_guide_kern(s,t,sys_L,'dp')-d_two_flat(s,t,sys_R,'dp');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; sysmatdd = chunkermat(chnkr,fkern_Dp,opts);
t1 = toc(start)
sysmatd = zeros(chnkr.npt);
sysmatsd = zeros(chnkr.npt); 


sysmat = [sysmatd,sysmats;-sysmatdd,-sysmatsd];
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
src.r = [-21;0];
src.n = [1;0];
tic;
rhs = [d_guide_kern(src,chnkr,sys_L,'s');-d_guide_kern(src,chnkr,sys_L,'sp')];
toc
tic;
% soln = sysinv(rhs);
soln = sysmat\rhs;
ddens = soln(1:chnkr.npt);
sdens = soln((chnkr.npt+1):end);
toc

%%
figure(2);clf
plot(real(chnkr.r(2,:)),imag(ddens),'.')
hold on
plot(real(chnkr.r(2,:)),imag(sdens),'.')
hold off
% 
% %%
% figure(3);clf
% iplot = (chnkr.r(2,:)>-2) & (chnkr.r(2,:)<0);
% plot(real(chnkr.r(2,iplot)),(abs(ddens(iplot))),'.')
% hold on
% plot(real(chnkr.r(2,iplot)),(abs(sdens(iplot))),'.')
% hold off
% xlim([-2.0,-1.99])

%%
% figure(4);clf
% iplot = (chnkr.r(2,:)>-2) & (chnkr.r(2,:)<1);
% plot(log10(chnkr.r(2,iplot)+2),log10(abs(ddens(iplot))),'.')
% hold on
% plot(log10(chnkr.r(2,iplot)+2),log10(abs(sdens(iplot))),'.')
% hold off
% % xlim([-6.0,-3])
% % ylim([0,5])




%%

Lplot = 20;
nplt = 160;

xts = linspace(-Lplot,Lplot,nplt);
yts = linspace(-Lplot,Lplot,nplt);
[X,Y] = meshgrid(xts,yts);
xtargs = X(:).'; ytargs = Y(:).';

targ = [];
targ.r = [xtargs;ytargs];
targ.n = [1;0]+0*targ.r;

il = (targ.r(1,:)<0);
ir = (targ.r(1,:)>0);

targl = [];
targl.r = targ.r(:,il);
targl.n = targ.n(:,il);

targr = [];
targr.r = targ.r(:,ir);
targr.n = targ.n(:,ir);

ntarg = size(targ.r,2);
uin = zeros(ntarg,1);
uscat = zeros(ntarg,1);

uin(il) = d_guide_kern(src,targl,sys_L,'s');
% uin(ir) = d_guide_kern(src,targr,sys_R,'s');
% uin(ir) = chnk.helm2d.kern(sys_L.zk,src,targr,'s');

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
figure(1); clf
subplot(1,3,1)
tools.pplot(X,Y,reshape(imag(uin),nplt,nplt),'$u_{\rm{in}}$')

subplot(1,3,2)
tools.pplot(X,Y,reshape(imag(uscat),nplt,nplt),'$u_{\rm{scat}}$')

subplot(1,3,3)
tools.pplot(X,Y,reshape(imag(uin-uscat),nplt,nplt),'imag($u_{\rm{tot}})$')

figure(5); clf
subplot(1,2,1)
tools.pplot(X,Y,reshape(imag(uin-uscat),nplt,nplt),'imag($u_{\rm{tot}})$')
chnkrL = real(sys_L.chnkr.r(:,:)); chnkrL(:,chnkrL(1,:)>0) = NaN;
chnkrR = real([cos(sys_R.chnkr.r(1,:));sys_R.chnkr.r(2,:)]); chnkrR(:,chnkrR(1,:)<0) = NaN;
    hold on
    plot(real(chnkrL(1,:)),real(chnkrL(2,:)),'r.','markersize',1)
    plot(real(chnkrR(1,:)),real(chnkrR(2,:)),'r.','markersize',1)
    plot(real(chnkr.r(1,:)),real(chnkr.r(2,:)),'k.','markersize',.5)
    hold off

subplot(1,2,2)
tools.pplot(X,Y,reshape((abs(uin-uscat)),nplt,nplt),'$|u_{\rm{tot}}|$')
hold on
    plot(real(chnkrL(1,:)),real(chnkrL(2,:)),'r.','markersize',1)
    plot(real(chnkrR(1,:)),real(chnkrR(2,:)),'r.','markersize',1)
    plot(real(chnkr.r(1,:)),real(chnkr.r(2,:)),'k.','markersize',.5)
    hold off

    % exportgraphics(gcf,'open_dirichlet_waveguide.pdf')

    %%
     figure(3)
h = pcolor(log10(abs(sysmat))); set(h,'edgecolor','none');
colorbar
title('$\log_{10}|K_{ij}|$','Interpreter','latex')
set(gca,'fontsize',16)
daspect([1,1,1])
set(gca,'TickLabelInterpreter','latex')
% exportgraphics(gcf,'d_waveguide_sysmat.pdf')