addpath(genpath('~/FastAlgs_backup/chunkie-1'))
addpaths_loc();

d = 8;
zk = 1;

nover_match = 1;
L_chnkr = 40;

% tube_sys = build_tube(d,zk,nover_match,L_chnkr);



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

il = (targ.r(1,:)<0);
ir = ~il;

targl = [];
targl.r = targ.r(:,il);
targl.n = targ.n(:,il);

targr = [];
targr.r = targ.r(:,ir);
targr.n = targ.n(:,ir);

ntarg = size(targ.r,2);


%% find and plot Green's function

src = []; src.r = [[4;0],[6;3],[7;-3]]; src.n = [1;0] + 0*src.r;


vals = tube_kern(src,targ,tube_sys);


figure(7); clf
% subplot(1,2,1)
tools.pplot(X,Y,reshape(imag(vals(:,3)),nplt,nplt),'$\Im(u_{\rm{tot}})$')
hold on
plot(real(tube_sys.pts(1,:)),real(tube_sys.pts(2,:)),'k.','markersize',1)
plot(real(tube_sys.chnkr.r(1,:)),real(tube_sys.chnkr.r(2,:)),'r.','markersize',1e-3)
hold off

