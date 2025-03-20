zk = 10;
zk2 = 1;
zk2 = zk;
test_acc = 0;
d = 4;
L = 40;

do = 1;
ctr = [3;0];

nover = 4;
sys = build_end_obs_sys(nover,d,zk,do,ctr,zk2);
sys_f = build_two_flats_sys(nover,d,zk);

Lplot = 4;
nplt = 300;

xts = linspace(-Lplot,Lplot,nplt);
yts = linspace(-Lplot,Lplot,nplt);
[X,Y] = meshgrid(xts,yts);
xtargs = X(:).'; ytargs = Y(:).';

targ = [];
targ.r = [xtargs;ytargs];
targ.n = [1;0]+0*targ.r;
%%
src = [];
src.r = [-1;1.5];
% src.r = [-1;10];
src.n = [1;0];
val = d_end_obs(src,targ,sys,'d');
val2 = d_two_flat(src,targ,sys_f,'d');
%%
figure(1)
tools.pplot(X,Y,reshape(imag(val2),nplt,nplt),'u')
% hold on
% plot(sys.pts.r(1,:),sys.pts.r(2,:),'k.')
% plot(sys.chnkr_o.r(1,:),sys.chnkr_o.r(2,:),'k.')
% hold off
% tools.pplot(X,Y,reshape(log10(abs(imag(val))),nplt,nplt),'u')
% tools.pplot(X,Y,reshape(abs(val)-abs(real(val)),nplt,nplt),'u')
% hold on
% plot(sys.chnkr,'.')
% hold off
% clim([-.3,.3])
% clim([-.01,.01])

figure(2)
tools.pplot(X,Y,reshape(abs(val),nplt,nplt),'u')
% hold on
% plot(sys.pts.r(1,:),sys.pts.r(2,:),'k.')
% plot(sys.chnkr_o.r(1,:),sys.chnkr_o.r(2,:),'k.')
% hold off

figure(3)
tools.pplot(X,Y,reshape(abs(val-val2),nplt,nplt),'u')