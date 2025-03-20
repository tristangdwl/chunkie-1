zk = 1;
test_acc = 0;
d = [8,4];
L = 40;

nover = 3;
sys = build_four_flats_sys(nover,d,zk);

Lplot = 5;
nplt = 200;

xts = linspace(-Lplot,Lplot,nplt);
yts = linspace(-Lplot,Lplot,nplt);
[X,Y] = meshgrid(xts,yts);
xtargs = X(:).'; ytargs = Y(:).';

targ = [];
targ.r = [xtargs;ytargs];
targ.n = [1;0]+0*targ.r;

src = [];
src.r = [-1;d(1)/4-.01];
% src.r = [-2;];
% src.r = [-1;10];
src.n = [1;0];
val = d_four_flat(src,targ,sys,'s');
%%
figure(1)
tools.pplot(X,Y,reshape(real(val),nplt,nplt),'real(u)')
hold on
scatter(sys.pts.r(1,:),sys.pts.r(2,:),'r.')
hold off
% tools.pplot(X,Y,reshape(log10(abs(imag(val))),nplt,nplt),'u')
% tools.pplot(X,Y,reshape(abs(val)-abs(real(val)),nplt,nplt),'u')
% hold on
% plot(sys.chnkr,'.')
% hold off
% clim([-.3,.3])
% clim([-.01,.01])

