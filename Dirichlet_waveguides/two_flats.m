zk = 1;
test_acc = 0;
d = 4;
L = 40;

nover = 5;
sys = build_two_flats_sys(nover,d,zk);

Lplot = 3;
nplt = 100;

xts = linspace(-Lplot,Lplot,nplt);
yts = linspace(-Lplot,Lplot,nplt);
[X,Y] = meshgrid(xts,yts);
xtargs = X(:).'; ytargs = Y(:).';

targ = [];
targ.r = [xtargs;ytargs];
targ.n = [1;0]+0*targ.r;

src = [];
src.r = [0;d/4];
% src.r = [-1;10];
src.n = [1;0];
val = d_two_flat(src,targ,sys,'dp');
%%
figure(1)
tools.pplot(X,Y,reshape(imag(val),nplt,nplt),'u')
% tools.pplot(X,Y,reshape(log10(abs(imag(val))),nplt,nplt),'u')
% tools.pplot(X,Y,reshape(abs(val)-abs(real(val)),nplt,nplt),'u')
% hold on
% plot(sys.chnkr,'.')
% hold off
% clim([-.3,.3])
% clim([-.01,.01])

