zk = 1;
test_acc = 0;
d = 4;
L = 40;

nover = 1;

% sys = build_dirichlet_sys(nover,d,L,zk);
sys = build_dirichlet_sys_cap(nover,d,0,zk);

Lplot = 30;
nplt = 200;

xts = linspace(-Lplot,Lplot,nplt);
yts = linspace(-Lplot,Lplot,nplt);
[X,Y] = meshgrid(xts,yts);
xtargs = X(:).'; ytargs = Y(:).';

targ = [];
targ.r = [xtargs;ytargs];
targ.n = [1;0]+0*targ.r;

src = [];
src.r = [-1;d/4];
% src.r = [-1;10];
src.n = [1;0];
% val = d_guide_kern(src,targ,sys,'s');
val = d_guide_kern_cap(src,targ,sys,'s');
%%
figure(1)
tools.pplot(X,Y,reshape(real(val),nplt,nplt),'u')
tools.pplot(X,Y,reshape(log10(abs(imag(val))),nplt,nplt),'u')
% tools.pplot(X,Y,reshape(abs(val)-abs(real(val)),nplt,nplt),'u')
% hold on
% plot(sys.chnkr,'.')
% hold off
% clim([-.3,.3])
% clim([-.01,.01])

if test_acc
src.r = [0;20];
[val,val_true] = d_guide_kern_err(src,targ,sys,'dp');
figure(2)
subplot(1,2,1)
tools.pplot(X,Y,reshape(imag(val),nplt,nplt),'u')

subplot(1,2,2)
tools.pplot(X,Y,reshape(log10(abs(val_true-val)),nplt,nplt),'u')
end



