%DEMO_SCATTER_SLIT
%
% Define a scattering problem on a slit-like domain and solve
%

iseed = 8675309;
rng(iseed);
addpaths_loc();

% incident wave and other problem definitions

src0 = [0.7;-5];
zk = 2.3 + 1i*0.001;
strengths = 1.0;

% define geometry

% chunkpoly is based on vertices

% 4 wide, 1 tall rectangle at origin

verts = [ [-2;-0.5],[2;-0.5],[2;0.5],[-2;0.5] ];

cparams = [];
cparams.eps = 1.0e-5;
pref = []; 
pref.k = 16;
start = tic; chnkr = chunkpoly(verts,cparams,pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

assert(checkadjinfo(chnkr) == 0);
chnkr = chnkr.refine(); chnkr = chnkr.sort();

% make 2 shifted copies of this chnkr and merge them

chnkr1 = chnkr;
chnkr1.r = chnkr1.r - [3;0]; % left rectangle
chnkr2 = chnkr;
chnkr2.r = chnkr2.r + [3;0]; % right rectangle


% OK, this is clunky. For now, using an array of chnkrs helps with
% plotting reasonably... the merged chunker doesn't have the data about
% which piece of the geometry each chunk belongs to.

chnkrs = [chnkr1,chnkr2];
chnkr = chunkermerge(chnkrs);

% plot geometry and data

figure(1)
clf
plot(chnkrs,'-b')
hold on
quiver(chnkrs,'r')
axis equal

%


% solve and visualize the solution

% build CFIE

fkern = @(s,t,stau,ttau) chnk.helm2d.kern(zk,s,t,stau,ttau,'C',1);
opdims(1) = 1; opdims(2) = 1;

opts = [];
opts.quadorder = 30;
start = tic; sysmat = chunkmat(chnkr,fkern,opts);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = 0.5*eye(chnkr.k*chnkr.nch) + sysmat;

% get the boundary data for a source located at the point above

kerns = @(s,t,sn,tn) chnk.helm2d.kern(zk,s,t,sn,tn,'s');
targs = chnkr.r; targs = reshape(targs,2,chnkr.k*chnkr.nch);
targstau = taus(chnkr); 
targstau = reshape(targstau,2,chnkr.k*chnkr.nch);

kernmats = kerns(src0,targs,[],targstau);
ubdry = -kernmats*strengths;

rhs = ubdry; rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

% evaluate at targets and plot

rmin = min(chnkr); rmax = max(chnkr);
xl = rmax(1)-rmin(1);
yl = rmax(2)-rmin(2);
nplot = 300;
xtarg = linspace(-6,6,nplot);
ytarg = linspace(-6,6,nplot);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);

start = tic;
chnkr2 = chnkr;
chnkr2 = chnkr2.makedatarows(1);
chnkr2.data(1,:) = sol(:);
optref = []; optref.nover = 4;
chnkr2 = chnkr2.refine(optref);
sol2 = chnkr2.data(1,:);
t1 = toc(start);

fprintf('%5.2e s : time to oversample boundary\n',t1)

%

start = tic; in = chunkerinflam(chnkr,targets); t1 = toc(start);
out = ~in;

fprintf('%5.2e s : time to find points in domain\n',t1)

% compute layer potential based on oversample boundary

wts2 = whts(chnkr2);

matfun = @(i,j) kernbyindexr(i,j,targets(:,out),chnkr2,wts2,fkern,opdims);
[pr,ptau,pw,pin] = proxy_square_pts();

pxyfun = @(rc,rx,cx,slf,nbr,l,ctr) proxyfunr(rc,rx,slf,nbr,l,ctr,chnkr2,wts2, ...
    fkern,opdims,pr,ptau,pw,pin);

xflam = chnkr2.r(:,:);

start = tic; F = ifmm(matfun,targets(:,out),xflam,200,1e-14,pxyfun); 
t1 = toc(start);
fprintf('%5.2e s : time for ifmm form (for plotting)\n',t1)
start = tic;
uscat = ifmm_mv(F,sol2(:),matfun); t1 = toc(start);
fprintf('%5.2e s : time for ifmm apply (for plotting)\n',t1)

uin = kerns(src0,targets(:,out),[],[])*strengths;
utot = uscat(:)+uin(:);

%%

maxin = max(abs(uin(:)));
maxsc = max(abs(uscat(:)));
maxtot = max(abs(utot(:)));

maxu = max(max(maxin,maxsc),maxtot);

%

figure(2)
clf
subplot(1,3,1)
zztarg = nan(size(xxtarg));
zztarg(out) = uin;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
hold on
plot(chnkrs,'g')
axis equal
axis tight
colormap(redblue)
caxis([-maxu,maxu])
title('$u_{in}$','Interpreter','latex','FontSize',24)


subplot(1,3,2)
zztarg = nan(size(xxtarg));
zztarg(out) = uscat;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
hold on
plot(chnkrs,'g')
axis equal
axis tight
colormap(redblue)
caxis([-maxu,maxu])
title('$u_{scat}$','Interpreter','latex','FontSize',24)

subplot(1,3,3)
zztarg = nan(size(xxtarg));
zztarg(out) = utot;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
hold on
plot(chnkrs,'g')
axis equal
axis tight
colormap(redblue)
caxis([-maxu,maxu])
title('$u_{tot}$','Interpreter','latex','FontSize',24)

