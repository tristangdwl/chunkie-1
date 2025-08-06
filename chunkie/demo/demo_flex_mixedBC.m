%DEMO_FREE_PLATE_SCATTER
%
% Define an interior scattering problem on a starfish-shaped domain and 
% solve
%

zk = 10;

% discretize domain

cparams = [];
cparams.eps = 1.0e-8;
cparams.nover = 0;
cparams.maxchunklen = 4.0/zk; % setting a chunk length helps when the
                              % frequency is known
                              
pref = []; 
pref.k = 16;
narms1 = 4;
narms2 = 3;
amp = 0.25*0.5*0;
start = tic; 
% chnkr_f = chunkerfunc(@(t) starfish(t,narms1,amp,[0;4],[],3),cparams,pref); 
chnkr_f = chunkerfunc(@(t) starfish(t,narms1,amp,[0;0],[],5),cparams,pref); 
chnkr_f = chunkerfunc(@(t) starfish(t,narms1,amp,[0;0],[],1),cparams,pref); 
% chnkr_f = chunkerfunc(@(t) starfish(t,narms1,amp,[0;0],[],1),cparams,pref); 
% chnkr_f = chunkerfunc(@(t) starfish(t,narms1,amp,[0;0],[],1.7),cparams,pref); 
% chnkr_c = chunkerfunc(@(t) starfish(t,narms2,amp,[],[],1),cparams,pref);% chnkr_c = refine(chnkr_c,struct('nover',1));
% chnkr_c = chunkerfunc(@(t) starfish(t,narms2,amp,[],[],1),cparams,pref); 

chnkr = merge([reverse(chnkr_f),chnkr_c]);
chnkr = merge([(chnkr_f),chnkr_c]);
chnkr = chnkr_c;chnkr = chnkr_f;
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

% plot geometry and data

figure(1)
clf
plot(chnkr,'-')
hold on
quiver(chnkr)
% quiver(chnkr_f)
% quiver(chnkr_c)
axis equal

%% Get system matrix and right hand side
ibc = 0;

start = tic;
if ibc == 2
sys = mix_sysmat(chnkr_f,chnkr_c,zk,1,-1);
elseif ibc == 1
sys = free_sysmat(chnkr,zk,1);
else
sys = clamped_sysmat(chnkr,zk,1); %sgn = -1 for exterior
end
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

% building RHS
src =[]; src.r = [1;2];
src =[]; src.r = [0.;0];src.r = [0;8.];
src.n = [1;0];
free_bc_kern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'free_plate_bcs');
clamp_bc_kern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_plate_bcs');

flex_free_kern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 's');

% free_bc_kern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_plate_bcs');
% clamp_bc_kern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_plate_bcs');

% flex_free_kern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'santi');

% free_bc_kern = @(s,t) kernel('z',2,1).eval(s,t);
% clamp_bc_kern = @(s,t) chnk.helm2d.kern(zk, s, t, 'c2trans',[0,1]);
% 
% flex_free_kern = @(s,t) chnk.helm2d.kern(zk, s, t, 's');

% get free plate BCs
rhs_f = -free_bc_kern(src,chnkr_f);
% get clamped plate BCs
rhs_c = -clamp_bc_kern(src,chnkr_c);
rhs = [rhs_f;rhs_c];

if ibc == 1
rhs = -free_bc_kern(src,chnkr);
elseif ibc == 0
rhs = -clamp_bc_kern(src,chnkr);
end

%% Solve and plot
% Solving linear system

% start = tic; sol = gmres(sys,rhs,[],1e-12,500); t1 = toc(start);
% fprintf('%5.2e s : time for dense gmres\n',t1)    
sol = sys\rhs;

% evaluate at targets and plot


rmin = min(chnkr); rmax = max(chnkr);
xl = rmax(1)-rmin(1);
yl = rmax(2)-rmin(2);
nplot = 100;% nplot = 150;
xtarg = linspace(rmin(1)-xl/10,rmax(1)+xl/10,nplot); 
ytarg = linspace(rmin(2)-yl/10,rmax(2)+yl/10,nplot);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);

start = tic; 
in = chunkerinterior(chnkr_c,{xtarg,ytarg}); out_c = ~in;
in = chunkerinterior(chnkr_f,{xtarg,ytarg}); out_f = ~in;
out = out_c & ~out_f;

in = chunkerinterior(chnkr,{xtarg,ytarg}); out = ~in; out=in;

t1 = toc(start);

fprintf('%5.2e s : time to find points in domain\n',t1)

start1 = tic;
if ibc == 2
uscat = mixed_eval(chnkr_f,chnkr_c,sol,targets(:,out),zk);
elseif ibc == 1
uscat = free_eval(chnkr,sol,targets(:,out),zk);
else
uscat = clamp_eval(chnkr,sol,targets(:,out),zk);
end
% uscat = free_eval(chnkr,sol,targets(:,out),zk);
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

%%
uin = flex_free_kern(src,struct("r",targets(:,out)));
utot = uscat(:)+uin(:);

maxu = max(abs(uin(:)));
maxu = max(abs(uscat(:)));

figure(2)
clf

t = tiledlayout(1,2,'TileSpacing','compact');

nexttile
zztarg = nan(size(xxtarg))+NaN*1i;
zztarg(out) = uin;
h=pcolor(xxtarg,yytarg,imag(zztarg)); h.FaceColor="interp";
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{\textrm{inc}}$','Interpreter','latex','FontSize',12)

nexttile
zztarg = nan(size(xxtarg))+NaN*1i;
zztarg(out) = uscat;
h=pcolor(xxtarg,yytarg,imag(zztarg)); h.FaceColor="interp";
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{\textrm{scat}}$','Interpreter','latex','FontSize',12)

colorbar

figure(4)
clf
t = tiledlayout(1,3,'TileSpacing','compact');

nexttile
zztarg = nan(size(xxtarg))+NaN*1i;
zztarg(out) = uin;
h=pcolor(xxtarg,yytarg,real(zztarg)); h.FaceColor="interp";
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{\textrm{inc}}$','Interpreter','latex','FontSize',12)

nexttile
zztarg = nan(size(xxtarg))+NaN*1i;
zztarg(out) = uscat;
h=pcolor(xxtarg,yytarg,real(zztarg)); h.FaceColor="interp";
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{\textrm{scat}}$','Interpreter','latex','FontSize',12)

nexttile
zztarg = nan(size(xxtarg))+NaN*1i;
zztarg(out) = utot;
h=pcolor(xxtarg,yytarg,real(zztarg)); h.FaceColor="interp";
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)

colorbar

figure(3);clf;%subplot(2,1,1)
zztarg = nan(size(xxtarg))+NaN*1i;
zztarg(out) = utot;
h=pcolor(xxtarg,yytarg,log10(abs(zztarg))); h.FaceColor="interp";
set(h,'EdgeColor','none')
% clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$\log_{10}$ error','Interpreter','latex','FontSize',12)
if ibc == 2
title(t,"Mixed BCs")
elseif ibc == 1
title(t,"Free BCs")
else
title(t,"Clamped BCs")
end

colorbar


figure(5);clf
zztarg = nan(size(xxtarg))+NaN*1i;
zztarg(out) = utot;
h=pcolor(xxtarg,yytarg,real(zztarg)); h.FaceColor="interp";
set(h,'EdgeColor','none')
clim([-maxu,maxu])
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)

colorbar

return
%%
tic;
h = 1e-3;
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:)+h; targets(2,:) = yytarg(:);
% uscat1 = mixed_eval(chnkr_f,chnkr_c,sol,targets(:,out),zk);
uscat1 = free_eval(chnkr,sol,targets(:,out),zk);
uin1 = flex_free_kern(src,struct("r",targets(:,out)));
utot1 = uscat1(:)+uin1(:);

targets(1,:) = xxtarg(:)-h; targets(2,:) = yytarg(:);
% uscat2 = mixed_eval(chnkr_f,chnkr_c,sol,targets(:,out),zk);
uscat2 = free_eval(chnkr,sol,targets(:,out),zk);
uin2 = flex_free_kern(src,struct("r",targets(:,out)));
utot2 = uscat2(:)+uin2(:);

targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:)+h;
% uscat3 = mixed_eval(chnkr_f,chnkr_c,sol,targets(:,out),zk);
uscat3 = free_eval(chnkr,sol,targets(:,out),zk);
uin3 = flex_free_kern(src,struct("r",targets(:,out)));
utot3 = uscat3(:)+uin3(:);

targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:)-h;
% uscat4 = mixed_eval(chnkr_f,chnkr_c,sol,targets(:,out),zk);
uscat4 = free_eval(chnkr,sol,targets(:,out),zk);
uin4 = flex_free_kern(src,struct("r",targets(:,out)));
utot4 = uscat4(:)+uin4(:);

uh = (utot1+utot2+utot3+utot4-4*utot)/h^2 + zk^2 *utot;
% uh = (uin1+uin2+uin3+uin4-4*uin)/h^2 + zk^2 *uin;
toc;
%%
figure(6);clf
zztarg = nan(size(xxtarg))+NaN*1i;
zztarg(out) = uh;
h=pcolor(xxtarg,yytarg,real(zztarg)); h.FaceColor="interp";
set(h,'EdgeColor','none')
% clim([-maxu,maxu]*5)
colormap(redblue);
hold on
plot(chnkr,'k','LineWidth',2)
axis equal tight
set(gca, "box","off","Xtick",[],"Ytick",[]);
title('$(\Delta+k^2) u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)

colorbar
% %%
% sol_f = sol(1:2*chnkr_f.npt);
% sys_fc = free2clamp_sysmat(chnkr_f,chnkr_c,zk);
% u2 = free_eval(chnkr_f,sol_f,chnkr_c,zk);
% norm(u2-sys_fc(1:2:end,:)*sol_f)

% %%
% ikern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'free_plate_eval');
% % ikern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_plate_eval');
% % ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');
% ikern = @(s,t) chnk.helm2d.kern(zk, s, t, 's');
% tic;
% amat = chunkerkernevalmat(chnkr, ikern,targets(:,out));
% toc;
% s = svd(amat);
% 
% figure(7);clf
% plot(log10(s))
% % %%
% 
% [U,S,V] = svd(amat);
% %%
% figure(6);clf
% zztarg = nan(size(xxtarg))+NaN*1i;
% zztarg(out) = U(:,1);
% h=pcolor(xxtarg,yytarg,real(zztarg)); h.FaceColor="interp";
% set(h,'EdgeColor','none')
% % clim([-maxu,maxu]*5)
% colormap(redblue);
% hold on
% plot(chnkr,'k','LineWidth',2)
% axis equal tight
% set(gca, "box","off","Xtick",[],"Ytick",[]);
% title('$(\Delta+k^2) u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)
% 
% colorbar


function sys = mix_sysmat(chnkr_f,chnkr_c,zk,sgn_f,sgn_c)
% build the block system matrix
    sys_ff = free_sysmat(chnkr_f,zk,sgn_f);
    sys_fc = free2clamp_sysmat(chnkr_f,chnkr_c,zk);
    sys_cf = clamp2free_sysmat(chnkr_f,chnkr_c,zk);
    sys_cc = clamped_sysmat(chnkr_c,zk,sgn_c);

    sys = [sys_ff,sys_cf;sys_fc,sys_cc];
    % sys = [sys_ff,sys_fc;sys_cf,sys_cc];
end


function sys = free_sysmat(chnkr,zk,sgn)
    % assembling system matrix for free boundary

    sys = clamped_sysmat(chnkr,zk,sgn);
end

function sys = free2clamp_sysmat(chnkr_f,chnkr_c,zk)
    % assemble matrix for free boundary talking to clamped boundary
    
    % defining free plate kernels
    fkern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'free_to_clamped'); 
    % building system matrix
    
    sys = chunkerkernevalmat(chnkr_f,fkern,chnkr_c);
end

function sys = clamp2free_sysmat(chnkr_f,chnkr_c,zk)
    % assemble matrix for clamped boundary talking to free boundary

    % defining clamped plate kernels
    fkern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_to_free'); 
    
    % building system matrix
    
    sys = chunkerkernevalmat(chnkr_c,fkern,chnkr_f);
end

function sys = clamped_sysmat(chnkr,zk,sgn)
    % assembling system matrix for clamped boundary
    
    fkern =  @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_plate');
    
    kappa = signed_curvature(chnkr);
    kappa = kappa(:);
    
    opts = [];
    opts.sing = 'log';
    
    sys = chunkermat(chnkr,fkern, opts);
    sys = sys - sgn*0.5*eye(2*chnkr.npt);
    sys(2:2:end,1:2:end) = sys(2:2:end,1:2:end) + sgn*kappa.*eye(chnkr.npt);
end

function uscat = mixed_eval(chnkr_f,chnkr_c,sol,targets,zk)
    % evaluate layer potentials one BC at a time
    sol_f = sol(1:2*chnkr_f.npt);
    sol_c = sol((1+2*chnkr_f.npt):end);

    uscat = free_eval(chnkr_f,sol_f,targets,zk);
    uscat = uscat+clamp_eval(chnkr_c,sol_c,targets,zk);

end

function uscat = free_eval(chnkr,sol,targets,zk)
    % evaluate layer potentials from free plate representation
    ikern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'free_plate_eval'); 
    
    uscat = chunkerkerneval(chnkr, ikern,sol,targets);
end

function uscat = clamp_eval(chnkr,sol,targets,zk)
    % evaluate layer potentials from clamped plate representation
    ikern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_plate_eval');   
    
    uscat = chunkerkerneval(chnkr, ikern,sol,targets);
end