function vals = d_four_flat(src,targ,sys,kop)


if strcmp(kop,'s')
    idp = 0;
    igrad = 0;
    u_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'s');
elseif strcmp(kop,'sp')
    idp = 0;
    igrad = 1;
    u_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'sprime');
elseif strcmp(kop,'d')
    idp = 1;
    igrad = 0;
    u_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');
elseif strcmp(kop,'dp')
    idp = 1;
    igrad = 1;
    u_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'dprime');
end


cosmap = @(s) struct("r",[cos(s.r(1,:));s.r(2,:)],"n",s.n);
cosmap2 = @(s) struct("r",[cos(s.r(1,:));s.r(2,:)]);
acosmap = @(s) struct("r",[acos(s.r(1,:));s.r(2,:)],"n",s.n);

if idp 
    r_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');
else
    r_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'s');
end

if igrad
    g_val = @(s,t) helm_gradsx(sys.zk, cosmap(s), t);
    g_val_a = @(s,t) helm_gradsx(sys.zk, cosmap(s), cosmap2(t));
else
    g_val = @(s,t) chnk.helm2d.kern(sys.zk,cosmap(s),t,'s');
    g_val_a = @(s,t) chnk.helm2d.kern(sys.zk,cosmap(s),cosmap2(t),'s');
end

nsrc = size(src.r,2);
ntarg = size(targ.r,2);

vals = zeros(ntarg,nsrc);

il = abs(targ.r(1,:)+1)<1;
ir = abs(targ.r(1,:)-1)<1;
ilo = ~il ;
iro = ~ ir;


targl = [];
targl.r = targ.r(:,il)+[1;0];
targl.n = targ.n(:,il);

targl = acosmap(targl);

targr = [];
targr.r = targ.r(:,ir)-[1;0];
targr.n = targ.n(:,ir);

targr = acosmap(targr);

targlo = [];
targlo.r = targ.r(:,ilo)+[1;0];
targlo.n = targ.n(:,ilo);

targro = [];
targro.r = targ.r(:,iro)-[1;0];
targro.n = targ.n(:,iro);

% tic;
for i = 1:nsrc
    srci = [];
    srci.r = src.r(:,i);
    srci.n = src.n(:,i);
 
    rhs = r_val(srci,sys.pts);
    sigma = sys.sysinv(rhs);

    w = sys.chnkr.weights;
    % vals(:,i) = u_val(srci,targ) - r_val(sys.pts,targ)*(sigma.*w(:));

    sigma_l = sigma(1:sys.chnkrl.npt);
    sigma_r = sigma((1+sys.chnkrl.npt):end);
    opts = []; opts.eps = 1e-8; opts.forcesmooth = false;

    vals_scat_l = chunkerkerneval(sys.chnkrl,g_val_a,sigma_l,targl.r,opts);
    vals_scat_lo = chunkerkerneval(sys.chnkrl,g_val,sigma_l,targlo.r,opts);

    vals_scat_r = chunkerkerneval(sys.chnkrr,g_val_a,sigma_r,targr.r,opts);
    vals_scat_ro = chunkerkerneval(sys.chnkrr,g_val,sigma_r,targro.r,opts);

    vals_scatl = zeros(ntarg,1);
    vals_scatr = zeros(ntarg,1);
    vals_scatl(il) = vals_scat_l;
    vals_scatr(ir) = vals_scat_r;
    vals_scatl(ilo) = vals_scat_lo;
    vals_scatr(iro) = vals_scat_ro;
    % vals_scat = chunkerkerneval(sys.chnkr,g_val,sigma,targ.r,opts);
    vals(:,i) = u_val(srci,targ) - vals_scatl -vals_scatr;
end
% toc
end

