function vals = d_two_flat(src,targ,sys,kop)

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
    r_val = @(s,t) chnk.helm2d.kern(sys.zk,s,cosmap(t),'d');
else
    r_val = @(s,t) chnk.helm2d.kern(sys.zk,s,cosmap(t),'s');
end

if igrad
    g_val = @(s,t) helm_gradsx(sys.zk, cosmap(s), t);
    g_val_a = @(s,t) helm_gradsx(sys.zk, cosmap(s), cosmap2(t));
else
    g_val = @(s,t) chnk.helm2d.kern(sys.zk,cosmap(s),t,'s');
    g_val_a = @(s,t) chnk.helm2d.kern(sys.zk,cosmap(s),cosmap2(t),'s');
end

nsrc = size(src.r(:,:),2);
ntarg = size(targ.r(:,:),2);

vals = zeros(ntarg,nsrc);

im = abs(targ.r(1,:))<=1;
io = ~im;

targm = [];
targm.r = targ.r(:,im);
targm.n = targ.n(:,im);

targm = acosmap(targm);

targo = [];
targo.r = targ.r(:,io);
targo.n = targ.n(:,io);

% tic;
for i = 1:nsrc
    srci = [];
    srci.r = src.r(:,i);
    srci.n = src.n(:,i);
 
    rhs = r_val(srci,sys.chnkr);
    sigma = sys.sysinv(rhs);
    opts = []; opts.eps = 1e-8; opts.forcesmooth = false;
    
    vals_scat_m = chunkerkerneval(sys.chnkr,g_val_a,sigma,targm.r,opts);

    vals_scat_o = chunkerkerneval(sys.chnkr,g_val,sigma,targo.r,opts);

    vals_scat = zeros(ntarg,1);
    vals_scat(im) = vals_scat_m;
    vals_scat(io) = vals_scat_o;
    % vals_scat = chunkerkerneval(sys.chnkr,g_val,sigma,targ.r,opts);
    vals(:,i) = u_val(srci,targ) - vals_scat;
end

dx = targ.r(1,:).' - src.r(1,:);
dy = targ.r(2,:).' - src.r(2,:);
r = sqrt(dx.^2 + dy.^2);
vals(r<1e-14) = 0;

% toc
end

