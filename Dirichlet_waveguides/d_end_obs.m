function vals = d_end_obs(src,targ,sys,kop)


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
    ro_val = @(s,t) leaky.kern_data_dp(sys.sys_o, s.r(:),s.n(:));
else
    r_val = @(s,t) chnk.helm2d.kern(sys.zk,s,cosmap(t),'s');
    ro_val = @(s,t) leaky.kern_data(sys.sys_o, s.r(:));
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


[iout, iiin] = wave_kern.flagout(sys.sys_o,targ.r);
iout = iout(:).';
iiin = iiin(:).';

% iout = true(ntarg,1).';
% iiin = false(ntarg,1).';

targout = targ.r(:,iout);
targin  = targ.r(:,iiin);

im = abs(targ.r(1,:))<=1 & iout;
io = ~im & iout;

targm = [];
targm.r = targ.r(:,im);
targm.n = targ.n(:,im);

targm = acosmap(targm);

targo = [];
targo.r = targ.r(:,io);
targo.n = targ.n(:,io);

if ntarg == 0 || nsrc == 0, return, end




% tic;
for i = 1:nsrc
    srci = [];
    srci.r = src.r(:,i);
    srci.n = src.n(:,i);
 
    rhs_f = r_val(srci,sys.chnkr);
    rhs_o = -ro_val(srci,sys.chnkr_o);
    rhs = [rhs_f;rhs_o];
    dens = sys.sysinv(rhs);
    sigma = dens(1:sys.chnkr.npt);
    soln_o = dens((sys.chnkr.npt+1):end);

    opts = []; opts.eps = 1e-8; opts.forcesmooth = false;
    
    vals_scat_m = chunkerkerneval(sys.chnkr,g_val_a,sigma,targm.r,opts);
    vals_scat_o = chunkerkerneval(sys.chnkr,g_val,sigma,targo.r,opts);

    opts.isp = 0;
    [uvals, tmp] = ...
        leaky.sys_apply(sys.sys_o, soln_o, targout, targin, ...
            iout, iiin, srci.r,srci.n,igrad+1,idp,opts);
    if igrad
        uvals = sum(tmp.*targ.n(:,:),1).';
    end

    vals_scat = zeros(ntarg,1);
    vals_scat(im) = vals_scat_m;
    vals_scat(io) = vals_scat_o;
    vals_scat(iout) = vals_scat(iout) - uvals(iout);
    vals_scat(iiin) = vals_scat(iiin) - uvals(iiin);
    % vals_scat = chunkerkerneval(sys.chnkr,g_val,sigma,targ.r,opts);
    % vals(:,i) = u_val(srci,targ) - vals_scat;
    vals(:,i) =  - vals_scat;
end
% toc
end

