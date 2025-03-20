function vals = d_obs_kern_lump(ssrc,dsrc,targ,sys)


cosmap = @(s) struct("r",[cos(s.r(1,:));s.r(2,:)],"n",s.n);
cosmap2 = @(s) struct("r",[cos(s.r(1,:));s.r(2,:)]);
acosmap = @(s) struct("r",[acos(s.r(1,:));s.r(2,:)],"n",s.n);

u_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'s');

ud_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');

rs_val = @(s,t) chnk.helm2d.kern(sys.zk,s,cosmap(t),'s');
rd_val = @(s,t) chnk.helm2d.kern(sys.zk,s,cosmap(t),'d');

rso_val = @(s,t) leaky.kern_data(sys.sys_o, s.r(:));
rdo_val = @(s,t) leaky.kern_data_dp(sys.sys_o, s.r(:),s.n(:));

g_val = @(s,t) chnk.helm2d.kern(sys.zk,cosmap(s),t,'s');
g_val_a = @(s,t) chnk.helm2d.kern(sys.zk,cosmap(s),cosmap2(t),'s');


[iout, iiin] = wave_kern.flagout(sys.sys_o,targ.r);
iout = iout(:).';
iiin = iiin(:).';
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


nssrc = size(ssrc.r,2);
ndsrc = size(dsrc.r,2);
ntarg = size(targ.r,2);

vals = zeros(ntarg,1);
rhs_f = zeros(sys.chnkr.npt,1);
rhs_o = zeros(2*sys.chnkr_o.npt,1);

% sum(im)
% tic;

for i = 1:nssrc
    srci = [];
    srci.r = ssrc.r(:,i);
    srci.n = ssrc.n(:,i);

    rhs_f = rhs_f + (rs_val(srci,sys.chnkr)).*ssrc.charge(i);
    rhs_o = rhs_o + rso_val(srci,sys.chnkr_o).*ssrc.charge(i);
    % vals = vals + (u_val(srci,targ)).*ssrc.charge(i);
    opts.isp = 2;
    [uvals, ~] = ...
        leaky.sys_apply(sys.sys_o, NaN, targout, targin, ...
            iout, iiin, srci.r,srci.n,1,0,opts);
    vals = vals + uvals.*ssrc.charge(i);

end

for i = 1:ndsrc
    srci = [];
    srci.r = dsrc.r(:,i);
    srci.n = dsrc.n(:,i);

    rhs_f = rhs_f + (rd_val(srci,sys.chnkr)).*dsrc.charge(i);
    rhs_o = rhs_o + rdo_val(srci,sys.chnkr_o).*dsrc.charge(i);
    % vals = vals + (ud_val(srci,targ)).*dsrc.charge(i);
    opts.isp = 2;
    [uvals, ~] = ...
        leaky.sys_apply(sys.sys_o, NaN, targout, targin, ...
            iout, iiin, srci.r,srci.n,1,1,opts);
    vals = vals + uvals.*dsrc.charge(i);
  
end
rhs = [rhs_f;rhs_o];
dens = sys.sysinv(rhs);
sigma = dens(1:sys.chnkr.npt);
soln_o = dens((sys.chnkr.npt+1):end);

opts = []; opts.eps = 1e-8;
vals_scatm = chunkerkerneval(sys.chnkr,g_val_a,sigma,targm.r,opts);
vals_scato = chunkerkerneval(sys.chnkr,g_val,sigma,targo.r,opts);

vals(im) = vals(im) - vals_scatm;
vals(io) = vals(io) - vals_scato;


opts.isp = 1;
[uvals, ~] = ...
    leaky.sys_apply(sys.sys_o, soln_o, targout, targin, ...
        iout, iiin, srci.r,srci.n,1,0,opts);
vals(iout) = vals(iout) + uvals(iout);
vals(iiin) = vals(iiin) + uvals(iiin);
% toc
end

