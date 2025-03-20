function vals = d_two_flat_kern_lump(ssrc,dsrc,targ,sys)


cosmap = @(s) struct("r",[cos(s.r(1,:));s.r(2,:)],"n",s.n);
cosmap2 = @(s) struct("r",[cos(s.r(1,:));s.r(2,:)]);
acosmap = @(s) struct("r",[acos(s.r(1,:));s.r(2,:)],"n",s.n);

u_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'s');

ud_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');

rs_val = @(s,t) chnk.helm2d.kern(sys.zk,s,cosmap(t),'s');
rd_val = @(s,t) chnk.helm2d.kern(sys.zk,s,cosmap(t),'d');

g_val = @(s,t) chnk.helm2d.kern(sys.zk,cosmap(s),t,'s');
g_val_a = @(s,t) chnk.helm2d.kern(sys.zk,cosmap(s),cosmap2(t),'s');


im = abs(targ.r(1,:))<=1;
io = ~im;

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
rhs = zeros(sys.chnkr.npt,1);

% sum(im)
% tic;

for i = 1:nssrc
    srci = [];
    srci.r = ssrc.r(:,i);
    srci.n = ssrc.n(:,i);

    rhs = rhs + (rs_val(srci,sys.chnkr)).*ssrc.charge(i);
    vals = vals + (u_val(srci,targ)).*ssrc.charge(i);

end

for i = 1:ndsrc
    srci = [];
    srci.r = dsrc.r(:,i);
    srci.n = dsrc.n(:,i);

    rhs = rhs + (rd_val(srci,sys.chnkr)).*dsrc.charge(i);
    vals = vals + (ud_val(srci,targ)).*dsrc.charge(i);
  
end

sigma = sys.sysinv(rhs);
opts = []; opts.eps = 1e-8;
vals_scatm = chunkerkerneval(sys.chnkr,g_val_a,sigma,targm.r,opts);
vals_scato = chunkerkerneval(sys.chnkr,g_val,sigma,targo.r,opts);

vals(im) = vals(im) - vals_scatm;
vals(io) = vals(io) - vals_scato;
% toc
end

