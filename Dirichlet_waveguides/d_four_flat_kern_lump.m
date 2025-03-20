function vals = d_four_flat_kern_lump(ssrc,dsrc,targ,sys)


cosmap = @(s) struct("r",[cos(s.r(1,:));s.r(2,:)],"n",s.n);
cosmap2 = @(s) struct("r",[cos(s.r(1,:));s.r(2,:)]);
acosmap = @(s) struct("r",[acos(s.r(1,:));s.r(2,:)],"n",s.n);

u_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'s');
ud_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');

% u_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'sprime');
% ud_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'dprime');

rs_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'s');
rd_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');

g_val = @(s,t) chnk.helm2d.kern(sys.zk,cosmap(s),t,'s');
g_val_a = @(s,t) chnk.helm2d.kern(sys.zk,cosmap(s),cosmap2(t),'s');

% g_val = @(s,t) helm_gradsx(sys.zk, cosmap(s), t);
% g_val_a = @(s,t) helm_gradsx(sys.zk, cosmap(s), cosmap2(t));

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

    rhs = rhs + (rs_val(srci,sys.pts)).*ssrc.charge(i);
    vals = vals + (u_val(srci,targ)).*ssrc.charge(i);

end

for i = 1:ndsrc
    srci = [];
    srci.r = dsrc.r(:,i);
    srci.n = dsrc.n(:,i);

    rhs = rhs + (rd_val(srci,sys.pts)).*dsrc.charge(i);
    vals = vals + (ud_val(srci,targ)).*dsrc.charge(i);
  
end

sigma = sys.sysinv(rhs);

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

vals = vals - vals_scatl-vals_scatr;
% toc
end

