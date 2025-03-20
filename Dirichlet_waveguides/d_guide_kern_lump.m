function vals = d_guide_kern_lump(ssrc,dsrc,targ,sys)


u_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'s');
ud_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');

% u_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'sprime');
% ud_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'dprime');

rs_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'s');
rd_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');

g_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');
% g_val = @(s,t) helm_graddx(sys.zk,s, t);

h = sys.d/2;

iu = targ.r(2,:)>h;
id = targ.r(2,:)<-h;
im = ~iu & ~id;

targu = [];
targu.r = targ.r(:,iu);
targu.n = targ.n(:,iu);

targm = [];
targm.r = targ.r(:,im);
targm.n = targ.n(:,im);

targd = [];
targd.r = targ.r(:,id);
targd.n = targ.n(:,id);

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
    if ssrc.r(2,i)>h
        srci_r = [];
        srci_r.r = [srci.r(1);h - (srci.r(2)-h)];
        srci_r.n = [srci.n(1);-srci.n(2)];
        vals(iu) = vals(iu)+(u_val(srci,targu) - u_val(srci_r,targu)).*ssrc.charge(i);

    elseif ssrc.r(2,i)<-h
        srci_r = [];
        srci_r.r = [srci.r(1);-h - (srci.r(2)+h)];
        srci_r.n = [srci.n(1);-srci.n(2)];

        vals(id) = vals(id)+(u_val(srci,targd) - u_val(srci_r,targd)).*ssrc.charge(i);
    else
        r_signs = [-1,-1,1,1];
        r_y1 = [h - (srci.r(2)-h),-h - (srci.r(2)+h)];
        r_y = [r_y1,h- (r_y1(2)-h),-h - (r_y1(1)+h)];

        srci_r = [];
        srci_r.r = [srci.r(1)+0*r_y;r_y];
        srci_r.n = [srci.n(1).*abs(r_signs);r_signs.*srci.n(2)]; 

        rhs = rhs + (rs_val(srci,sys.chnkr)+rs_val(srci_r,sys.chnkr)*r_signs.').*ssrc.charge(i);
        vals(im) = vals(im) + (u_val(srci,targm)+ u_val(srci_r,targm)*r_signs.').*ssrc.charge(i);
    end
end

for i = 1:ndsrc
    srci = [];
    srci.r = dsrc.r(:,i);
    srci.n = dsrc.n(:,i);
    if dsrc.r(2,i)>h
        srci_r = [];
        srci_r.r = [srci.r(1);h - (srci.r(2)-h)];
        srci_r.n = [srci.n(1);-srci.n(2)];
        vals(iu) = vals(iu)+(ud_val(srci,targu) - ud_val(srci_r,targu)).*dsrc.charge(i);

    elseif dsrc.r(2,i)<-h
        srci_r = [];
        srci_r.r = [srci.r(1);-h - (srci.r(2)+h)];
        srci_r.n = [srci.n(1);-srci.n(2)];

        vals(id) = vals(id)+(ud_val(srci,targd) - ud_val(srci_r,targd)).*dsrc.charge(i);
    else
        r_signs = [-1,-1,1,1];
        r_y1 = [h - (srci.r(2)-h),-h - (srci.r(2)+h)];
        r_y = [r_y1,h- (r_y1(2)-h),-h - (r_y1(1)+h)];

        srci_r = [];
        srci_r.r = [srci.r(1)+0*r_y;r_y];
        srci_r.n = [srci.n(1).*abs(r_signs);r_signs.*srci.n(2)]; 

        rhs = rhs + (rd_val(srci,sys.chnkr)+rd_val(srci_r,sys.chnkr)*r_signs.').*dsrc.charge(i);
        vals(im) = vals(im) + (ud_val(srci,targm)+ ud_val(srci_r,targm)*r_signs.').*dsrc.charge(i);
    end
end

sigma = sys.sysinv(rhs);
opts = []; opts.eps = 1e-8;
vals_scat = chunkerkerneval(sys.chnkr,g_val,sigma,targm.r,opts);

vals(im) = vals(im) - vals_scat;

% toc
end

