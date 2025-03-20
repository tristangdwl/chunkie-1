
function [val,val_true] = d_guide_kern_err(src,targ,sys,kop)


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

if idp 
    r_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');
else
    r_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'s');
end

if igrad
    g_val = @(s,t) helm_graddx(sys.zk,s,t);
else
    g_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');
end


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

nsrc = size(src.r,2);
ntarg = size(targ.r,2);

val = zeros(ntarg,nsrc);
val_true = zeros(ntarg,nsrc);

for i = 1:nsrc
    srci = [];
    srci.r = src.r(:,i);
    srci.n = src.n(:,i);

    r_signs = 0*[-1,-1,1,1];
    r_y1 = [h - (srci.r(2)-h),-h - (srci.r(2)+h)];
    r_y = [r_y1,h- (r_y1(2)-h),-h - (r_y1(1)+h)];

    srci_r = [];
    srci_r.r = [srci.r(1)+0*r_y;r_y];
    srci_r.n = [srci.n(1).*abs(r_signs);r_signs.*srci.n(2)]; 

    rhs = r_val(srci,sys.chnkr)+r_val(srci_r,sys.chnkr)*r_signs.';

    sigma = sys.sysinv(rhs);
    opts = []; opts.eps = 1e-8; opts.forcesmooth = false;
    val(im,i) = chunkerkerneval(sys.chnkr,g_val,sigma,targm.r,opts)-u_val(srci_r,targm)*r_signs.';

    val_true(im,i) = u_val(srci,targm);

    val(iu,i) = 0;
    val(id,i) = 0;
end


end

