function vals = d_guide_kern_cap(src,targ,sys,kop)


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

% [iouts, iiins] = wave_kern.flagout(sys,src.r);
% [ioutt, iiint] = wave_kern.flagout(sys,targ.r);

ioutt = abs(targ.r(2,:))>sys.d/2;
iiint = ~ioutt;

iouts = abs(src.r(2,:))>sys.d/2;
iiins = ~iouts;


targi = [];
targi.r = targ.r(:,iiint);
targi.n = targ.n(:,iiint);


targo = [];
targo.r = targ.r(:,ioutt);
targo.n = targ.n(:,ioutt);

nsrc = size(src.r,2);
ntarg = size(targ.r,2);

vals = zeros(ntarg,nsrc);

h = sys.d/2;

% sum(im)
% tic;
for i = 1:nsrc
    srci = [];
    srci.r = src.r(:,i);
    srci.n = src.n(:,i);
    

    r_signs = [-1,-1,1,1];
    r_y1 = [h - (srci.r(2)-h),-h - (srci.r(2)+h)];
    r_y = [r_y1,h- (r_y1(2)-h),-h - (r_y1(1)+h)];
    r_signs = r_signs.*(abs(r_y) < h);
    srci_r = [];
    srci_r.r = [srci.r(1)+0*r_y;r_y];
    srci_r.n = [srci.n(1).*abs(r_signs);r_signs.*srci.n(2)]; 

    rhs = r_val(srci,sys.chnkr);%+r_val(srci_r,sys.chnkr)*r_signs.';
    opts = []; opts.eps = 1e-8;
    if iouts(i)
    sigma = sys.sysinvo(rhs);
    vals_scat = chunkerkerneval(sys.chnkr,g_val,sigma,targo.r,opts);
    % valstot = u_val(srci,targo)+u_val(srci_r,targo)*r_signs.' - vals_scat;
    valstot = u_val(srci,targo) - vals_scat;
    else
    sigma = sys.sysinvi(rhs);
    vals_scat = chunkerkerneval(sys.chnkr,g_val,sigma,targi.r,opts);
    % valstot = u_val(srci,targi)+u_val(srci_r,targi)*r_signs.' - vals_scat;
    valstot = u_val(srci,targi)- vals_scat;
    end
    
    


    if iouts(i)
        vals(iiint,i) = 0;
        vals(ioutt,i) = valstot;
    else
        vals(iiint,i) = valstot;
        vals(ioutt,i) = 0;
    end
end
% toc
end
