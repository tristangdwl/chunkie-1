function vals = d_guide_kern_lump_cap(ssrc,dsrc,targ,sys)


u_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'s');

ud_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');

rs_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'s');
rd_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');

g_val = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');

h = sys.d/2;



[ioutt, iiint] = wave_kern.flagout(sys,targ.r);

targi = [];
targi.r = targ.r(:,iiint);
targi.n = targ.n(:,iiint);

targo = [];
targo.r = targ.r(:,ioutt);
targo.n = targ.n(:,ioutt);

nssrc = size(ssrc.r,2);
ndsrc = size(dsrc.r,2);
ntarg = size(targ.r,2);

vals = zeros(ntarg,1);
rhsi = zeros(sys.chnkr.npt,1);
rhso = zeros(sys.chnkr.npt,1);

% sum(im)
% tic;
[iouts, iiins] = wave_kern.flagout(sys,ssrc.r);
for i = 1:nssrc
    srci = [];
    srci.r = ssrc.r(:,i);
    srci.n = ssrc.n(:,i);

    % r_signs = [-1,-1,1,1];
    % r_y1 = [h - (srci.r(2)-h),-h - (srci.r(2)+h)];
    % r_y = [r_y1,h- (r_y1(2)-h),-h - (r_y1(1)+h)];
    % r_signs = r_signs.*(abs(r_y) < h);
    % 
    % srci_r = [];
    % srci_r.r = [srci.r(1)+0*r_y;r_y];
    % srci_r.n = [srci.n(1).*abs(r_signs);r_signs.*srci.n(2)]; 

    if iouts(i)
        % rhso = rhso + (rs_val(srci,sys.chnkr)+rs_val(srci_r,sys.chnkr)*r_signs.').*ssrc.charge(i);
        % vals(ioutt) = vals(ioutt) + (u_val(srci,targo)+ u_val(srci_r,targo)*r_signs.').*ssrc.charge(i);
        rhso = rhso + (rs_val(srci,sys.chnkr)).*ssrc.charge(i);
        vals(ioutt) = vals(ioutt) + (u_val(srci,targo)).*ssrc.charge(i);
    else
        % rhsi = rhsi + (rs_val(srci,sys.chnkr)+rs_val(srci_r,sys.chnkr)*r_signs.').*ssrc.charge(i);
        % vals(iiint) = vals(iiint) + (u_val(srci,targi)+ u_val(srci_r,targi)*r_signs.').*ssrc.charge(i);
        rhsi = rhsi + (rs_val(srci,sys.chnkr)).*ssrc.charge(i);
        vals(iiint) = vals(iiint) + (u_val(srci,targi)).*ssrc.charge(i);
    end
end

[iouts, iiins] = wave_kern.flagout(sys,dsrc.r);
for i = 1:ndsrc
    srci = [];
    srci.r = dsrc.r(:,i);
    srci.n = dsrc.n(:,i);

    % r_signs = [-1,-1,1,1];
    % r_y1 = [h - (srci.r(2)-h),-h - (srci.r(2)+h)];
    % r_y = [r_y1,h- (r_y1(2)-h),-h - (r_y1(1)+h)];
    % r_signs = r_signs.*(abs(r_y) < h);
    % 
    % srci_r = [];
    % srci_r.r = [srci.r(1)+0*r_y;r_y];
    % srci_r.n = [srci.n(1).*abs(r_signs);r_signs.*srci.n(2)]; 

    if iouts(i)
        % rhso = rhso + (rd_val(srci,sys.chnkr)+rd_val(srci_r,sys.chnkr)*r_signs.').*dsrc.charge(i);
        rhso = rhso + (rd_val(srci,sys.chnkr)).*dsrc.charge(i);
        % vals(ioutt) = vals(ioutt) + (ud_val(srci,targo)+ ud_val(srci_r,targo)*r_signs.').*dsrc.charge(i);
        vals(ioutt) = vals(ioutt) + (ud_val(srci,targo)).*dsrc.charge(i);
    else
        % rhsi = rhsi + (rd_val(srci,sys.chnkr)+rd_val(srci_r,sys.chnkr)*r_signs.').*dsrc.charge(i);
        % vals(iiint) = vals(iiint) + (ud_val(srci,targi)+ ud_val(srci_r,targi)*r_signs.').*dsrc.charge(i);
        rhsi = rhsi + (rd_val(srci,sys.chnkr)).*dsrc.charge(i);
        vals(iiint) = vals(iiint) + (ud_val(srci,targi)).*dsrc.charge(i);
    end
end

sigmai = sys.sysinvi(rhsi);
sigmao = sys.sysinvo(rhso);
opts = []; opts.eps = 1e-8;
vals_scati = chunkerkerneval(sys.chnkr,g_val,sigmai,targi.r,opts);
vals_scato = chunkerkerneval(sys.chnkr,g_val,sigmao,targo.r,opts);

vals(iiint) = vals(iiint) - vals_scati;
vals(ioutt) = vals(ioutt) - vals_scato;
% toc
end

