
function [dgrad] = helm_graddx(zk, s, t)
    src = s.r;
    targ = t.r;
    [~, ~, hess] = chnk.helm2d.green(zk, src, targ);
    [~, nt] = size(targ);
    [~, ns] = size(src);

    srcnorm = s.n;
    nx = repmat(srcnorm(1,:), nt, 1);
    ny = repmat(srcnorm(2,:), nt, 1);


    dgrad = -(hess(:,:,1).*nx + hess(:,:,2).*ny);
end