function [sgrad] = helm_gradsx(zk, s, t)
    src = s.r;
    targ = t.r;
    [~, grad] = chnk.helm2d.green(zk, src, targ);
    [~, nt] = size(targ);
    [~, ns] = size(src);
    sgrad= grad(:,:,1);
end