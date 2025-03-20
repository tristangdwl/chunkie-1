function vals = tube_kern(src,targ,sys)

chnkr = sys.chnkr;
sys_L = sys.sys_L;
sys_R = sys.sys_R;
sysmat = sys.sysmat;


il = (targ.r(1,:)<chnkr.r(1,1));
ir = ~il;

targl = [];
targl.r = targ.r(:,il);
targl.n = targ.n(:,il);

targr = [];
targr.r = targ.r(:,ir);
targr.n = targ.n(:,ir);

ntarg = size(targ.r,2);
nsrc = size(src.r,2);

vals = zeros(ntarg,nsrc);

for i = 1:nsrc
    srci = [];
    srci.r = src.r(:,i);
    srci.n = src.n(:,i);

    if srci.r(1,1) < 0, disp('invalid source location'), end

rhs = -[d_two_flat(srci,chnkr,sys_R,'s');-d_two_flat(srci,chnkr,sys_R,'sp')];


soln = sysmat\rhs;
ddens = soln(1:chnkr.npt);
sdens = soln(chnkr.npt+(1:chnkr.npt));


%%
uin = zeros(ntarg,1);
uscat = zeros(ntarg,1);

uin(ir) = d_two_flat(srci,targr,sys_R,'s');
%%

w = chnkr.weights; w = w(:);
ssrc = [];
ssrc.r = chnkr.r(:,:);
ssrc.n = chnkr.n(:,:);
ssrc.charge = sdens.*w;

dsrc = [];
dsrc.r = chnkr.r(:,:);
dsrc.n = chnkr.n(:,:);
dsrc.charge = ddens.*w;

uscat(il) = d_guide_kern_lump(ssrc,dsrc,targl,sys_L);
uscat(ir) = d_two_flat_kern_lump(ssrc,dsrc,targr,sys_R);

vals(:,i) = uin - uscat;

end



end