function comp_data = comp_far(sys,L_comp,tol)
iempty = 0;
if isfield(sys,'iempty'), iempty= sys.iempty; end

if iempty
    comp_data = [];
else
    tic;
    skel_out = get_outskel(sys,tol,L_comp);
    skel_in  = get_inskel(sys,tol,L_comp);
    toc
    
    tic;
    comp_data = [];
    comp_data.f2f = comp_far2far(sys,skel_out,skel_in);
    toc

end
end

function comp_f2f_dat= comp_far2far(sys,skel_out,skel_in)
comp_f2f_dat = [];
comp_f2f_dat.ssrc = skel_out.ssrc;
comp_f2f_dat.dsrc = skel_out.dsrc;
comp_f2f_dat.nd = skel_out.nd;
comp_f2f_dat.ns = skel_out.ns;

comp_f2f_dat.sdat = skel_in.sdat;
comp_f2f_dat.spdat = skel_in.spdat;
comp_f2f_dat.nsdat = skel_in.nsdat;
comp_f2f_dat.nspdat = skel_in.nspdat;

comp_f2f_dat.dat2str = skel_out.IX*sys.sysinv(skel_in.IT);
end

function skel_out = get_outskel(sys,tol,L_comp)
if nargin < 3, L_comp = 40; end

chnkr = vert_chnkr(2,L_comp);

chnkr = move(chnkr,[],[-1/sys.zk1 + min(real(sys.chnkr.r(1,:)));0]);
targ = chnkr.r(:,:);


nsrc = sys.chnkr.npt;


C = zeros(size(targ,2),2*nsrc);

opts = [];
opts.eps = tol;
opts.isp = 1;
opts.forcesmooth = true;
skern = @(s,t) chnk.helm2d.kern(sys.zk1,s,t,'s');
C(:,(nsrc+1):end) = chunkerkernevalmat(sys.chnkr,skern,targ,opts);

dkern = @(s,t) chnk.helm2d.kern(sys.zk1, s, t, 'd');
C(:,1:nsrc) = chunkerkernevalmat(sys.chnkr,dkern,targ,opts);

[sk2,rd2,X] = id(C,tol);
nskel2 = length(sk2);
IX = zeros(2*nsrc,nskel2);
IX(sk2,:) = eye(nskel2);
IX(rd2,:) = X.';
IX = IX.';

nd = sum(sk2<=nsrc);
ns = sum(sk2> nsrc);
[sk2,I] = sort(sk2);
IX = IX(I,:);

% dat = @(src) leaky.kern_data_dp(sys, src,[1;0]);
% getent = @(A,i) A(i);
% data_skel = @(src) getent(dat(src),sk);
% e = C(:,sk2)*SXT*data_skel([1;10]);

w = sys.chnkr.weights;
ssrc = [];
ssrc.r = sys.chnkr.r(:,sk2(sk2>nsrc)-nsrc);
ssrc.n = sys.chnkr.n(:,sk2(sk2>nsrc)-nsrc);
IX(sk2>nsrc,:) = w(sk2(sk2>nsrc)-nsrc).'.*IX(sk2>nsrc,:);

dsrc = [];
dsrc.r = sys.chnkr.r(:,sk2(sk2<=nsrc));
dsrc.n = sys.chnkr.n(:,sk2(sk2<=nsrc));
dsrc.w = w(sk2(sk2<=nsrc));
IX(sk2<=nsrc,:) = w(sk2(sk2<=nsrc)).'.*IX(sk2<=nsrc,:);

skel_out = [];
skel_out.IX = IX;
skel_out.ssrc = ssrc;
skel_out.dsrc = dsrc;
skel_out.nd = nd;
skel_out.ns = ns;
% skel_out.is = (sk2>nsrc);
% skel_out.id = (sk2<=nsrc);



end

function skel_in = get_inskel(sys,tol,L_comp)
if nargin < 3, L_comp = 10; end
if nargin < 2, tol = 1e-10; end

chnkr = vert_chnkr(2,L_comp);
chnkr = move(chnkr,[],[-1/sys.zk1 + min(real(sys.chnkr.r(1,:)));0]);
src = chnkr.r(:,:);




nsrc = sys.chnkr.npt;

proxmat = leaky.kern_data_full(sys, src).';

[sk, rd, T] = id(proxmat,tol);
nskel = length(sk);




IT = zeros(2*nsrc,nskel);
IT(sk,:) = eye(nskel);
IT(rd,:) = T.';
% 
% src2 = [0;30];
% dat = leaky.kern_data(sys, src2).';
% dat_skel = dat(:,sk).';
% 
% vecnorm(dat(:) - IT*dat_skel)



[sk,I] = sort(sk);
IT = IT(:,I);
nsdat  = sum(sk<=nsrc);
nspdat = sum(sk> nsrc);

sdat = [];
sdat.r = sys.chnkr.r(:,sk(sk<=nsrc));

spdat = [];
spdat.r = sys.chnkr.r(:,sk(sk>nsrc)-nsrc);
spdat.n = sys.chnkr.n(:,sk(sk>nsrc)-nsrc);


skel_in =[];
skel_in.sdat = sdat;
skel_in.spdat = spdat;
skel_in.IT = IT;

skel_in.nsdat = nsdat;
skel_in.nspdat = nspdat;
% skel_in.is1 = (sk<=nsrc);
% skel_in.isp1 = (sk>nsrc);


end
