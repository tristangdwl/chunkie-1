function sys = build_tube(d,zk,nover_match,L_chnkr,nover_tube,nover_slit)

if nargin < 5, nover_tube = 2; end
if nargin < 6, nover_slit = 3; end


sys_L = build_dirichlet_sys(nover_tube,d,L_chnkr,zk);
sys_R = build_two_flats_sys(nover_slit,d,zk);



chnkr = vert_chnkr(nover_match,L_chnkr);
chnkr = sort(chnkr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   .  .  .  single double layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
fkern_S = @(s,t) d_guide_kern(s,t,sys_L,'s')-d_two_flat(s,t,sys_R,'s');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
sysmats = chunkermat(chnkr,fkern_S,opts);
% sysmats = fkern_S(chnk_struc,chnk_struc).*wts(:).';
t1 = toc(start)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   .  .  .  derivative double layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fkern_Sp = @(s,t) d_guide_kern(s,t,sys_L,'sp')-d_two_flat(s,t,sys_R,'sp');

fkern_Dp = @(s,t) d_guide_kern(s,t,sys_L,'dp')-d_two_flat(s,t,sys_R,'dp');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; 
sysmatdd = chunkermat(chnkr,fkern_Dp,opts);
% sysmatdd = fkern_Dp(chnk_struc,chnk_struc).*wts(:).';
t1 = toc(start)
sysmatd = zeros(chnkr.npt);
sysmatsd = zeros(chnkr.npt); 


sysmat = [sysmatd,sysmats;-sysmatdd,-sysmatsd];
sysmat = sysmat-eye(size(sysmat,1));


occ = 1024;
rank_or_tol = 1e-8;

xxflam = repmat(real(chnkr.r(1,:)),1,2);
xyflam = repmat(real(chnkr.r(2,:)),1,2);
xflam = [xxflam(:).';xyflam(:).'];
tic;
F = rskelf(sysmat,xflam,occ,rank_or_tol); 
tflam = toc
sysinv = @(rhs) rskelf_sv(F,rhs);


sys = [];
sys.d = d;
sys.zk = zk;
sys.sys_L = sys_L;
sys.sys_R = sys_R;

sys.chnkr = chnkr;
sys.sysmat = sysmat;
sys.sysinv = sysinv;

pts_L = real(sys_L.chnkr.r(:,:)); pts_L = pts_L(:,pts_L(1,:)<0);
pts_R = real([cos(sys_R.chnkr.r(1,:));sys_R.chnkr.r(2,:)]); pts_R = pts_R(:,pts_R(1,:)>0);

sys.pts = [pts_L,pts_R];

end

