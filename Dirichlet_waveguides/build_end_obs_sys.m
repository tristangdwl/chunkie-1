function sys = build_end_obs_sys(nover,d,zk,do,ctr,zk2)



nch = 2.^(nover+2);
cparams = []; cparams.ta = 0; cparams.tb = pi;
chnkr_tmp = chunkerfuncuni(@flatinterface,nch,cparams);
chnkr_tmp = sort(chnkr_tmp);

chnkr1 = move(chnkr_tmp,[0; d/2]);
chnkr2 = move(chnkr_tmp,[0;-d/2]);
chnkr = merge([chnkr1,chnkr2]);

opts = []; opts.eps = 1e-8;
fkern_self = @(s,t) 0.25*1i*besselh(0,zk*sqrt((cos(s.r(1,:))-cos(t.r(1,:).')).^2+(s.r(2,:)-t.r(2,:).').^2));
sysmat_f = chunkermat(chnkr,fkern_self,opts);




cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 0;
cparams.maxchunklen = 4.0/max(zk,zk2)*2^(1-nover); % setting a chunk length helps when the
                              % frequency is known
                              
pref = []; 
pref.k = 16;
start = tic;
narms = 3;
chnkr_o = chunkerfunc(@(t) ellcfun(t,do/2,do),cparams,pref); 
% chnkr_o = chunkerfunc(@(t) starfish(t,narms,.5),cparams,pref); 
t1 = toc(start);
chnkr_o = move(chnkr_o,-ctr,[],[],do);

chnkr_o.n = normals(chnkr_o);


% chnkr = move(chnkr,[-7;0],[],[],[]);

fprintf('%5.2e s : time to build geo\n',t1)
optscomp = []; optscomp.type = 'flam'; optscomp.keep_mat = true;
sys_o = leaky.build_chan_mat(zk,zk2,chnkr_o,do,optscomp);
sys_o.shape = 'obstacle';

sysmat_o = sys_o.sysmat;

cosmap = @(s) struct("r",[cos(s.r(1,:));s.r(2,:)],"n",s.n);
cosmap2 = @(s) struct("r",[cos(s.r(1,:));s.r(2,:)]);
acosmap = @(s) struct("r",[acos(s.r(1,:));s.r(2,:)],"n",s.n);


fkern_s = @(s,t) chnk.helm2d.kern(zk,cosmap(s),t,'s',1);
fkern_sp = @(s,t) chnk.helm2d.kern(zk,cosmap(s),t,'sprime',1);


sysmat_of = -[smth_mat(fkern_s,chnkr,chnkr_o); -smth_mat(fkern_sp,chnkr,chnkr_o)];

fkern_s = @(s,t) chnk.helm2d.kern(zk,s,cosmap(t),'s',1);
fkern_d = @(s,t) chnk.helm2d.kern(zk,s,cosmap(t),'d',1);
sysmat_fo = [smth_mat(fkern_s,chnkr_o,chnkr), smth_mat(fkern_d,chnkr_o,chnkr)];


sysmat = [sysmat_f,sysmat_fo; sysmat_of,sysmat_o];

sys = [];
sys.sys_o = sys_o;
sys.chnkr_o = chnkr_o;
sys.do = do;
sys.zko = zk2;

pts = [cos(chnkr.r(1,:)); chnkr.r(2,:)];
sys.pts = struct("r",[pts]);


sys.d = d;
sys.chnkr = chnkr;
sys.zk = zk;

occ = 1024; rank_or_tol = opts.eps;
xxflam = repmat(real(chnkr_o.r(1,:)),1,2);
xyflam = repmat(real(chnkr_o.r(2,:)),1,2);
xflam = [real(chnkr.r(1,:)), xxflam(:).';real(chnkr.r(2,:)),xyflam(:).'];

tic;
F = rskelf(sysmat,xflam,occ,rank_or_tol); 
tflam = toc
% 
sys.sysinv = @(rhs) rskelf_sv(F,rhs);
sys.F = F;


end


function sysmat = smth_mat(fkern,chnkr1,chnkr2)
wts = chnkr1.weights; wts = wts(:).';

sysmat  = fkern(chnkr1,chnkr2).*wts;

end

function [r,d,d2] = ellcfun(t,a,b) 
% parameterization of an ellipse
r = [a*cos(t(:).');b*sin(t(:).')];
d = [-a*sin(t(:).');b*cos(t(:).')];
d2 = [-a*cos(t(:).');-b*sin(t(:).')];

end