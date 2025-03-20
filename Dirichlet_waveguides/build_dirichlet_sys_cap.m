function sys = build_dirichlet_sys_cap(nover,d,igeom,zk)

if igeom == 1
cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 0;
cparams.maxchunklen = 4.0/zk*2^(-nover); % setting a chunk length helps when the
                              % frequency is known
                              
pref = []; 
pref.k = 16;
chnkr = chunkerfunc(@(t) ellcfun(t,d,d/2),cparams,pref); 
else

vertsL = [0.805875576036866   0.386
0.557027649769585   0.386
0.490207373271889   0.272262773722628
% 0.589285714285714   0.164233576642336
% 0.563940092165899   0.123357664233577
0.448732718894009   0.225547445255475
0.398041474654378   0.394890510948905
0.236751152073733   0.161313868613139
0.195276497695852   0.237226277372263
0.324308755760369   0.409489051094891
0.294354838709677   0.435766423357664
0.063940092165899   0.386131386861314
0.080069124423963   0.464963503649635
0.236751152073733   0.491240875912409
0.324308755760369   0.555474452554745
0.301267281105991   0.625547445255475
0.057027649769585   0.678102189781022
0.061635944700461   0.745255474452555
0.309605937048494   0.704520252752197
0.301267281105991   0.789051094890511
0.255184331797235   0.853284671532847
0.262096774193548   0.937956204379562
0.377304147465438   0.818248175182482
% 0.375000000000000   0.540875912408759
% 0.485599078341014   0.462043795620438
% 0.517857142857143   0.590510948905110
% 0.508640552995392   0.724817518248175
% 0.467165898617512   0.876642335766423
% 0.529377880184332   0.829927007299270
% 0.598502304147465   0.943795620437956
% 0.628456221198157   0.905839416058394
% 0.557027649769585   0.759854014598540
0.577764976958525   0.54
0.805875576036866   0.54].';


vertsL(1,:) = vertsL(1,:) - vertsL(1,1);
vertsL(2,:) = vertsL(2,:) - vertsL(2,1);
xmax = max(abs(vertsL(1,:)));
ymax = max(abs(vertsL(2,:)));

vertsL(1,:) = ymax/xmax*vertsL(1,:);

d_old = vertsL(2,end);
d_new = d;
vertsL(2,:) = (d_new/d_old)*(vertsL(2,:)-d_old/2);
vertsL(1,:) = (d_new/d_old)*vertsL(1,:);

vertsR = fliplr(vertsL); 
vertsR(1,:) = -vertsR(1,:); vertsR = vertsR(:,2:end-1);

verts = fliplr([vertsL,vertsR]);


cparams = [];
cparams.eps = 1.0e-3;
                              
pref = []; 
pref.k = 16;
cparams.autowidths = true;
cparams.autowidthsfac = .4;
chnkr = chunkerpoly(verts,cparams,pref);
opts = [];
opts.lvlr = 'a';
opts.maxchunklen = 4.0/zk*2^(-nover);
chnkr = refine(chnkr,opts);
chnkr.n = normals(chnkr);

end

chnkr.n = normals(chnkr);

sys = [];
sys.d = d;
sys.zk = zk;
sys.chnkr = chnkr;

fkern = @(s,t) chnk.helm2d.kern(sys.zk,s,t,'d');
opts = []; opts.eps = 1e-8;
sysmati = chunkermat(chnkr,fkern,opts)-eye(chnkr.npt)/2;
sysmato = chunkermat(chnkr,fkern,opts)+eye(chnkr.npt)/2;
occ = 1024; rank_or_tol = opts.eps;
tic;
Fi = rskelf(sysmati,real(chnkr.r(:,:)),occ,rank_or_tol); 
Fo = rskelf(sysmato,real(chnkr.r(:,:)),occ,rank_or_tol); 
tflam = toc
% 
sys.sysinvi = @(rhs) rskelf_sv(Fi,rhs);
sys.sysinvo = @(rhs) rskelf_sv(Fo,rhs);
sys.Fi = Fi;
sys.Fo = Fo;

sys.shape = 'cap';
end

function [r,d,d2] = ellcfun(t,a,b) 
% parameterization of an ellipse
r = [a*cos(t(:).');b*sin(t(:).')];
d = [-a*sin(t(:).');b*cos(t(:).')];
d2 = [-a*cos(t(:).');-b*sin(t(:).')];

end