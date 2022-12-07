%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   .  .  .  builds a simple pentagonal chunkergraph 
%            and tests the interior Helmholtz Dirichlet problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

vertsx = [0,1,1,-1,-1,-1/sqrt(2),0,1/sqrt(2)];
vertsy = [0,0,2, 2, 0, 1/sqrt(2),2/sqrt(2),1/sqrt(2)];
verts = [vertsx;vertsy];

edge2verts = [-1, 1, 0, 0, 0, 0, 0, 0; ...
               0,-1, 1, 0, 0, 0, 0, 0; ...
               0, 0,-1, 1, 0, 0, 0, 0; ...
               0, 0, 0,-1, 1, 0, 0, 0; ...
               1, 0, 0, 0,-1, 0, 0, 0; ...
              -1, 0, 0, 0, 0, 1, 0, 0; ...
               0, 0, 0, 0, 0,-1, 1, 0; ...
               0, 0, 0, 0, 0, 0,-1, 1; ...
               1, 0, 0, 0, 0, 0, 0,-1];
edge2verts = sparse(edge2verts);


fchnks    = {};

prefs      = [];
prefs.chsmall = 1d-4;
[cgrph] = chunkgraphinit(verts,edge2verts,fchnks,prefs);

vstruc = procverts(cgrph);
rgns = findregions(cgrph);
cgrph = balance(cgrph);

zk = 1.0;
fkern = @(s,t) chnk.helm2d.kern(zk,s,t,'d');

opts = [];
[sysmat] = chunkermat(cgrph,fkern,opts);
sysmat = sysmat - eye(size(sysmat,2))/2;

% generate some targets...

xs = -1:0.01:1;
ys =  0:0.01:2;
[X,Y] = meshgrid(xs,ys);
targs = [X(:).';Y(:).'];

srcinfo = [];
srcinfo.sources = cgrph.r(:,:);
w = weights(cgrph);
n = normals(cgrph);

% a quick hack to find the interior points

srcinfo.dipstr = w(:).';
srcinfo.dipvec = n(:,:); 
eps = 1E-8;
pg  = 0;
pgt = 1;
[U] = lfmm2d(eps,srcinfo,pg,targs,pgt);
U = U.pottarg;
inds = find(abs(U-2*pi)<pi/10);

%%%%%%%%%%%%%%%%%%
% generate the right hand side

x0 = 1.3;
y0 = 0.9;

srcinfo = [];
srcinfo.sources = cgrph.r(:,:);
w = weights(cgrph);
n = normals(cgrph);

s = [];
s.r = [x0;y0];
t = [];
t.r = cgrph.r(:,:);
rhs = chnk.helm2d.kern(zk,s,t,'s');
dens = sysmat\rhs;

srcinfo.dipstr = (w(:).*dens(:)).';
srcinfo.dipvec = n(:,:); 
fints = hfmm2d(eps,zk,srcinfo,pg,targs(:,inds),pgt);

t.r = targs(:,inds);
true_sol = chnk.helm2d.kern(zk,s,t,'s');

usol = zeros(size(targs,2),1);
uerr = usol;
usol(inds) = fints.pottarg;
uerr(inds) = usol(inds)-true_sol;
usol = reshape(usol,size(X));
uerr = reshape(uerr,size(X));