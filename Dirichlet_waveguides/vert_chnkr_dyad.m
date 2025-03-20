function chnkr = vert_chnkr_dyad(nover,h,L)
if nargin < 2
    L = 40;
elseif isempty(L)
    L = 40;
end



cparams.ta = -L-8;
cparams.tb =  -h;
cparams.chsmall = 1e-2;
cparams.ifclosed = 0;
cpars = [];
cpars.L = L;
cpars.c2 = 5;
cpars.c1 = 20;
nch = ceil(48*2^(nover-1));
chnkrd = chunkerfunc(@(t) complexx2(t,cpars),...
   cparams);

cparams.ta = -h;
cparams.tb =  h;
cpars = [];
cpars.L = L;
cpars.c2 = 5;
cpars.c1 = 20;
nch = ceil(48*2^(nover-1));
chnkrm = chunkerfunc(@(t) complexx2(t,cpars),...
   cparams);

cparams.ta = h;
cparams.tb = L+8;
cpars = [];
cpars.L = L;
cpars.c2 = 5;
cpars.c1 = 20;
nch = ceil(48*2^(nover-1));
chnkru = chunkerfunc(@(t) complexx2(t,cpars),...
   cparams);


chnkr = merge([chnkrd,chnkrm,chnkru]);

chnkr.r  = flipud(chnkr.r);
chnkr.d  = flipud(chnkr.d);
chnkr.d2 = flipud(chnkr.d2);
chnkr.n = normals(chnkr);


end