function chnkr = vert_chnkr2(nover,d,L)
if nargin < 2
    L = 40;
elseif isempty(L)
    L = 40;
end

cparams.ta = d/2;
cparams.tb =  L+8;
cpars = [];
cpars.L = L;
cpars.c2 = 5;
cpars.c1 = 20;
nch = ceil(22*2^(nover-1));
chnkru = chunkerfuncuni(@(t) complexx2(t,cpars),...
   nch,cparams);

cparams.ta = -L-8;
cparams.tb =  -d/2;
cpars = [];
cpars.L = L;
cpars.c2 = 5;
cpars.c1 = 20;
nch = ceil(22*2^(nover-1));
chnkrd = chunkerfuncuni(@(t) complexx2(t,cpars),...
   nch,cparams);


chnkr = merge([chnkrd,chnkru]);

chnkr.r  = flipud(chnkr.r);
chnkr.d  = flipud(chnkr.d);
chnkr.d2 = flipud(chnkr.d2);
chnkr.n = normals(chnkr);


end