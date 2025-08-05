phi = 0.3;
src = []; src.r = [0;0]; src.n = [cos(phi);sin(phi)]; src.d = [cos(phi+pi/2);sin(phi+pi/2)];

targ = []; targ.r = [1;1]; targ.n = [1;1]/sqrt(2); targ.d = [1;-1]/sqrt(2);
kappa = 0.282842712474619;


ifr2logr = 0;
k = 1.2; h = 1e-3;

[valm,gradm,hessm,der3m,der4m,der5m] = chnk.flex2d.helmdiffgreen(k-h,src.r,targ.r,ifr2logr);
[valp,gradp,hessp,der3p,der4p,der5p] = chnk.flex2d.helmdiffgreen(k+h,src.r,targ.r,ifr2logr);

[valdk,graddk,hessdk,der3dk,der4dk,der5dk] = chnk.flex2d.helmdkdiffgreen(k,src.r,targ.r,ifr2logr);

norm((valp-valm)/2/h - valdk*2*k)/norm(valp(:))
norm((gradp(:)-gradm(:))/2/h - graddk(:)*2*k)/norm(gradp(:))
norm((hessp(:)-hessm(:))/2/h - hessdk(:)*2*k)/norm(hessp(:))
norm((der3p(:)-der3m(:))/2/h - der3dk(:)*2*k)/norm(der3p(:))
norm((der4p(:)-der4m(:))/2/h - der4dk(:)*2*k)/norm(der4p(:))
norm((der5p(:)-der5m(:))/2/h - der5dk(:)*2*k)/norm(der5p(:))