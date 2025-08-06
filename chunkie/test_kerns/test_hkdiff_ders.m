src = []; src.r = [0;0]; 

targ = []; targ.r = [1;0.9]; 


ifr2logr = 0;
k = 1.2; h = 1e-3;

[valm,gradm,hessm,der3m,der4m,der5m] = chnk.flex2d.helmdiffgreen(k-h,src.r,targ.r,ifr2logr);
[valp,gradp,hessp,der3p,der4p,der5p] = chnk.flex2d.helmdiffgreen(k+h,src.r,targ.r,ifr2logr);

[valdk,graddk,hessdk,der3dk,der4dk,der5dk] = chnk.flex2d.helmdkdiffgreen(k,src.r,targ.r,ifr2logr);

assert(norm(-(valp-valm)/2/h - valdk*2*k)/norm(valp(:))<1e-5)
assert(norm(-(gradp(:)-gradm(:))/2/h - graddk(:)*2*k)/norm(gradp(:))<1e-5)
assert(norm(-(hessp(:)-hessm(:))/2/h - hessdk(:)*2*k)/norm(hessp(:))<1e-5)
assert(norm(-(der3p(:)-der3m(:))/2/h - der3dk(:)*2*k)/norm(der3p(:))<1e-5)
assert(norm(-(der4p(:)-der4m(:))/2/h - der4dk(:)*2*k)/norm(der4p(:))<1e-5)
assert(norm(-(der5p(:)-der5m(:))/2/h - der5dk(:)*2*k)/norm(der5p(:))<1e-5)



targ = []; targ.r = [0.2;0.3]; 

ifr2logr = 0;
k = 1; h = 1e-3;

[valm,gradm,hessm,der3m,der4m,der5m] = chnk.flex2d.helmdiffgreen(k-h,src.r,targ.r,ifr2logr);
[valp,gradp,hessp,der3p,der4p,der5p] = chnk.flex2d.helmdiffgreen(k+h,src.r,targ.r,ifr2logr);

[valdk,graddk,hessdk,der3dk,der4dk,der5dk] = chnk.flex2d.helmdkdiffgreen(k,src.r,targ.r,ifr2logr);

assert(norm(-(valp-valm)/2/h - valdk*2*k)/norm(valp(:))<1e-5)
assert(norm(-(gradp(:)-gradm(:))/2/h - graddk(:)*2*k)/norm(gradp(:))<1e-5)
assert(norm(-(hessp(:)-hessm(:))/2/h - hessdk(:)*2*k)/norm(hessp(:))<1e-5)
assert(norm(-(der3p(:)-der3m(:))/2/h - der3dk(:)*2*k)/norm(der3p(:))<1e-5)
assert(norm(-(der4p(:)-der4m(:))/2/h - der4dk(:)*2*k)/norm(der4p(:))<1e-5)
assert(norm(-(der5p(:)-der5m(:))/2/h - der5dk(:)*2*k)/norm(der5p(:))<1e-5)

src = []; src.r = [0;0]; 

targ = []; targ.r = [1;1]; 


ifr2logr = 0;
k = 1.2; h = 1e-3;

[valm,gradm,hessm,der3m,der4m,der5m] = chnk.flex2d.helmdkdiffgreen(k-h,src.r,targ.r,ifr2logr);
% [val,grad,hess,der3,der4,der5] = chnk.flex2d.helmdiffgreen(k,src.r,targ.r,ifr2logr);
[valp,gradp,hessp,der3p,der4p,der5p] = chnk.flex2d.helmdkdiffgreen(k+h,src.r,targ.r,ifr2logr);

[valdk,graddk,hessdk,der3dk,der4dk,der5dk] = chnk.flex2d.helmdkdkdiffgreen(k,src.r,targ.r,ifr2logr);

assert(norm(-0.5*(valp-valm)/2/h - valdk*2*k)/norm(valp(:))<1e-5)
assert(norm(-0.5*(gradp(:)-gradm(:))/2/h - graddk(:)*2*k)/norm(gradp(:))<1e-5)
assert(norm(-0.5*(hessp(:)-hessm(:))/2/h - hessdk(:)*2*k)/norm(hessp(:))<1e-5)
assert(norm(-0.5*(der3p(:)-der3m(:))/2/h - der3dk(:)*2*k)/norm(der3p(:))<1e-5)
assert(norm(-0.5*(der4p(:)-der4m(:))/2/h - der4dk(:)*2*k)/norm(der4p(:))<1e-5)
assert(norm(-0.5*(der5p(:)-der5m(:))/2/h - der5dk(:)*2*k)/norm(der5p(:))<1e-5)

% norm((valp-2*val+valm)/h^2 - valdk*2*k)/norm(valp(:))
% norm((gradp(:)-2*grad(:)+gradm(:))/h^2 - graddk(:)*2*k)/norm(gradp(:))



%%
k = 1.3;
h = 2e-4;
targ = []; targ.r = [-3;2];
% targ.r = [1;0];
targh = []; targh.r = targ.r + [h,0,-h].*[1;0]; 

[val,grad,hess,der3,der4,der5] = chnk.flex2d.helmdkdiffgreen(k,src.r,targh.r,0);
sixth = chnk.flex2d.green_helmsq_sixth(k,src.r,targ.r);

uh = (der5(1,:) - der5(3,:))/2/h;
% c = mean(uh./ sixth(1,1:6))
% c= 25/36; c =1;
norm(uh - sixth(1,1:6))

% uh = (der4(1,:) - der4(3,:))/2/h;
% norm(uh - der5(2,1:5))

%%
h = 1e-3;
targ = []; targ.r = [1;0.5]; targ.r = [1;0];
targh = []; targh.r = targ.r + [[0;0],[h;0],[-h;0],[0;h],[0;-h]]; 
val = chnk.flex2d.helmdkdiffgreen(k,src.r,targh.r,0);
val2 = chnk.helm2d.green(k,src.r,targ.r);

uh = (sum(val(2:end)-val(1)))/h^2 + k^2*val(1);
norm(val2 - uh)

%%
val = chnk.flex2d.helmdkdkdiffgreen(k,src.r,targh.r,0);
val2 = chnk.flex2d.helmdkdiffgreen(k,src.r,targ.r,0);

uh = (sum(val(2:end)-val(1)))/h^2 + k^2*val(1);
norm(val2 - uh)





