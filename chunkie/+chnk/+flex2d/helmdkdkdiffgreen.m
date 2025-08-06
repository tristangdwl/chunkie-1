function [val,grad,hess,der3,der4,der5] = helmdkdkdiffgreen(k,src,targ,ifr2logr)
%HELMDKDIFFGREEN evaluate the difference of the 
% Helmholtz Green function and the Laplace Green function
% for the given sources and targets, i.e. 
%
% G(x,y) = d^2/d(k^2)^2 [(i/4 H_0^(1)(k|x-y|) + 1/(2 pi) log(|x-y|))]
%
% or the difference of the Helmholtz and Laplace Green funcions 
% and k^2 r^2 log r/ 8 pi (a constant times the biharmonic Green function)
% i.e. 
%
% G(x,y) =  d^2/d(k^2)^2[i/4 H_0^(1)(k|x-y|) + 1/(2 pi) log(|x-y|) + ...
%                    - k^2/(8*pi) |x-y|^2 log(|x-y|)]
%
% where H_0^(1) is the principal branch of the Hankel function
% of the first kind. This routine avoids numerical cancellation
% when |k||x-y| is small.
%
% - grad(:,:,1) has G_{x1}, grad(:,:,2) has G_{x2}
% - hess(:,:,1) has G_{x1x1}, hess(:,:,2) has G_{x1x2}, 
% hess(:,:,3) has G_{x2x2}
% - der3 has the third derivatives in the order G_{x1x1x1}, G_{x1x1x2}, 
% G_{x1x2x2}, G_{x2x2x2}
% - der4 has the fourth derivatives in the order G_{x1x1x1x1}, 
% G_{x1x1x1x2}, G_{x1x1x2x2}, G_{x1x2x2x2}, G_{x2x2x2x2}
%
% derivatives are on the *target* variables
%
% input:
%
% src - (2,ns) array of source locations
% targ - (2,nt) array of target locations
% k - wave number, as above
%
% optional input:
%
% ifr2logr - boolean, default: false. If true, also subtract off the 
%             k^2/(8pi) r^2 log r kernel

if nargin < 4
    ifr2logr = false;
end

r2logrfac = 1;
if ifr2logr
    r2logrfac = 0;
end

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

dx = xt-xs;
dy = yt-ys;

dx2 = dx.*dx;
dx3 = dx2.*dx;
dx4 = dx3.*dx;
dx5 = dx4.*dx;

dy2 = dy.*dy;
dy3 = dy2.*dy;
dy4 = dy3.*dy;
dy5 = dy4.*dy;

r2 = dx2 + dy2;
r = sqrt(r2);
rm1 = 1./r;
rm2 = rm1.*rm1;
rm3 = rm1.*rm2;
rm4 = rm1.*rm3;
rm5 = rm1.*rm4;

% get value and r derivatives
      
[g0,g1,g21,g3,g4,g5] = diff_h0log_and_rders(k,r,r2logrfac);
g0 = g0/2; g1 = g1/2; g21 = g21/2; g3 = g3/2; g4 = g4/2; g5 = g5/2;
%     evaluate potential and derivatives

if nargout > 0
    val = g0;  
end
if nargout > 1
    grad(:,:,1) = dx.*g1.*rm1;
    grad(:,:,2) = dy.*g1.*rm1;
end
if nargout > 2
    hess(:,:,1) = dx2.*g21.*rm2+g1.*rm1;
    hess(:,:,2) = dx.*dy.*g21.*rm2;
    hess(:,:,3) = dy2.*g21.*rm2+g1.*rm1;
end
if nargout > 3
    der3(:,:,1) = (dx3.*g3+3*dy2.*dx.*g21.*rm1).*rm3;
    der3(:,:,2) = dx2.*dy.*(g3.*rm3-3*g21.*rm4) + ...
             dy.*g21.*rm2;
    der3(:,:,3) = dx.*dy2.*(g3.*rm3-3*g21.*rm4) + ...
             dx.*g21.*rm2;
    der3(:,:,4) = (dy3.*g3+3*dx2.*dy.*g21.*rm1).*rm3;
end

if nargout > 4
    der4(:,:,1) = (dx4.*(g4-6*g3.*rm1+15*g21.*rm2)).*rm4 + ...
             (6*dx2.*(g3-3*g21.*rm1)).*rm3 + ...
             3*g21.*rm2;
    der4(:,:,2) = (dx3.*dy.*(g4-6*g3.*rm1+15*g21.*rm2)).*rm4 + ...
             (3*dx.*dy.*(g3-3*g21.*rm1)).*rm3;
    der4(:,:,3) = dx2.*dy2.*(g4-6*g3.*rm1+15*g21.*rm2).*rm4 + ...
             g3.*rm1 - 2*g21.*rm2;
    der4(:,:,4) = dx.*dy3.*(g4-6*g3.*rm1+15*g21.*rm2).*rm4 + ...
             3*dx.*dy.*(g3-3*g21.*rm1).*rm3;
    der4(:,:,5) = dy4.*(g4-6*g3.*rm1+15*g21.*rm2).*rm4 + ...
             6*dy2.*(g3-3*g21.*rm1).*rm3 + ...
             3*g21.*rm2;
end

if nargout > 5
    der5(:,:,1) = (dx5.*g5+10*dy2.*dx3.*g4.*rm1 + ...
          (15*dy4.*dx-30*dy2.*dx3).*g3.*rm2 + ...
         (60*dy2.*dx3-45*dy4.*dx).*g21.*rm3).*rm5;
    der5(:,:,2) = (dy.*dx4.*g5+(6*dy3.*dx2-4*dy.*dx4).*g4.*rm1 + ...
      (3*dy5+12*dy.*dx4-30*dy3.*dx2).*g3.*rm2 + ...
     (72*dy3.*dx2-9*dy5-24*dy.*dx4).*g21.*rm3).*rm5;
    der5(:,:,3) = (dy2.*dx3.*g5+(3*dy4.*dx-6*dy2.*dx3+dx5).*g4.*rm1 + ...
      (27*dy2.*dx3-15*dy4.*dx-3*dx5).*g3.*rm2 + ...
     (36*dy4.*dx-63*dy2.*dx3+6*dx5).*g21.*rm3).*rm5;
    der5(:,:,4) = (dx2.*dy3.*g5+(3*dx4.*dy-6*dx2.*dy3+dy5).*g4.*rm1 + ...
      (27*dx2.*dy3-15*dx4.*dy-3*dy5).*g3.*rm2 + ...
     (36*dx4.*dy-63*dx2.*dy3+6*dy5).*g21.*rm3).*rm5;
    der5(:,:,5) = (dx.*dy4.*g5+(6*dx3.*dy2-4*dx.*dy4).*g4.*rm1 + ...
      (3*dx5+12*dx.*dy4-30*dx3.*dy2).*g3.*rm2 + ...
     (72*dx3.*dy2-9*dx5-24*dx.*dy4).*g21.*rm3).*rm5;
    der5(:,:,6) = (dy5.*g5+10*dx2.*dy3.*g4.*rm1 + ...
      (15*dx4.*dy-30*dx2.*dy3).*g3.*rm2 + ...
     (60*dx2.*dy3-45*dx4.*dy).*g21.*rm3).*rm5;
end

end

function [g0,g1,g21,g3,g4,g5] = diff_h0log_and_rders(k,r,r2logrfac)
% g0 = g
% g1 = g'
% g21 = g'' - g'/r
%
% maybe later:
% g321 = g''' - 3*g''/r + 3g'/r^2
% g4321 = g'''' - 6*g'''/r + 15*g''/r^2 - 15*g'/r^3

io4 = 1i*0.25;
o2p = 1/(2*pi);


kr = k*r;

h0 = besselh(0,1,kr);
h1 = besselh(1,1,kr);
h2 = besselh(2,1,kr);
h3 = besselh(3,1,kr);

dh0da = -h1.*r/2/k;
dh1da = (h0-h2).*r/4/k;

dh0daa = r.*h1./4/k^3 - r.^2.*(h0-h2)/8/k^2;
dh1daa = (-r/16/k^3) .* (2*h0 + 3*k*r.*h1 - 2*h2 -k*r.*h3);

rm1 = 1./r;
rm2 = rm1.*rm1;
rm3 = rm1.*rm2;
rm4 = rm1.*rm3;
rm5 = rm1.*rm4;

r2fac = (1-r2logrfac)*0.25*o2p;
logr = log(r);
g0 = io4*dh0daa - r2fac*r.*r.*logr;
g1 = -k*io4*dh1daa - 2*io4*dh1da/2/k + io4*h1/k^3/4;
g21 = -2*io4*dh0da -k*k*io4*dh0daa + k*io4*dh1daa.*rm1 + 2*io4*dh1da.*rm1/2/k - io4*h1.*rm1/4/k^3 - ...
    g1.*rm1;
g3 = k*k*io4*dh0daa.*rm1 + 2*io4*dh0da.*rm1 + io4*(1/k/2+1/k/2-(k*k-2*rm2)/k^3/4).*h1 + io4*(k+(k*k-2*rm2)/k/2).*dh1da ...
    + io4*((k*k-2*rm2)/2/k + k).*dh1da + io4*k*(k*k-2*rm2).*dh1daa;
g4 = io4*((-(3*rm2-k*k)/4/k^3 - 1/k).*(2*h1.*rm1-k*h0) + 2*((3*rm2-k*k)/2/k-k).*(2*dh1da.*rm1-k*dh0da-h0/2/k)+k*(3*rm2-k*k).*(2*dh1daa.*rm1-k*dh0daa-2*dh0da/2/k+h0/4/k^3));

g5 = io4*((-4*rm1).*h0+2*(12*rm3-4*k^2*rm1).*dh0da+(12*k^2*rm3-2*k^4*rm1).*dh0daa  ...
    +(-(-24*rm4+7*k*k*rm2-k^4)/4/k^3 + (7*rm2-2*k^2)/2/k + (7*rm2-2*k^2)/2/k -2*k).*h1+2*((-24*rm4+7*k*k*rm2-k^4)/2/k + k*(7*rm2-2*k^2)).*dh1da+k*(-24*rm4+7*k*k*rm2-k^4).*dh1daa);

end

