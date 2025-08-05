function submat= kern_alt(zk,srcinfo,targinfo,type,varargin)
%FLEX2D.KERN_ALT flexural wave kernels in 2D
% 
% Syntax: submat = chnk.flex2d.kern_alt(zk,srcinfo,targingo,type,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). Note that the normal information is obtained
% by taking the perpendicular to the provided tangential deriviative
% info and normalizing  
%  
% Kernels for the equation (Detla + k^2)^2, based on the Green's function:
%         G(x,y)=(-i/8k) H_1^(1)(k|x-y|)|x-y|
%


  
src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);


if strcmpi(type, 's') % flexural wave single layer

   val = chnk.flex2d.hkdiffgreen(zk,src,targ);  
   submat = 1/(2*zk^2).*val;

end

