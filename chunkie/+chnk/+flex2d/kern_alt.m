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

%%% SOURCES

if strcmpi(type, 's') % flexural wave single layer

   submat = chnk.flex2d.helmdkdiffgreen(zk,src,targ);  

end

%%% BCs

% boundary conditions applied to a point source
if strcmpi(type, 'clamped_plate_bcs')
    nxtarg = targinfo.n(1,:).'; 
    nytarg = targinfo.n(2,:).';  
    submat = zeros(2*nt,ns);
    
    [val, grad] = chnk.flex2d.helmdkdiffgreen(zk, src, targ);
    
    firstbc = val ;
    secondbc = (grad(:, :, 1).*nxtarg + grad(:, :, 2).*nytarg);
   
    submat(1:2:end,:) = firstbc;
    submat(2:2:end,:) = secondbc;
end

% boundary conditions applied to a point source
if strcmpi(type, 'free_plate_bcs')
    targnorm = targinfo.n;
    targtang = targinfo.d;
    
    [val, ~, hess, third] = chnk.flex2d.helmdkdiffgreen(zk, src, targ);

    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    
    dx1 = repmat((targtang(1,:)).',1,ns);
    dy1 = repmat((targtang(2,:)).',1,ns);
    
    ds1 = sqrt(dx1.*dx1+dy1.*dy1); 
    
    tauxtarg = dx1./ds1;
    tauytarg = dy1./ds1;
        
    firstbc = (hess(:, :, 1).*(nxtarg.*nxtarg + tauxtarg.*tauxtarg) +  ...
        hess(:, :, 2).*(2*nxtarg.*nytarg + 2*tauxtarg.*tauytarg) + ...
        hess(:, :, 3).*(nytarg.*nytarg + tauytarg.*tauytarg)) + ...
        zk^2*val;
    
    secondbc = (third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
        third(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 4).*(nytarg.*nytarg.*nytarg)) +...
        (third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg+ tauytarg.*tauytarg.*nxtarg) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*nytarg))+...
        zk^2*(grad(:, :, 1).*nxtarg + grad(:, :, 2).*nytarg);

    submat = zeros(2*nt,ns);
    submat(1:2:end,:) = firstbc;
    submat(2:2:end,:) = secondbc;

end

%%% 


% kernels for the clamped plate integral equation
if strcmpi(type, 'clamped_plate')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);
   
   [~, ~, hess, third] = chnk.flex2d.helmdkdiffgreen(zk, src, targ); 
   [~, ~, ~, ~, fourth] = chnk.flex2d.helmdkdiffgreen(zk, src, targ, true);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;
   
   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   rn = rx.*nx + ry.*ny;
   rtau = rx.*taux + ry.*tauy;
   ntargtau = nxtarg.*taux + nytarg.*tauy;

   rntarg = rx.*nxtarg + ry.*nytarg;

   K11 = -(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny)) - ...
       3*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy)); 

   K12 = -(1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

   K21 = -(1/(2*zk^2).*(fourth(:, :, 1).*(nx.*nx.*nx.*nxtarg) + fourth(:, :, 2).*(nx.*nx.*nx.*nytarg + 3*nx.*nx.*ny.*nxtarg) + ...
          fourth(:, :, 3).*(3*nx.*nx.*ny.*nytarg + 3*nx.*ny.*ny.*nxtarg) + fourth(:, :, 4).*(3*nx.*ny.*ny.*nytarg +ny.*ny.*ny.*nxtarg)+...
          fourth(:, :, 5).*(ny.*ny.*ny.*nytarg)) ) - ...
          (3/(2*zk^2).*(fourth(:, :, 1).*(nx.*taux.*taux.*nxtarg)+ fourth(:, :, 2).*(nx.*taux.*taux.*nytarg + 2*nx.*taux.*tauy.*nxtarg + ny.*taux.*taux.*nxtarg) +...
          fourth(:, :, 3).*(2*nx.*taux.*tauy.*nytarg + ny.*taux.*taux.*nytarg + nx.*tauy.*tauy.*nxtarg + 2*ny.*taux.*tauy.*nxtarg) + ...
          fourth(:, :, 4).*(nx.*tauy.*tauy.*nytarg +2*ny.*taux.*tauy.*nytarg + ny.*tauy.*tauy.*nxtarg) +...
          fourth(:, :, 5).*(ny.*tauy.*tauy.*nytarg))) + ...
          1/pi.*(-3*rn.*rntarg./(r2.^2) + 4.*(rn.^3).*rntarg./(r2.^3) + 3*(rn.*rtau.*ntargtau)./ (r2.^2));

   K22 = -(1/(2*zk^2).*(third(:,:, 1).*(nx.*nx.*nxtarg) +third(:, :, 2).*(nx.*nx.*nytarg + 2*nx.*ny.*nxtarg) + third(:, :, 3).*(2*nx.*ny.*nytarg + ny.*ny.*nxtarg)+...
         third(:, :,4).*(ny.*ny.*nytarg))) + ...
         (1/(2*zk^2).*(third(:,:, 1).*(taux.*taux.*nxtarg) +third(:, :, 2).*(taux.*taux.*nytarg + 2*taux.*tauy.*nxtarg) + third(:, :, 3).*(2*taux.*tauy.*nytarg + tauy.*tauy.*nxtarg)+...
         third(:, :,4).*(tauy.*tauy.*nytarg)));

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;
end


% clamped plate kernels for plotting
if strcmpi(type, 'clamped_plate_eval')

    submat = zeros(nt,2*ns);

    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    ds = sqrt(dx.*dx+dy.*dy);

    taux = dx./ds;
    tauy = dy./ds;

    [~, ~, hess, third] = chnk.flex2d.helmdkdiffgreen(zk, src, targ);           % Hankel part

    K1 = -(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny)) - ...
       3*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K2 =  -(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))+...
          (hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy)); % -G_{ny ny}  + G_{tauy tauy}

    submat(:,1:2:end) = K1;
    submat(:,2:2:end) = K2;

end

% clamped plate kernels for plotting
if strcmpi(type, 'free_plate_eval')

    submat = zeros(nt,2*ns);

    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    ds = sqrt(dx.*dx+dy.*dy);

    taux = dx./ds;
    tauy = dy./ds;

    [~, ~, hess, third] = chnk.flex2d.helmdk2diffgreen(zk, src, targ);           % Hankel part

    K1 = -(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny)) - ...
       3*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K2 =  -(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))+...
          (hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy)); % -G_{ny ny}  + G_{tauy tauy}

    submat(:,1:2:end) = K1;
    submat(:,2:2:end) = K2;

end




end