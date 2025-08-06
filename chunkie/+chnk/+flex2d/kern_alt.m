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
    
    [val, grad, hess, third] = chnk.flex2d.helmdkdiffgreen(zk, src, targ);

    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);

            
    firstbc = hess(:, :, 1) + hess(:, :, 3) + zk^2*val;
    
    secondbc = nxtarg.*(third(:,:,1) + third(:,:,3)) + ...
        nytarg.*(third(:,:,2) + third(:,:,4)) + ...
        zk^2*(grad(:, :, 1).*nxtarg + grad(:, :, 2).*nytarg);

    submat = zeros(2*nt,ns);
    submat(1:2:end,:) = firstbc;
    submat(2:2:end,:) = secondbc;

end

%%% INTEGRAL REPRESENTATIONS

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

% free plate kernels for plotting
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

    [~, ~, hess, third] = chnk.flex2d.helmdkdkdiffgreen(zk, src, targ);           % Hankel part

    K1 = -(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny)) - ...
       3*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K2 =  -(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))+...
          (hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy)); % -G_{ny ny}  + G_{tauy tauy}

    submat(:,1:2:end) = K1;
    submat(:,2:2:end) = K2;

end

%%% INTEGRAL EQUATIONS - DIAGONAL BLOCKS

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
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy)); % G_{ny ny ny} + 3G_{ny tauy tauy}

   K12 = -((hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny)))+...
          ((hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

   K21 = -(fourth(:, :, 1).*(nx.*nx.*nx.*nxtarg) + fourth(:, :, 2).*(nx.*nx.*nx.*nytarg + 3*nx.*nx.*ny.*nxtarg) + ...
          fourth(:, :, 3).*(3*nx.*nx.*ny.*nytarg + 3*nx.*ny.*ny.*nxtarg) + fourth(:, :, 4).*(3*nx.*ny.*ny.*nytarg +ny.*ny.*ny.*nxtarg)+...
          fourth(:, :, 5).*(ny.*ny.*ny.*nytarg)) - ...
          3*(fourth(:, :, 1).*(nx.*taux.*taux.*nxtarg)+ fourth(:, :, 2).*(nx.*taux.*taux.*nytarg + 2*nx.*taux.*tauy.*nxtarg + ny.*taux.*taux.*nxtarg) +...
          fourth(:, :, 3).*(2*nx.*taux.*tauy.*nytarg + ny.*taux.*taux.*nytarg + nx.*tauy.*tauy.*nxtarg + 2*ny.*taux.*tauy.*nxtarg) + ...
          fourth(:, :, 4).*(nx.*tauy.*tauy.*nytarg +2*ny.*taux.*tauy.*nytarg + ny.*tauy.*tauy.*nxtarg) +...
          fourth(:, :, 5).*(ny.*tauy.*tauy.*nytarg)) + ...
          -1/pi.*(-3*rn.*rntarg./(r2.^2) + 4.*(rn.^3).*rntarg./(r2.^3) + 3*(rn.*rtau.*ntargtau)./ (r2.^2)); % G_{nx ny ny ny} + 3G_{nx ny tauy tauy}

   K22 = -(third(:,:, 1).*(nx.*nx.*nxtarg) +third(:, :, 2).*(nx.*nx.*nytarg + 2*nx.*ny.*nxtarg) + third(:, :, 3).*(2*nx.*ny.*nytarg + ny.*ny.*nxtarg)+...
         third(:, :,4).*(ny.*ny.*nytarg)) + ...
         (third(:,:, 1).*(taux.*taux.*nxtarg) +third(:, :, 2).*(taux.*taux.*nytarg + 2*taux.*tauy.*nxtarg) + third(:, :, 3).*(2*taux.*tauy.*nytarg + tauy.*tauy.*nxtarg)+...
         third(:, :,4).*(tauy.*tauy.*nytarg)); % -G_{nx ny ny}  + G_{nx tauy tauy}

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;
end


% kernels for the free plate integral equation % SAME AS CLAMPED PLATE
if strcmpi(type, 'free_plate')
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
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy)); % G_{ny ny ny} + 3G_{ny tauy tauy}

   K12 = -((hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny)))+...
          ((hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

   K21 = -(fourth(:, :, 1).*(nx.*nx.*nx.*nxtarg) + fourth(:, :, 2).*(nx.*nx.*nx.*nytarg + 3*nx.*nx.*ny.*nxtarg) + ...
          fourth(:, :, 3).*(3*nx.*nx.*ny.*nytarg + 3*nx.*ny.*ny.*nxtarg) + fourth(:, :, 4).*(3*nx.*ny.*ny.*nytarg +ny.*ny.*ny.*nxtarg)+...
          fourth(:, :, 5).*(ny.*ny.*ny.*nytarg)) - ...
          3*(fourth(:, :, 1).*(nx.*taux.*taux.*nxtarg)+ fourth(:, :, 2).*(nx.*taux.*taux.*nytarg + 2*nx.*taux.*tauy.*nxtarg + ny.*taux.*taux.*nxtarg) +...
          fourth(:, :, 3).*(2*nx.*taux.*tauy.*nytarg + ny.*taux.*taux.*nytarg + nx.*tauy.*tauy.*nxtarg + 2*ny.*taux.*tauy.*nxtarg) + ...
          fourth(:, :, 4).*(nx.*tauy.*tauy.*nytarg +2*ny.*taux.*tauy.*nytarg + ny.*tauy.*tauy.*nxtarg) +...
          fourth(:, :, 5).*(ny.*tauy.*tauy.*nytarg)) + ...
          -1/pi.*(-3*rn.*rntarg./(r2.^2) + 4.*(rn.^3).*rntarg./(r2.^3) + 3*(rn.*rtau.*ntargtau)./ (r2.^2)); % G_{nx ny ny ny} + 3G_{nx ny tauy tauy}

   K22 = -(third(:,:, 1).*(nx.*nx.*nxtarg) +third(:, :, 2).*(nx.*nx.*nytarg + 2*nx.*ny.*nxtarg) + third(:, :, 3).*(2*nx.*ny.*nytarg + ny.*ny.*nxtarg)+...
         third(:, :,4).*(ny.*ny.*nytarg)) + ...
         (third(:,:, 1).*(taux.*taux.*nxtarg) +third(:, :, 2).*(taux.*taux.*nytarg + 2*taux.*tauy.*nxtarg) + third(:, :, 3).*(2*taux.*tauy.*nytarg + tauy.*tauy.*nxtarg)+...
         third(:, :,4).*(tauy.*tauy.*nytarg)); % -G_{nx ny ny}  + G_{nx tauy tauy}

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;
end


%%% INTEGRAL EQUATIONS - OFF DIAGONAL BLOCKS

% free plate rep in clamped plate BCs
if strcmpi(type, 'free_to_clamped')           
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;

   targnorm = targinfo.n;

   [~,~,hess,third,fourth] = chnk.flex2d.helmdkdkdiffgreen(zk,src,targ); 

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds; 
   tauy = dy./ds;

    K11 = -(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny)) - ...
       3*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K12 =  -(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))+...
          (hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy)); % -G_{ny ny}  + G_{tauy tauy}

    K1x = -(fourth(:, :, 1).*(nx.*nx.*nx) + fourth(:, :, 2).*(3*nx.*nx.*ny) +...
       fourth(:, :, 3).*(3*nx.*ny.*ny) + fourth(:, :, 4).*(ny.*ny.*ny)) - ...
       3*(fourth(:, :, 1).*(nx.*taux.*taux) + fourth(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       fourth(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + fourth(:, :, 4).*(ny.*tauy.*tauy));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    ioff = 1;
    K1y = -(fourth(:, :, 1+ioff).*(nx.*nx.*nx) + fourth(:, :, 2+ioff).*(3*nx.*nx.*ny) +...
       fourth(:, :, 3+ioff).*(3*nx.*ny.*ny) + fourth(:, :, 4+ioff).*(ny.*ny.*ny)) - ...
       3*(fourth(:, :, 1+ioff).*(nx.*taux.*taux) + fourth(:, :, 2+ioff).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       fourth(:, :, 3+ioff).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + fourth(:, :, 4+ioff).*(ny.*tauy.*tauy));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K2x =  -(third(:, :, 1).*(nx.*nx) + third(:, :, 2).*(2*nx.*ny) + third(:, :, 3).*(ny.*ny))+...
          (third(:, :, 1).*(taux.*taux) + third(:, :, 2).*(2*taux.*tauy) + third(:, :, 3).*(tauy.*tauy)); % -G_{ny ny}  + G_{tauy tauy}

    ioff = 1;
    K2y =  -(third(:, :, 1+ioff).*(nx.*nx) + third(:, :, 2+ioff).*(2*nx.*ny) + third(:, :, 3+ioff).*(ny.*ny))+...
          (third(:, :, 1+ioff).*(taux.*taux) + third(:, :, 2+ioff).*(2*taux.*tauy) + third(:, :, 3+ioff).*(tauy.*tauy)); % -G_{ny ny}  + G_{tauy tauy}

    K21 = K1x.*nxtarg + K1y.*nytarg;

    K22 = K2x.*nxtarg + K2y.*nytarg;

  submat = zeros(2*nt,2*ns);

  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;

end


% clamped plate rep in free plate BCs
if strcmpi(type, 'clamped_to_free')           
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;

   targnorm = targinfo.n;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds; 
   tauy = dy./ds;


    [~, ~, hess, third,fourth,fifth] = chnk.flex2d.helmdkdiffgreen(zk, src, targ);           % Hankel part
    sixth = chnk.flex2d.green_helmsq_sixth(zk,src,targ);

    K1 = -(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny)) - ...
       3*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K1xx = -(fifth(:, :, 1).*(nx.*nx.*nx) + fifth(:, :, 2).*(3*nx.*nx.*ny) +...
       fifth(:, :, 3).*(3*nx.*ny.*ny) + fifth(:, :, 4).*(ny.*ny.*ny)) - ...
       3*(fifth(:, :, 1).*(nx.*taux.*taux) + fifth(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       fifth(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + fifth(:, :, 4).*(ny.*tauy.*tauy));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    ioff = 2;
    K1yy = -(fifth(:, :, 1+ioff).*(nx.*nx.*nx) + fifth(:, :, 2+ioff).*(3*nx.*nx.*ny) +...
       fifth(:, :, 3+ioff).*(3*nx.*ny.*ny) + fifth(:, :, 4+ioff).*(ny.*ny.*ny)) - ...
       3*(fifth(:, :, 1+ioff).*(nx.*taux.*taux) + fifth(:, :, 2+ioff).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       fifth(:, :, 3+ioff).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + fifth(:, :, 4+ioff).*(ny.*tauy.*tauy));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K11 = K1xx + K1yy + zk^2*K1;

    K1xxx = -(sixth(:, :, 1).*(nx.*nx.*nx) + sixth(:, :, 2).*(3*nx.*nx.*ny) +...
       sixth(:, :, 3).*(3*nx.*ny.*ny) + sixth(:, :, 4).*(ny.*ny.*ny)) - ...
       3*(sixth(:, :, 1).*(nx.*taux.*taux) + sixth(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       sixth(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + sixth(:, :, 4).*(ny.*tauy.*tauy));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    ioff = 1;
    K1xxy = -(sixth(:, :, 1+ioff).*(nx.*nx.*nx) + sixth(:, :, 2+ioff).*(3*nx.*nx.*ny) +...
       sixth(:, :, 3+ioff).*(3*nx.*ny.*ny) + sixth(:, :, 4+ioff).*(ny.*ny.*ny)) - ...
       3*(sixth(:, :, 1+ioff).*(nx.*taux.*taux) + sixth(:, :, 2+ioff).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       sixth(:, :, 3+ioff).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + sixth(:, :, 4+ioff).*(ny.*tauy.*tauy));  % G_{ny ny ny} + 3G_{ny tauy tauy}
    ioff = 2;
    K1xyy = -(sixth(:, :, 1+ioff).*(nx.*nx.*nx) + sixth(:, :, 2+ioff).*(3*nx.*nx.*ny) +...
       sixth(:, :, 3+ioff).*(3*nx.*ny.*ny) + sixth(:, :, 4+ioff).*(ny.*ny.*ny)) - ...
       3*(sixth(:, :, 1+ioff).*(nx.*taux.*taux) + sixth(:, :, 2+ioff).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       sixth(:, :, 3+ioff).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + sixth(:, :, 4+ioff).*(ny.*tauy.*tauy));  % G_{ny ny ny} + 3G_{ny tauy tauy}
    ioff = 3;
    K1yyy = -(sixth(:, :, 1+ioff).*(nx.*nx.*nx) + sixth(:, :, 2+ioff).*(3*nx.*nx.*ny) +...
       sixth(:, :, 3+ioff).*(3*nx.*ny.*ny) + sixth(:, :, 4+ioff).*(ny.*ny.*ny)) - ...
       3*(sixth(:, :, 1+ioff).*(nx.*taux.*taux) + sixth(:, :, 2+ioff).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       sixth(:, :, 3+ioff).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + sixth(:, :, 4+ioff).*(ny.*tauy.*tauy));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K1x = -(fourth(:, :, 1).*(nx.*nx.*nx) + fourth(:, :, 2).*(3*nx.*nx.*ny) +...
       fourth(:, :, 3).*(3*nx.*ny.*ny) + fourth(:, :, 4).*(ny.*ny.*ny)) - ...
       3*(fourth(:, :, 1).*(nx.*taux.*taux) + fourth(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       fourth(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + fourth(:, :, 4).*(ny.*tauy.*tauy));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    ioff = 1;
    K1y = -(fourth(:, :, 1+ioff).*(nx.*nx.*nx) + fourth(:, :, 2+ioff).*(3*nx.*nx.*ny) +...
       fourth(:, :, 3+ioff).*(3*nx.*ny.*ny) + fourth(:, :, 4+ioff).*(ny.*ny.*ny)) - ...
       3*(fourth(:, :, 1+ioff).*(nx.*taux.*taux) + fourth(:, :, 2+ioff).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       fourth(:, :, 3+ioff).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + fourth(:, :, 4+ioff).*(ny.*tauy.*tauy));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K21 = nxtarg.*(K1xxx+K1xyy + zk^2*K1x) + nytarg.*(K1xxy + K1yyy + zk^2*K1y); 

    K2 =  -(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))+...
          (hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy)); % -G_{ny ny}  + G_{tauy tauy}

    K2xx =  -(fourth(:, :, 1).*(nx.*nx) + fourth(:, :, 2).*(2*nx.*ny) + fourth(:, :, 3).*(ny.*ny))+...
          (fourth(:, :, 1).*(taux.*taux) + fourth(:, :, 2).*(2*taux.*tauy) + fourth(:, :, 3).*(tauy.*tauy)); % -G_{ny ny}  + G_{tauy tauy}

    ioff = 2;
    K2yy =  -(fourth(:, :, 1 + ioff).*(nx.*nx) + fourth(:, :, 2+ioff).*(2*nx.*ny) + fourth(:, :, 3+ioff).*(ny.*ny))+...
          (fourth(:, :, 1+ioff).*(taux.*taux) + fourth(:, :, 2+ioff).*(2*taux.*tauy) + fourth(:, :, 3+ioff).*(tauy.*tauy)); % -G_{ny ny}  + G_{tauy tauy}

    K12 = K2xx + K2yy + zk^2*K2;

    K2x =  -(third(:, :, 1).*(nx.*nx) + third(:, :, 2).*(2*nx.*ny) + third(:, :, 3).*(ny.*ny))+...
          (third(:, :, 1).*(taux.*taux) + third(:, :, 2).*(2*taux.*tauy) + third(:, :, 3).*(tauy.*tauy)); % -G_{ny ny}  + G_{tauy tauy}

    ioff = 1;
    K2y =  -(third(:, :, 1+ioff).*(nx.*nx) + third(:, :, 2+ioff).*(2*nx.*ny) + third(:, :, 3+ioff).*(ny.*ny))+...
          (third(:, :, 1+ioff).*(taux.*taux) + third(:, :, 2+ioff).*(2*taux.*tauy) + third(:, :, 3+ioff).*(tauy.*tauy)); % -G_{ny ny}  + G_{tauy tauy}


    K2xxx =  -(fifth(:, :, 1).*(nx.*nx) + fifth(:, :, 2).*(2*nx.*ny) + fifth(:, :, 3).*(ny.*ny))+...
          (fifth(:, :, 1).*(taux.*taux) + fifth(:, :, 2).*(2*taux.*tauy) + fifth(:, :, 3).*(tauy.*tauy)); % -G_{ny ny}  + G_{tauy tauy}

    ioff = 1;
    K2xxy =  -(fifth(:, :, 1 + ioff).*(nx.*nx) + fifth(:, :, 2+ioff).*(2*nx.*ny) + fifth(:, :, 3+ioff).*(ny.*ny))+...
          (fifth(:, :, 1+ioff).*(taux.*taux) + fifth(:, :, 2+ioff).*(2*taux.*tauy) + fifth(:, :, 3+ioff).*(tauy.*tauy)); % -G_{ny ny}  + G_{tauy tauy}

    ioff = 2;
    K2xyy =  -(fifth(:, :, 1 + ioff).*(nx.*nx) + fifth(:, :, 2+ioff).*(2*nx.*ny) + fifth(:, :, 3+ioff).*(ny.*ny))+...
          (fifth(:, :, 1+ioff).*(taux.*taux) + fifth(:, :, 2+ioff).*(2*taux.*tauy) + fifth(:, :, 3+ioff).*(tauy.*tauy)); % -G_{ny ny}  + G_{tauy tauy}
    
    ioff = 3;
    K2yyy =  -(fifth(:, :, 1 + ioff).*(nx.*nx) + fifth(:, :, 2+ioff).*(2*nx.*ny) + fifth(:, :, 3+ioff).*(ny.*ny))+...
          (fifth(:, :, 1+ioff).*(taux.*taux) + fifth(:, :, 2+ioff).*(2*taux.*tauy) + fifth(:, :, 3+ioff).*(tauy.*tauy)); % -G_{ny ny}  + G_{tauy tauy}

    K22 = nxtarg.*(K2xxx+K2xyy + zk^2*K2x) + nytarg.*(K2xxy + K2yyy + zk^2*K2y); 

  submat = zeros(2*nt,2*ns);

  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;

end




end