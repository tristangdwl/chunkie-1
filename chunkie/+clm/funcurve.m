function [r,d,d2] = funcurve(t,icurve,cpars,icase)
%%funcurve
% return position, first and second derivatives of four curve segments
% in the embedded eye layered media problem.
%
% Inputs:
% t - paramter values to evaluate these quantities
%
% Outputs:
% r - coordinates
% d - first derivatives w.r.t. t
% d2 - second derivatives w.r.t. t

if icase == 2
  [r,d,d2] = clm.complexx(t,cpars.L,cpars.c1,cpars.c2);
elseif icase == 4
  [r,d,d2] = clm.funcurve4(t,icurve,cpars);  
elseif icase == 6
  [r,d,d2] = clm.funcurve6(t,icurve,cpars);
end

end