function [r, d, d2] = autodiff_parameterization(fcurve,ts,nout)

if nout == 1
    auto_fcurve = @(t) auto_der(fcurve,t);
elseif nout == 2
    auto_fcurve = @(t) auto_sec_der(fcurve,t);
else
    auto_fcurve = fcurve;
end

nt = length(ts);
r = zeros(2,nt);
d = zeros(2,nt);
d2 = zeros(2,nt);
for i = 1:nt
    t = dlarray(ts(i));
    [r(:,i), d(:,i), d2(:,i)] = dlfeval(auto_fcurve,t);
end

end

function [r, d, d2] = auto_der(fcurve,t)
    r = fcurve(t);
    d = dlarray(zeros(2,1));
    d(1) = dlgradient(r(1),t,EnableHigherDerivatives=true,RetainData=true);
    d(2) = dlgradient(r(2),t,EnableHigherDerivatives=true,RetainData=true);
    d2 = dlarray(zeros(2,1));
    d2(1) = dlgradient(d(1),t,RetainData=true);
    d2(2) = dlgradient(d(2),t);
end

function [r, d, d2] = auto_sec_der(fcurve,t)
    [r,d] = fcurve(t);
    d2 = dlarray(zeros(2,1));
    d2(1) = dlgradient(d(1),t,RetainData=true);
    d2(2) = dlgradient(d(2),t);
end
