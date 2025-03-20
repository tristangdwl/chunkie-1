function val = comp_obs(src, targ, sys, kernstr, opts)


if nargin < 5
    opts = []; opts.eps = 1e-12;
end  

opts.isp = 2;

if strcmp(kernstr,'s')
    % f_inc  = @(s,t) wave_kern.s(s, t, sys, opts);
    idp = 0;
    igrad = 0;
elseif strcmp(kernstr,'d')
    % f_inc  = @(s,t) wave_kern.d(s, t, sys, opts);
    idp = 1;
    igrad = 0;
elseif strcmp(kernstr,'sp')
    % f_inc  = @(s,t) wave_kern.sp(s, t, sys, opts);
    idp = 0;
    igrad = 1;
    kernstr = 'sprime';
elseif strcmp(kernstr,'dp')
    % f_inc  = @(s,t) wave_kern.dp(s, t, sys, opts);
    idp = 1;
    igrad = 1;
    kernstr = 'dprime';
end

iempty = 0;
if isfield(sys,'iempty'), iempty= sys.iempty; end

if iempty
    val = 0;
else
    val = compress.apply_comp_far(sys.comp_data.f2f,sys,src,targ,idp,igrad);
end

val = val + chnk.helm2d.kern(sys.zk1,src,targ,kernstr);




end