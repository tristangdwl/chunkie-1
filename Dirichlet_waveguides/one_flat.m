zk = 1;


n = 100;
scheb = cos((2*(1:n)-1)/2/n*pi);
T = cos((0:(n-1)) .*acos(scheb.'));

Bn = -(0:(n-1))/pi;
Bn(1) = -1/pi/log(2);
Bn = Bn(:);


src = [2;0];

targ = [scheb;0*scheb];

f = -chnk.lap2d.green(src,targ);

fc = T\f;
% plot(log10(abs(fc)),'.')

% xtest = -.5;
% ftest = chnk.lap2d.green(src,[xtest;0])
% cos((0:(n-1)) .*acos(xtest))*fc



uc = fc.*Bn;

% us = T*uc./(sqrt(1-scheb.^2).');
% plot(scheb,us)

% %%
% S = zeros(n,n);
% for i = 1:n
%     x = zeros(n,1); x(i) = 1;
%     S(:,i) = T\one_flat_pot(x,targ,xleg,wleg);
% end


%%
Lplot = 1.5;
nplt = 200;
xts = linspace(-Lplot,Lplot,nplt);
[X,Y] = meshgrid(xts,xts);
targs = [X(:).';Y(:).'];

[xleg,wleg] = lege.exps(1000);
% xleg = linspace(-1,1,100).';
% wleg = 2/(length(xleg)-1)*ones(length(xleg),1); wleg([1,end])=wleg([1,end])/2; 
%%
% 
% uscat = one_flat_pot(uc,targs,xleg,wleg);
% uin = chnk.lap2d.green(src,targs);
% 
% %%
% utot = uin+uscat;
% 
% figure(1)
% % h = pcolor(X,Y,reshape(abs(uscat),nplt,nplt)); set(h,'edgecolor','none');
% h = pcolor(X,Y,reshape(abs(utot),nplt,nplt)); set(h,'edgecolor','none');
% colorbar
% 
% %%
% uscat2 = one_flat_pot(uc,targ,xleg,wleg);
% f2 = chnk.lap2d.green(src,targ);
% norm(uscat2(:)-f2(:))
% % figure(2)
% % plot(scheb,[uscat2,f2].','.')
% 

%%

src = [0;.5];
zk = 1;


% n = 100;
% scheb = cos((2*(1:n)-1)/2/n*pi);
% T = cos((0:(n-1)) .*acos(scheb.'));
% 

% 
% 
% [xleg,wleg] = lege.exps(1000);
% targ = [scheb;0*scheb];
% f = -chnk.helm2d.green(zk,src,targ);
% fc = T\f;
% 
% Sk = zeros(n,n);
% for i = 1:n
%     x = zeros(n,1); x(i) = 1;
%     Sk(:,i) = T\one_flat_pot_helm(zk,x,targ,xleg,wleg);
% end
% 
% szk_c = Sk\fc;


nch = 32;
cparams = []; cparams.ta = 0; cparams.tb = pi;
chnkr = chunkerfuncuni(@flatinterface,nch,cparams);
chnkr = sort(chnkr);
fkern_self = @(s,t) 0.25*1i*besselh(0,zk*abs(cos(s.r(1,:))-cos(t.r(1,:).')));
Sk_mat = chunkermat(chnkr,fkern_self);

fkern_self2 = @(s,t) 1/2/pi*log(abs(cos(s.r(1,:))-cos(t.r(1,:).')));
Sk_mat2 = chunkermat(chnkr,fkern_self2);


targ2 = [cos(chnkr.r(1,:));0*chnkr.r(1,:)];
fchnk = -chnk.helm2d.green(zk,src,[cos(chnkr.r(1,:));0*chnkr.r(1,:)]);
% szk = Sk_mat\fchnk;

n = chnkr.npt;
Bn = -(0:(n-1))/pi;
Bn(1) = -1/pi/log(2);
Bn = Bn(:); Bn = 2*pi*Bn;

Tchnk = cos((0:(n-1)) .*acos(cos(chnkr.r(1,:)).'));

B = Tchnk*(Bn.*inv(Tchnk));

Sk_precomp = Tchnk*(Bn.*(Tchnk\Sk_mat));

% szk = Sk_mat\fchnk;

% fchnk2 = -chnk.lap2d.green(src,[cos(chnkr.r(1,:));0*chnkr.r(1,:)]);
% s = Sk_mat2\fchnk2;
% s2 = B*fchnk2;



A = Sk_mat2\Sk_mat;
[cond(Sk_mat),cond(Sk_precomp),cond(A)]

szk = Sk_mat\fchnk;
szk2 = Sk_precomp\(B*fchnk);
szk3 = A\(B*fchnk);


%%
% uzk_c = fc.*Bn;

% uscat = one_flat_pot_helm(zk,szk_c,targs,xleg,wleg);

cdist = @(s,t) sqrt((cos(s.r(1,:))-t.r(1,:).').^2+t.r(2,:).'.^2);
fkern_plot = @(s,t) 0.25*1i*besselh(0,zk*cdist(s,t));
opts=[];opts.forcesmooth = false;
uscat = chunkerkerneval(chnkr,fkern_plot,szk,targs,opts);
uscat3 = chunkerkerneval(chnkr,fkern_plot,szk2,targs,opts);
uscat4 = chunkerkerneval(chnkr,fkern_plot,szk3,targs,opts);
uin = chnk.helm2d.green(zk,src,targs);
utot = uin+uscat;

%%
figure(3)
subplot(1,2,1)
h = pcolor(X,Y,reshape(imag(utot),nplt,nplt)); set(h,'edgecolor','none');
colorbar
subplot(1,2,2)
% h = pcolor(X,Y,reshape(imag(uscat),nplt,nplt)); set(h,'edgecolor','none');
h = pcolor(X,Y,reshape((abs(utot)),nplt,nplt)); set(h,'edgecolor','none');
colorbar

%%
% uscat2 = one_flat_pot_helm(zk,szk_c,targ,xleg,wleg);
% norm(uscat2(:)-f(:))


function u = one_flat_pot(uc,targ,xleg,wleg)
    n = length(uc);
    ntarg = size(targ,2);


    u = zeros(ntarg,1);

    % T = cos((0:(n-1)) .*acos(xleg));
    % sigma_w = T*uc./(sqrt(1-xleg.^2)).*wleg;

    sigma_w2 = cos((0:(n-1))*pi.*(xleg+1)/2)*uc.*wleg;

    for i = 1:ntarg
        % u(i) = sum(log((xleg-targ(1,i)).^2+targ(2,i)^2).*sigma_w);
        u(i) = sum(log((cos(pi*(xleg+1)/2)-targ(1,i)).^2+targ(2,i)^2).*sigma_w2)/2*pi;
    end

    u = u/2;
end


function u = one_flat_pot_helm(zk,uc,targ,xleg,wleg)
    n = length(uc);
    ntarg = size(targ,2);


    u = zeros(ntarg,1);

    % T = cos((0:(n-1)) .*acos(xleg));
    % sigma_w = T*uc./(sqrt(1-xleg.^2)).*wleg;

    sigma_w2 = cos((0:(n-1))*pi.*(xleg+1)/2)*uc.*wleg;

    for i = 1:ntarg
        % u(i) = sum(log((xleg-targ(1,i)).^2+targ(2,i)^2).*sigma_w);
        u(i) = sum(besselh(0,zk*sqrt((cos(pi*(xleg+1)/2)-targ(1,i)).^2+targ(2,i)^2)).*sigma_w2)/2*pi;
    end

    u = .25*1i*u;
end