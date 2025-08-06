%Test jumps for clamped reps
%
% Define an interior scattering problem on a starfish-shaped domain and 
% solve
%

zk = 1.1;

% discretize domain

cparams = [];
cparams.eps = 1.0e-8;
cparams.nover = 0;
cparams.maxchunklen = 4.0/zk; % setting a chunk length helps when the
                              % frequency is known
                              
pref = []; 
pref.k = 16;
narms1 = 4;
amp = 0.25*0.5;
start = tic; 
chnkr = chunkerfunc(@(t) starfish(t,narms1,amp,[0;0],[],3),cparams,pref); 

% figure(1); clf
% plot(chnkr)
% hold on 
% quiver(chnkr)

dens1 = sin(chnkr.r(1,:)) + exp(2i*chnkr.r(2,:));

h = 0.00000005;
rng(101);
ind = randi(chnkr.npt);
eval_pts = [];
eval_pts.r = chnkr.r(:,ind) + h*[-3 -2 -1 1 2 3].*chnkr.n(:,ind);

ikern1 = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_plate_eval_test_1'); 
ikern2 = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_plate_eval_test_2'); 

vals1 = chunkerkerneval(chnkr,ikern1,dens1,eval_pts);
vals2 = chunkerkerneval(chnkr,ikern2,dens1,eval_pts);

dudn1 = ( vals2(6) - vals2(4) ) /2/h;
dudn2 = ( vals2(3) - vals2(1) ) /2/h;

diff1 = abs(vals1(4) - vals1(3) - (dens1(ind)));
diff2 = abs(dudn1 - dudn2 - dens1(ind));

assert(abs(diff1) < 1e-6)
assert(abs(diff2) < 1e-6)