%Test jumps for free plate reps
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

figure(1); clf
plot(chnkr)
hold on 
quiver(chnkr)

dens1 = sin(chnkr.r(1,:)) + exp(2i*chnkr.r(2,:));

h = 0.0004;
rng(101);
ind = randi(chnkr.npt);
ind = 4;

targh_lap_out = [];
targh_lap_out.r = chnkr.r(:,ind) + h*[3,2,1,2,2].*chnkr.n(:,ind)+h*[0,0,0,1,-1].*chnkr.d(:,ind);

targh_lap_in = [];
targh_lap_in.r = chnkr.r(:,ind) - h*[1,2,3,2,2].*chnkr.n(:,ind)+h*[0,0,0,1,-1].*chnkr.d(:,ind);

ikern1 = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'free_plate_eval_test_1'); 
ikern2 = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'free_plate_eval_test_2'); 

out1 = chunkerkerneval(chnkr,ikern1,dens1,targh_lap_out);
in1 = chunkerkerneval(chnkr,ikern1,dens1,targh_lap_in);

out_bc1 = (out1(1)-4*out1(2)+out1(3)+out1(4)+out1(5))/h/h + zk^2*out1(2);
in_bc1 = (in1(1)-4*in1(2)+in1(3)+in1(4)+in1(5))/h/h + zk^2*in1(2);

diff1 = abs(out_bc1 - in_bc1 - (dens1(ind)));

assert(abs(diff1) < 1e-2)

nsten = [5,4,3,2,1]; dsten = [0,0,0,0,0];

targh_lap.r = chnkr.r(:,ind) + h*nsten.*chnkr.n(:,ind)+h*chnkr.d(:,ind);
uh_a = chunkerkerneval(chnkr,ikern2,dens1,targh_lap);
out_ntt_a = (-1/12*uh_a(1,:)+2/3*uh_a(2,:)-2/3*uh_a(4,:)+1/12*uh_a(5,:))/h;

targh_lap.r = chnkr.r(:,ind)+ h*nsten.*chnkr.n(:,ind);
uh_b = chunkerkerneval(chnkr,ikern2,dens1,targh_lap);
out_ntt_b = (-1/12*uh_b(1,:)+2/3*uh_b(2,:)-2/3*uh_b(4,:)+1/12*uh_b(5,:))/h;
out_nnn = (0.5*uh_b(1,:)-uh_b(2,:)+uh_b(4,:)-0.5*uh_b(5,:))/h^3;
out_n = out_ntt_b;

targh_lap.r = chnkr.r(:,ind) + h*nsten.*chnkr.n(:,ind)-h*chnkr.d(:,ind);
uh_c = chunkerkerneval(chnkr,ikern2,dens1,targh_lap);
out_ntt_c = (-1/12*uh_c(1,:)+2/3*uh_c(2,:)-2/3*uh_c(4,:)+1/12*uh_c(5,:))/h;

out_ntt = (out_ntt_a - 2*out_ntt_b + out_ntt_c)/h^2;

out_fd = out_nnn + out_ntt + zk^2*out_n;

nsten = -fliplr(nsten);

targh_lap.r = chnkr.r(:,ind) + h*nsten.*chnkr.n(:,ind)+h*chnkr.d(:,ind);
uh_a= chunkerkerneval(chnkr,ikern2,dens1,targh_lap);
in_ntt_a = (-1/12*uh_a(1,:)+2/3*uh_a(2,:)-2/3*uh_a(4,:)+1/12*uh_a(5,:))/h;

targh_lap.r = chnkr.r(:,ind) + h*nsten.*chnkr.n(:,ind);
uh_b= chunkerkerneval(chnkr,ikern2,dens1,targh_lap);
in_ntt_b = (-1/12*uh_b(1,:)+2/3*uh_b(2,:)-2/3*uh_b(4,:)+1/12*uh_b(5,:))/h;
in_nnn = (0.5*uh_b(1,:)-uh_b(2,:)+uh_b(4,:)-0.5*uh_b(5,:))/h^3;
in_n = in_ntt_b;

targh_lap.r = chnkr.r(:,ind) + h*nsten.*chnkr.n(:,ind)-h*chnkr.d(:,ind);
uh_c = chunkerkerneval(chnkr,ikern2,dens1,targh_lap);
in_ntt_c = (-1/12*uh_c(1,:)+2/3*uh_c(2,:)-2/3*uh_c(4,:)+1/12*uh_c(5,:))/h;

in_ntt = (in_ntt_a - 2*in_ntt_b + in_ntt_c)/h^2;

in_fd = in_nnn + in_ntt + zk^2*in_n;

diff2 = abs(out_fd - in_fd - dens1(ind));

assert(abs(diff2) < 1e-2)