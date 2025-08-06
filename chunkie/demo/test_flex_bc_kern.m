

phi = 0.3;
src = []; src.r = [0;0]; src.n = [cos(phi);sin(phi)]; src.d = [cos(phi+pi/2);sin(phi+pi/2)];

targ = []; targ.r = [1;1]; targ.n = [1;1]/sqrt(2); targ.d = [1;-1]/sqrt(2);
kappa = 0.282842712474619;

% targ = []; targ.r = [1;1]; targ.n = [1;0]; targ.d = [1;0];
% targ = []; targ.r = [1;1]; targ.n = [0;1]; targ.d = [0;-1];

targ.d2 = [0.1;0.3];

zk = 0.9; nu = 0.3;
% nu = 2;
h = 1e-3;

ifree = 0;

%% checking all kernels in the clamped BCs

targh_n = [];
targh_n.r = targ.r + h*[1,0,-1].*targ.n;
targh_n.n = targ.n + 0*[1,0,-1];
targh_n.d = targ.d + 0*[1,0,-1];

ikern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'free_plate_eval'); 
ikern_1 = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'free_to_clamped'); 

u1 = ikern_1(src,targ);

uh= ikern(src,targh_n); u1_fd = (uh(1,:)-uh(3,:))/2/h;

assert(norm(u1(1,:)-uh(2,:))<1e-9)
assert(norm(u1(2,:) - u1_fd)<1e-6)

ikern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_plate_eval'); 
ikern_1 = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_plate'); 

u1 = ikern_1(src,targ);

uh= ikern(src,targh_n); u1_fd = (uh(1,:)-uh(3,:))/2/h;

assert(norm(u1(1,:)-uh(2,:))<1e-9)
assert(norm(u1(2,:) - u1_fd)<1e-6)

ikern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 's'); 
ikern_1 = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_plate_bcs'); 

u1 = ikern_1(src,targ);

uh= ikern(src,targh_n); u1_fd = (uh(1,:)-uh(3,:))/2/h;

assert(norm(u1(1,:)-uh(2,:))<1e-9)
assert(norm(u1(2,:) - u1_fd)<1e-6)


%% checking all kernels in the free plate BC 1

targh_lap = [];
targh_lap.r = targ.r + h*[1,0,-1,0,0].*targ.n+h*[0,0,0,1,-1].*targ.d;
targh_lap.n = targ.n + 0*[1,0,-1,0,0];
targh_lap.d = targ.d + 0*[1,0,-1,0,0];
targh_lap.d2 = targ.d2 + 0*[1,0,-1,0,0];

ikern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_plate_eval'); 
ikern_2 = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_to_free');

u2 = ikern_2(src,targ);

uh= ikern(src,targh_lap); 
u2_fd = (uh(1,:)-4*uh(2,:)+uh(3,:)+uh(4,:)+uh(5,:))/h/h + zk^2*uh(2,:);

assert(norm(u2(1,:) - u2_fd)<1e-6)

ikern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'free_plate_eval'); 
ikern_2 = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'free_plate'); 

u2 = ikern_2(src,targ);

uh= ikern(src,targh_lap); 
u2_fd = (uh(1,:)-4*uh(2,:)+uh(3,:)+uh(4,:)+uh(5,:))/h/h + zk^2*uh(2,:);

assert(norm(u2(1,:) - u2_fd)<1e-6)

ikern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 's'); 
ikern_2 = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'free_plate_bcs', nu); 

u2 = ikern_2(src,targ);

uh= ikern(src,targh_lap); 
u2_fd = (uh(1,:)-4*uh(2,:)+uh(3,:)+uh(4,:)+uh(5,:))/h/h + zk^2*uh(2,:);

assert(norm(u2_fd-u2(1,:))<1e-6)


%% checking all kernels in the free plate BC 2

targh_lap = [];
nsten = [2,1,0,-1,-2]; dsten = [0,0,0,0,0];
targh_lap.r = targ.r + h*nsten.*targ.n+h*dsten.*targ.d;
targh_lap.n = targ.n + 0*nsten;
targh_lap.d = targ.d + 0*nsten;

ikern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_plate_eval'); 
ikern_3 = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'clamped_to_free');

u3 = ikern_3(src,targ);

uh= ikern(src,targh_lap); 

u3_fd_nnn = (0.5*uh(1,:)-uh(2,:)+uh(4,:)-0.5*uh(5,:))/h^3;

targh_lap.r = targ.r + h*nsten.*targ.n+h*targ.d;
uh_a= ikern(src,targh_lap);
u3_fd_ntt_a = (-1/12*uh_a(1,:)+2/3*uh_a(2,:)-2/3*uh_a(4,:)+1/12*uh_a(5,:))/h;

targh_lap.r = targ.r + h*nsten.*targ.n;
uh_b= ikern(src,targh_lap);
u3_fd_ntt_b = (-1/12*uh_b(1,:)+2/3*uh_b(2,:)-2/3*uh_b(4,:)+1/12*uh_b(5,:))/h;

targh_lap.r = targ.r + h*nsten.*targ.n-h*targ.d;
uh_c = ikern(src,targh_lap);
u3_fd_ntt_c = (-1/12*uh_c(1,:)+2/3*uh_c(2,:)-2/3*uh_c(4,:)+1/12*uh_c(5,:))/h;

u3_fd_ntt = (u3_fd_ntt_a - 2*u3_fd_ntt_b + u3_fd_ntt_c)/h^2;

uh_n = u3_fd_ntt_b;

u3_fd = u3_fd_nnn + u3_fd_ntt + zk^2*uh_n;

assert(norm(u3(2,:) - u3_fd)<2e-6)

ikern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'free_plate_eval'); 
ikern_3 = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'free_plate'); 

u3 = ikern_3(src,targ);

uh= ikern(src,targh_lap); 

u3_fd_nnn = (0.5*uh(1,:)-uh(2,:)+uh(4,:)-0.5*uh(5,:))/h^3;

targh_lap.r = targ.r + h*nsten.*targ.n+h*targ.d;
uh_a= ikern(src,targh_lap);
u3_fd_ntt_a = (-1/12*uh_a(1,:)+2/3*uh_a(2,:)-2/3*uh_a(4,:)+1/12*uh_a(5,:))/h;

targh_lap.r = targ.r + h*nsten.*targ.n;
uh_b= ikern(src,targh_lap);
u3_fd_ntt_b = (-1/12*uh_b(1,:)+2/3*uh_b(2,:)-2/3*uh_b(4,:)+1/12*uh_b(5,:))/h;

targh_lap.r = targ.r + h*nsten.*targ.n-h*targ.d;
uh_c = ikern(src,targh_lap);
u3_fd_ntt_c = (-1/12*uh_c(1,:)+2/3*uh_c(2,:)-2/3*uh_c(4,:)+1/12*uh_c(5,:))/h;

u3_fd_ntt = (u3_fd_ntt_a - 2*u3_fd_ntt_b + u3_fd_ntt_c)/h^2;

uh_n = u3_fd_ntt_b;

u3_fd = u3_fd_nnn + u3_fd_ntt + zk^2*uh_n;

assert(norm(u3(2,:) - u3_fd)<8e-5)

ikern = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 's'); 
ikern_3 = @(s,t) chnk.flex2d.kern_alt(zk, s, t, 'free_plate_bcs', nu); 

u3 = ikern_3(src,targ);

uh= ikern(src,targh_lap); 

u3_fd_nnn = (0.5*uh(1,:)-uh(2,:)+uh(4,:)-0.5*uh(5,:))/h^3;

targh_lap.r = targ.r + h*nsten.*targ.n+h*targ.d;
uh_a= ikern(src,targh_lap);
u3_fd_ntt_a = (-1/12*uh_a(1,:)+2/3*uh_a(2,:)-2/3*uh_a(4,:)+1/12*uh_a(5,:))/h;

targh_lap.r = targ.r + h*nsten.*targ.n;
uh_b= ikern(src,targh_lap);
u3_fd_ntt_b = (-1/12*uh_b(1,:)+2/3*uh_b(2,:)-2/3*uh_b(4,:)+1/12*uh_b(5,:))/h;

targh_lap.r = targ.r + h*nsten.*targ.n-h*targ.d;
uh_c = ikern(src,targh_lap);
u3_fd_ntt_c = (-1/12*uh_c(1,:)+2/3*uh_c(2,:)-2/3*uh_c(4,:)+1/12*uh_c(5,:))/h;

u3_fd_ntt = (u3_fd_ntt_a - 2*u3_fd_ntt_b + u3_fd_ntt_c)/h^2;

uh_n = u3_fd_ntt_b;

u3_fd = u3_fd_nnn + u3_fd_ntt + zk^2*uh_n;

assert(norm(u3(2,:) - u3_fd)<1e-6)