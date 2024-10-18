%AUTODIFF_parameterizationTest test the routines for autodifferentiating
%parameterization;
%

clearvars; close all;
seed = 8675309;
rng(seed);
addpaths_loc();

fcurve = @starfish;

nt = 100;
ts = 2*pi*rand(nt,1);

[rt, dt, d2t] = fcurve(ts);

tic;
[r_1, d_1, d2_1] = autodiff_parameterization(fcurve,ts,1);
t1 = toc;
fprintf('%.e second autodiffs per second\n',nt/t1)


err_r1 = norm(r_1-rt,inf);
err_d1 = norm(d_1-dt,inf);
err_d21 = norm(d2_1-d2t,inf);
assert(max([err_r1, err_d1, err_d21]) < 1e-13)

tic;
[r_2, d_2, d2_2] = autodiff_parameterization(fcurve,ts,2);
t2 = toc;
fprintf('%.e first autodiffs per second\n',nt/t2)

err_r2 = norm(r_2-rt,inf);
err_d2 = norm(d_2-dt,inf);
err_d22 = norm(d2_2-d2t,inf);
assert(max([err_r2, err_d2, err_d22]) < 1e-13)

