%% operator and param descriptions

% alphas occupancy
% ed prob of bit flip, 0.5% in paper
% em misalignment error, 0.25% in paper
% ep dark count, 5e-8 in paper

% delta- prob of phase error, maximized by optimization
% Xdelta- trace(Xdelta*rho) is the objective of optimization
% h, helper module
% Ki(j) is the jth basis vector in i's space (alice or bob)
% m number of d-bits
% rho_ trivial matrix after quantum channel
% rhoA - reduced density matrix function at Alice
% rho_triv- trivial matrix with beamsplitter loss, without bit or phase errs.
% vars- [ep, ed, em a t]
% vals- [values of vars]
% px_plus probabilities of monitoring plus detector where (a,b) index is
% ath state at alice sent and detection at time b at bob's detector. 
% px_minus probabilities of monitoring minus detector

%% set opt parameters 
m=1
D=2
err = 0.5
name = sprintf("m%dD%d_%0.2f_err",m,D,err);

addpath(genpath('../YALMIP-master')) %yalmip dir
addpath(genpath('../mosek')) %mosek dir


alphas = 1*sqrt(logspace(0,-3,10));
lossDBs = 1:12;

clear v;
v.ed = err*1e-2;
v.ep = err*1e-7;
v.em = err*0.5e-2;
v.D = D;
v.m = m;


%% run optimization
tic
res = optimize_alphas(v,alphas,lossDBs);
toc
beep
% rate_opts (i,j): ith alpha and jth loss. rate per time slot
% rate_trivs:
% rel_deltas: optimized_phase_err/gain
% rel_errs: pd_err/gainsol
% save results
res.alphas = alphas;
res.lossDBs = lossDBs;
res.vals = [v.ep,v.ed,v.em,0,0];
res.vars = ['ep','ed','em','a','t'];
res.D = D;
res.m = m;

%% plot results using h.rps for rate per second, and h.rpb for rate per block
h = helper;
figure(5); 
o = h.rps(res,'ps',1); plot(res.lossDBs,o); 
grid on; set(gca,'YScale','log');
title('rate per second')

figure(6);
o = h.rpb(res); plot(res.lossDBs,o); 
grid on; set(gca,'YScale','log');
 title('rate per block')


%% save res object
% save(sprintf("res_%s_%s.mat",name,datetime('now','Format','yyyy.MM.dd.HHmm')),"res")
