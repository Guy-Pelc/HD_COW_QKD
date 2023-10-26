function [o,p] = gen_constraints(v)
a_ = v.a;
t_ = v.t;
m = v.m;
D = v.D;
ed_ = v.ed;
em_ = v.em;
ep_ = v.ep;


nA = 2^(m*D);
nB = m*D+2;
nE = nB;
aux = nB;
vac = nB-1;

h = helper;
kB = @(i) h.zket(nB,i); %unit vector at Bob
kA = @(i) h.zket(nA,i); %unit vector at Alice

%% generate rho_triv  
syms a %a, occupancy; t, transmittance amplitude
assume(0<a & a<1);
ref_ = sqrt(1-t_^2); %r, reflectance
alpha =@(n) sqrt((1/D).^n.*((D-1)./D).^(D-n));
p0 = exp(-a^2); sp0 = sqrt(p0);
p1 = a^2*p0; sp1 = sqrt(p1);

psi = 0;
for x = 0:nA-1
    n=count(dec2bin(x),'1');
    psi_B = sp0^n*kB(vac);
    for i = 1:(m*D)
        if bitand(x,2^(i-1))~=0
            psi_B = psi_B+kB(i)*sp1*sp0^(n-1);
        end
    end
    psi_B = psi_B+sqrt(1-psi_B'*psi_B)*kB(aux);
    % kA(x+1)
    psi = psi+alpha(n)*kron(kA(x+1),kron(subs(psi_B,a,a_*t_),subs(psi_B,a,a_*ref_)));

end

psi = double(psi);
rho_ABE = psi*psi';

rho_triv = h.reducedRho1(rho_ABE,nA*nB,nE); %trivial matrix at alice,bob

%% generate rhoA_triv

o.rhoA = @(rho) h.reducedRho1(rho,nA,nB);
rhoA_triv = zeros(nA,nA); %
alphas = zeros(nA);
for x = 1:nA
    for y = 1:nA
        nx=count(dec2bin(x-1),'1'); 
        ny=count(dec2bin(y-1),'1');
        n=count(dec2bin(bitxor(x-1,y-1)),'1');
        alphas(x,y) = alpha(nx)*alpha(ny);
        rhoA_triv(x,y) = n;
    end
end
p.rhoA_triv = alphas.*double(subs(sp0,a_)).^rhoA_triv;


%% dataline observables (works for arbitrary d,m)
Md = @(d,ep) (1-ep)^(D*m)*kB(d)*kB(d)'+...
           ep*(1-ep)^(D*m-1)*(kB(vac)*kB(vac)'); 

Md_corr= @(d,ep) kron(h.Ad(d,D,m),Md(d,ep));
Md_err = @(d,ep) kron(h.Ad_not(d,D,m),Md(d,ep));
Md_inc= @(d,ep) kron(eye(nA)-h.Ad(d,D,m)-h.Ad_not(d,D,m),Md(d,ep)); % prob of alice inc when bob measures a single photon at time d

% monitoring line observables
kB_plus = @(c) 1/sqrt(2)*(helper.zket(nB,c)+helper.zket(nB,helper.mod1(c-1,D*m)));
kB_minus = @(c) 1/sqrt(2)*(helper.zket(nB,c)-helper.zket(nB,helper.mod1(c-1,D*m)));

Mc_plus  = @(c,ep,em) (1-ep)^(D*m)*((1-em)*kB_plus(c)*kB_plus(c)'+em*kB_minus(c)*kB_minus(c)') + ep*(1-ep)^(D*m-1)*(kB(vac)*kB(vac)');
Mc_minus = @(c,ep,em) (1-ep)^(D*m)*((1-em)*kB_minus(c)*kB_minus(c)'+em*kB_plus(c)*kB_plus(c)') + ep*(1-ep)^(D*m-1)*(kB(vac)*kB(vac)');
rho_ = h.ch_bitflip(rho_triv,D,m,ed_);

% get probabilities
[p.res_corr,p.res_err] = h.get_dataline_probs(rho_,@(d)Md_err(d,ep_),@(d)Md_corr(d,ep_),false,m,D);
[p.res_plus,p.res_minus] = h.get_monitoring_probs(rho_,@(c)Mc_plus(c,ep_,em_),@(c)Mc_minus(c,ep_,em_),false,m,D);
p.res_inc = h.get_inc_probs(rho_,@(d)Md_inc(d,ep_),m,D); %prob of Alice inc is when bob measures a single photon in a block
% get obseravbles (errorless since eve corrects for errors)
o.Md_corr = @(d) Md_corr(d,0);
o.Md_err = @(d) Md_err(d,0);
o.Mc_plus = @(c) Mc_plus(c,0,0);
o.Mc_minus = @(c) Mc_minus(c,0,0);
o.rho = rho_;


%% generate Xdelta 

kA_block = @(i) helper.zket(2^(D),i);
X = 1/sqrt(D)*fft(eye(D));

Xdelta = 0;
for v=1:m
    for d = 1:D

        d_bob = h.mod1(D+2-d,D); %bob and alice are anticorrelated in kspace!
        kAx= 0;
        kBx = 0;
        % get bob's part and single alice block
        for i = 1:D                    
            kAx = kAx+ X(i,d)*kA_block(1+2^(i-1));
            kBx = kBx + X(i,d_bob)*kB(i+D*(v-1));
        end
        B = kBx*kBx';
        A = 1;
        % wrap with identities for alice's unused blocks
        for v_= 1:m
            if v_==v
                A = kron(kAx*kAx',A);
            else
                A = kron(eye(2^D),A);
            end
        end
        Xdelta_d= kron(A,B);
        Xdelta = Xdelta+ Xdelta_d;   
    end
end
o.Xdelta = Xdelta;

end