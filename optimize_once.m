function res = optimize_once(v,alpha,lossDB)
arguments
    v
    alpha = []
    lossDB = []
end
    h = helper;
    D = v.D;
    m = v.m;
if ~isempty(lossDB)
    v.t = sqrt(10^(-lossDB/10));
end
if ~isempty(alpha)
    v.a = alpha;
end


    options = sdpsettings('solver','mosek','verbose',0, 'showprogress', 0, 'cachesolvers', 1,'dualize',1,'debug',0);
    dim_rho = 2^(D*m)*(D*m+2);
    rho = sdpvar(dim_rho, dim_rho,'hermitian','complex');

    [o,p] = gen_constraints(v);
    
    F = [rho >= 0,                          % fully positive density matrix
        trace(rho) == 1,                    % density matrix trace constraint (1)
        ];
    %set params
    %alice constraints
    rhoA_opt = o.rhoA(rho);    
    F = [F,rhoA_opt-p.rhoA_triv==0];
    
    %dataline constraints
    [res_corr_opt,res_err_opt] = h.get_dataline_probs(rho,o.Md_err,o.Md_corr,true,m,D);
    F = [F,res_err_opt-p.res_err==0,res_corr_opt-p.res_corr==0];
       
    %monitoring line constraints    
    [res_plus_opt,res_minus_opt] = h.get_monitoring_probs(rho,o.Mc_plus,o.Mc_minus,true,m,D);    
    F = [F,res_plus_opt-p.res_plus==0, res_minus_opt-p.res_minus == 0];
    
    
    %calc results
    gain = sum(p.res_corr(:)+p.res_err(:));
    objective = h.pr(o.Xdelta,rho);
     sol = optimize(F,objective,options);
     rho_opt = value(rho);
         
    delta = real(gain-h.pr(o.Xdelta,rho_opt));
    
        
    rel_err = sum(p.res_err(:))/gain;
    rel_delta = delta/gain;   
    
    if max(rel_err,rel_delta)>(1-1/D) %we require h2(delta,D) is monotonic, 
        rate_per_photon = 0;
    else
        rate_per_photon = max(log2(D)-h.h2(rel_err,D)-h.h2(rel_delta,D),0);
    end

    res.v = v;
    res.rel_delta = rel_delta;
    res.rel_err = rel_err;
    res.gain = gain;
    res.rpp = rate_per_photon;
    res.rho = rho_opt;
    res.o =o;
    res.p = p;
end