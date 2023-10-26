function res = optimize_alphas(v,alphas,lossDBs)
    %% calc opt rates
    n_alphas =numel(alphas);
    n_loss = numel(lossDBs);
    h = helper;
    D = v.D;
    m = v.m;

    %%
    rates_per_photon = zeros(n_alphas,n_loss);
    gains= zeros(n_alphas,n_loss);
    inc_probs= zeros(n_alphas,n_loss);
    rel_deltas = zeros(n_alphas,n_loss);
    
    full = cell(n_alphas,n_loss);
     
    parfor j = 1:n_loss
%     parfor j = 1:n_loss
    for i = 1:n_alphas
        v_ =v;
        
        fprintf('%d,%d\n',i,j)
        lossDB =lossDBs(j);
        v_.a= alphas(i);
        v_.t = sqrt(10^(-lossDB/10));
        cur = optimize_once(v_);

        full{i,j} = cur;
        rates_per_photon(i,j) = cur.rpp;
        gains(i,j) = cur.gain;
        inc_probs(i,j) = cur.p.res_inc;
        rel_deltas(i,j) = cur.rel_delta;

    end
    end

    res.full = full;
    res.rates_per_photon = rates_per_photon;
    res.gains = gains;
    res.inc_probs = inc_probs;
    res.rel_deltas = rel_deltas;
end