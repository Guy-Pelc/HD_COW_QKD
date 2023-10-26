classdef helper
    methods(Static)
        function res = floor1(a,n)
        res = floor((a-1)/n)+1;
        end
        %% mod(a,n) starting at 1 and ending at n
        function b = mod1(a,n)
        b = mod(a-1,n)+1;
        end
        function y = reducedRho1(op, dim1, dim2)
        rho1= 0;
        for i=1:dim2
            zj = kron(eye(dim1),helper.zket(dim2,i)); 
            rho1 = rho1 + zj' * op * zj;
        end
        y = rho1;
        end
        function [y] = reducedRho2(op, dim1, dim2)
        rho2= 0;
        for i=1:dim1
            zj = kron(helper.zket(dim1,i),eye(dim2)); 
            rho2 = rho2 + zj' * op * zj;
        end
        y = rho2;
        end
        %% FUNCTION NAME: zket
        % It outpus a ket vector $\ket{j}$ in the computational basis (z basis)
        % when the index = j, counting from 1 to dimension = dim.
        function outputVec = zket(dim, index)
            idMatrix=eye(dim);
            outputVec =idMatrix(:,index);
        end

    %%
%         function rho2 = ch_phaseflip(rho,nA,nB,em)
%             h = helper;
% %             if nA~=4; error('m not supported');end
%         Mphase0= sqrt(1-em)*eye(nA*nB);
%         rho2 = rho;
%         m = log2(nA);
% 
%         for i = 1:2*m
%             rho2 = Mphase0*rho2*Mphase0'+h.Mphase_c(i,nA,nB,em)*rho2*h.Mphase_c(i,nA,nB,em)'; % flip phase of pulses i,i+2
%         end
% %         rho2 = Mphase0*rho*Mphase0'+h.Mphase_c(1,nA,nB,em)*rho*h.Mphase_c(1,nA,nB,em)'; %flip phase pulses 1,3
% %         rho2 = Mphase0*rho2*Mphase0'+h.Mphase_c(2,nA,nB,em)*rho2*h.Mphase_c(2,nA,nB,em)'; % flip phase pulses 2,4
%         end
%         % pauli like operator sigma_z over pulses c,c+2 mod(2*m)
%         function res = Mphase_c(c,nA,nB,em)        
%             h = helper;
%             m = log2(nA);
%             res = eye(nB);
%             res(h.mod1(c+2,2*m),h.mod1(c+2,2*m)) = -1;
%             res = sqrt(em)*kron(eye(nA),res);
%         end
   
        function [rho2] = ch_bitflip(rho,D,m,ed)
            nA = 2^(m*D);
            nB = m*D+2;            
            h = helper;            
            kB = @(i) h.zket(nB,i);

            next = @(d) h.mod1(d+1,D*m);
            Md = @(d) kron(eye(nA), eye(nB)- kB(d)*kB(d)'-kB(next(d))*kB(next(d))'+kB(d)*kB(next(d))'+kB(next(d))*kB(d)');
            rho2 = rho*(1-0.5*D*m*ed);
            for i = 1:D*m-1 %change to D*m-1 to remove cyclicity
                rho2 = rho2+0.5*ed*Md(i)*rho*Md(i)';
            end                                     
        end
        %%
%         function [res_corr,res_err] = get_dataline_probs(rho,Md_err,Md_corr)
%             res_err = helper.pr(Md_err(1),rho);
%             res_corr = helper.pr(Md_corr(1),rho);
%         end
        
        function [res_corr,res_err] = get_dataline_probs(rho,Md_err,Md_corr,is_opt,m,D)
            h = helper;           
            if is_opt
            res_corr = sdpvar(D*m,1,'full');
            res_err= sdpvar(D*m,1,'full');
            else
            res_corr  = zeros([D*m,1]);
            res_err= zeros([D*m,1]);
            end
            for d = 1:D*m
                res_err(d) = h.pr(Md_err(d),rho);
                res_corr(d) = h.pr(Md_corr(d),rho);
            end
        end

        function res = get_inc_probs(rho,Md_inc,m,D)
            res = 0;
            for d = 1:m*D
                res = res+helper.pr(Md_inc(d),rho);
            end            
        end

        function [res_plus,res_minus]=get_monitoring_probs(rho,Mc_plus,Mc_minus,is_opt,m,D) %% see 20220720 derivation
            h = helper;
            nA = 2^(D*m);
            kA = @(i) h.zket(nA,i); %unit vector at Alice

            if is_opt
                res_plus = sdpvar(4,D*m,'full');
                res_minus = sdpvar(4,D*m,'full');
            else
                res_plus  = zeros([4,D*m]);
                res_minus = zeros([4,D*m]);
            end
            for a = 1:4
                    a0 = mod(a-1,2);
                    a1 = mod(bitshift(a-1,-1),2);
                for c = 1:D*m 
                            mask = 2^mod(c-1,D*m)+2^mod(c-2,D*m);
                            target =  a1*2^mod(c-1,D*m)+a0*2^mod(c-2,D*m);
                            Ax = 0;
                            for x= 0:nA-1
                                r = bitand(x,mask);
                                r=bitxor(r,target);
                                if r==0
%                                     dec2bin(x)
                                    Ax = Ax+kA(x+1)*kA(x+1)';
                                end
                            end
%                             Ax
                        res_plus(a,c) = h.pr(kron(Ax,Mc_plus(c)),rho);
                        res_minus(a,c) = h.pr(kron(Ax,Mc_minus(c)),rho);
                end
            end
             res_plus = res_plus(:,2:end);%remove cyclicity
             res_minus = res_minus(:,2:end);%remove cyclicity
        end

        %%
        % get prob of measuring A at state rho
        function res = pr(A,rho)
            rho_t = rho.';
            res = A(:).'*rho_t(:);
%             res = trace(A*rho);
        end
        %%

%%  Ad(d), Ad_not, given bob measures at time d, project onto alice subspace possible states that lead to correct/incorrect key contribution
% see derivation 20220719
%works for any dimension D, any block size m
function res = Ad(d,D,m)

b = 1+floor((d-1)/D); % get block number
res = 1;
for b_ = 1:m
    if b_ == b
        res = kron(helper.Ad_(helper.mod1(d,D),D),res);
    else
        res = kron(eye(2^(D)),res);
    end
end
end

function res = Ad_not(d,D,m)
    b = 1+floor((d-1)/D); % get block number
    res = 1;
    for b_ = 1:m
        if b_ == b
            res = kron(helper.Ad_not_(helper.mod1(d,D),D),res);
        else
            res = kron(eye(2^(D)),res);
        end
    end
end
function res = Ad_(d,D)
%     if m>1; error('m not supported');end 
    kA = @(i) helper.zket(2^(D),i);
    idx = 1+2^(d-1);
    res = kA(idx)*kA(idx)';
end



function res = Ad_not_(d,D)    
    kA = @(i) helper.zket(2^D,i);
    res = 0;
    for d_ = 1:D
    if d_ == d; continue;end
    idx = 1+2^(d_-1);
    res = res + kA(idx)*kA(idx)';
    end
end

        %%
        function res= h2(Q,D)
            arguments
                Q
                D=2
            end
            if (Q==0 || Q==1); res=0;
            else
            res = -Q.*log2(Q/(D-1))-(1-Q).*log2(1-Q); 
            end
            
        end

        function res = subs(expr, vals)
            syms ep ed em a t
            vars= [ep, ed, em, a, t];
            res = double(subs(expr,vars,vals));
        end
        function [] = print_bin_vec(kA)
            for i = 1:numel(kA)
                if kA(i)~=0
                    fprintf('%0.2f,%s\n',kA(i),dec2bin(i-1,log2(numel(kA))));
                end
            end
        end

        function [] = figl(n,t)
            arguments
                n
                t='';
            end
                figure(n); clf; hold on; grid on; set(gca,'YScale','log');
                title(t)
        end
        function [res,full] = rps(o,is_plot,options)
            arguments
                o
                is_plot = false
                options.ps =1; %pulse separation
                options.T = 10e-6
            end
            ps =options.ps;
            T = options.T;
            
            tau = 2e-9;
%             photon_per_sec = @(D,m,gain,inc_prob) 1./(T+tau*(D*m+ps)./(inc_prob+gain)); %modified 6.2.23
            photon_per_sec = @(D,m,gain,inc_prob) 1./(T.*(1+inc_prob./gain)+tau*(D*m+ps)./gain); %modified 6.2.23
            if ~isfield(o,'m');o.m=1;end
            if ~isfield(o,'D');o.D=2;end
            full = o.rates_per_photon.*photon_per_sec(o.D,o.m,o.gains,o.inc_probs);
            res = max(full,[],1);
            l=o.lossDBs;
            if is_plot; plot(l,res,'linewidth',1.5);ylim([1e3 2e5]); end
        
        end
        function [res,full] = rps_upper(o,T) %not generic!
            arguments
                o
                T = 10e-6
            end
            mu = o.alphas.^2; mu = repmat(mu',1,13);
            V=1-5e-3;
            Q=0.25e-2;
            d=o.D;
%             T = 10e-6; 
            tau = 2e-9;
            
            x = abs(exp(-mu./2).*sqrt(V)-sqrt(1-exp(-mu))*sqrt(1-V)).^2;
            S = @(rho) -(rho .* log2(rho));
            y = 1-(d-1).*Q;
            holevo=  Q.*(d-1).*log2(d) +...
                    S(y./d * ((d-1).*x+1)) +...
                    (d-1).*S(y./d.*(1-x)) -...
                    S(y);
            rpp = log2(d)+(d-1).*Q.*log2(Q)+y.*(log2(y))-holevo;
            
            
            pps = @(D,gain,inc_prob) 1./(T.*(1+inc_prob./gain)+tau*(D+1)./gain); %modified 6.2.23
            full = rpp.*pps(o.D,o.gains,o.inc_probs);
            res = max(full,[],1);
        end

        function res = rpp(o)   
            
            res= max(o.rates_per_photon,[],1);
            l=o.lossDBs;
            plot(l,res,'linewidth',1.5);
        end
        
        function [res,full] = rpb(o,scale,is_plot)
        arguments
            o
            scale=1;
            is_plot = false;
        end
            full = o.gains.*real(o.rates_per_photon);
            l=o.lossDBs;
            res = max(full,[],1);
            if is_plot; plot(l,res*scale,'linewidth',1.5); end
        end
        function [K,k] = get_tomograpy(rhoA,m,D,is_sym)
           arguments
          rhoA
          m
          D
        is_sym = true
        end
            N= m*D;
            nA = 2^(m*D);
            nB = m*D+2;
            h = helper;

            % get pauli matrices
            Sx = [0,1;1,0];
            Sy = 1j*[0,-1;1,0];
            Sz = [1,0;0,-1];
            S = eye(2);
            S = cat(3,S,Sx);
            S = cat(3,S,Sy);
            S = cat(3,S,Sz);

            % get projection operators, stokes param
            K = zeros(nA*nB,nA*nB,4^N);
%             K = zeros(nA,nA,4^N);
            if is_sym
                k = sym('a',[4^N,1]);
            else
                k = zeros(4^N,1);
            end
            

            %             k = syms
            for i = 1:4^N
                idx= dec2base(i-1,4,N);
                cur_K = 1;
                for j = 1:N
                    cur_S = str2double(idx(j))+1;
                    cur_K = kron(cur_K,S(:,:,cur_S));
                end
%                 cur_K = kron(cur_K,eye(nB));
                k(i) = h.pr(rhoA,cur_K);
%                 K(:,:,i) = cur_K;
                K(:,:,i) = kron(cur_K,eye(nB));
%                 K = cat(3,K,);
            end
            % get stokes params
            
%             for i=1:4^N
%                 k(i)= real(h.pr(rho,K(:,:,i)));
%             end
        
        end
        function [F] = get_F_gain(Md_err,Md_corr,m,D)
            F = 0;
            for i = 1:m*D
                F = F+Md_err(i);
                F = F+Md_corr(i);                
            end
        end
    end    
end
