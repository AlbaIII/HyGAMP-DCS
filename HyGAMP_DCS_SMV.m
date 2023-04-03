function [ X_est, nmse, mse, se_x_avg] = HyGAMP_DCS_SMV( Y, A, P, pa, g, tao_w, X_0, iter, choice, EM, SNR0 )
%HYGAMP_DCS_SMV Summary of this function goes here
%   Detailed explanation goes here

[L, T] = size(Y);
[~, N]   = size(A);
damp     = 1;
beta     = L/N;

faixsr       = @(u,s,rho) u.*exp(-u)./(1+(1./rho-1).*(1+s).*exp(-s.*u));
V2           = @(tau,g,rho,fai) rho.*g.*tau./(g+tau)+rho.*(g.^2)./(g+tau).*(1-fai); 

if EM == 1
    tao_w  = 10^(-(SNR0)/10)/L;
    pa     = 0.1;
    P      = zeros(2,2);
    P(2,2) = pa;
    P(1,2) = 1 - P(2,2);
    P(2,1) = (1-P(2,2))*pa/(1-pa);
    P(1,1) = 1 - P(2,1);
end

epi     = 0.01;
mu      = zeros(N,T);
% g       = g*ones(1,T);
X_est   = zeros(N,T);
nmse    = zeros(iter,1);
mse     = zeros(iter,1);
tao_r1 = zeros(iter,1);
tao_r2 = zeros(iter,1);
tao_x1 = zeros(iter,1);

% case of vector variance
z     = zeros(L,T);
tao_z = zeros(L,T);
x_hat = zeros(N,T);
tao_p = zeros(L,T);
p_hat = zeros(L,T);
s_hat = zeros(L,T);
tao_s = zeros(L,T);
tao_r = zeros(N,T);
r_hat = zeros(N,T);
x_bar = zeros(N,T);
tao_x = pa*g;
tao_x_avg = zeros(iter,1);
se_x = tao_x;
se_r = zeros(N,T);
se_x_avg = zeros(iter,1);


LLR_tnm    = log(pa/(1-pa))*ones(N,T);
rau        = 1 - 1./(1+exp(LLR_tnm));
pai_nt_r   = pa*ones(N,T);
pai_tn_l   = pa*ones(N,T);


for it = 1 : iter
    
    if it ~= 1
        damp = 1;
    else
        damp = 1;
    end
    
%     GAMP procedure
%     factor node update
    z     = A*x_hat;
    tao_p = damp*(abs(A).^2 * tao_x) + (1-damp)*tao_p;
    p_hat = z - tao_p.*s_hat;
    
    z_0   = (tao_p.*Y+tao_w.*p_hat) ./ (tao_p+tao_w);
    tao_z = (tao_p.*tao_w) ./ (tao_p+tao_w);
    s_hat = damp*((z_0-p_hat)./tao_p) + (1-damp)*s_hat;
    tao_s = damp*((1-tao_z./tao_p)./tao_p) + (1-damp)*tao_s;
    
%     variable node update
    x_bar = damp*x_hat + (1-damp)*x_bar;
    tao_r = 1 ./ ((abs(A').^2)*tao_s); 
    r_hat = x_bar + tao_r.*(A'*s_hat);
    
%     Assuming that the channel coefficeints satisfy 
%     Bernoulli-Gaussian distributions.
    veta  = (g.*tao_r) ./ (g+tao_r);
    gamma = veta .* (r_hat./tao_r);
    delta = 1./(tao_r) - 1./(g+tao_r);
    pai   = 1 ./ (1+((1-rau)./rau).*((g+tao_r)./(tao_r)).*exp(-(abs(r_hat).^2).*delta));
    x_hat = pai.*gamma;
    tao_x = pai.*((1-pai).*(abs(gamma).^2)+veta);
    
    nmse(it) = 10*log10((norm(X_0-x_hat, 'fro').^2) / (norm(X_0,'fro').^2));
  
    se_x_t = ones(N,1)*mean(se_x,1);
    se_r = tao_w + se_x_t/beta;
    for iT = 1 : T
        for iN = 1 : N
            faix = @(u) faixsr(u,g(iN,iT)./se_r(iN,iT),rau(iN,iT));
            fai1 = quad(faix,0,30);
            se_x(iN,iT) = V2(se_r(iN,iT),g(iN,iT),rau(iN,iT),fai1);
        end
    end
    se_x_avg(it) = 10*log10(mean(mean(se_x,2),1)) + 10;
    
    if choice == 1
        
        LLR_mnt_l = log(tao_r)-log(tao_r+g)-abs(r_hat-mu).^2./(tao_r+g)+abs(r_hat).^2./tao_r;
        for i1 = 1 : N
            for i2 = 1 : T
                if LLR_mnt_l(i1,i2) > 200
                 	LLR_mnt_l(i1,i2) = 200;
                end
            end
        end    
        pai_mnt_l = 1 - 1./(1+exp(LLR_mnt_l));
        
        for iT = 1 : T
            if iT == 1
                pai_nt_r(:,iT) = pa*ones(N,1);
            else
                pai_nt_r(:,iT) = (P(2,1).*(1-pai_mnt_l(:,iT-1)).*(1-pai_nt_r(:,iT-1)) + P(2,2).*pai_mnt_l(:,iT-1).*pai_nt_r(:,iT-1)) ...
                                 ./ ((1-pai_mnt_l(:,iT-1)).*(1-pai_nt_r(:,iT-1))+pai_mnt_l(:,iT-1).*pai_nt_r(:,iT-1));
            end
        end

        for iT = T : -1 : 1
            if iT == T
                p_temp11 = (P(1,1)+P(1,2)).*(1-pai_mnt_l(:,iT)) + (P(2,2)+P(2,1)).*pai_mnt_l(:,iT);
                pai_tn_l(:,iT) = (P(1,2)*(1-pai_mnt_l(:,iT)) + P(2,2)*pai_mnt_l(:,iT)) ./ p_temp11;
            else
                p_temp21 = (P(1,1)+P(1,2)).*(1-pai_tn_l(:,iT+1)).*(1-pai_mnt_l(:,iT)) + (P(2,2)+P(2,1)).*pai_tn_l(:,iT+1).*pai_mnt_l(:,iT);
                pai_tn_l(:,iT) = (P(1,2).*(1-pai_tn_l(:,iT+1)).*(1-pai_mnt_l(:,iT)) + P(2,2).*pai_tn_l(:,iT+1).*pai_mnt_l(:,iT)) ./ ...
                                 p_temp21;
            end
        end
        
        for iT = 1 : T
            if iT ~= T
                LLR_tnm(:,iT) = log(pai_nt_r(:,iT)./(1-pai_nt_r(:,iT)))+log(pai_tn_l(:,iT+1)./(1-pai_tn_l(:,iT+1)));
            else
                LLR_tnm(:,iT) = log(pai_nt_r(:,iT)./(1-pai_nt_r(:,iT)));
            end
        end
        rau = 1 - 1./(1+exp(LLR_tnm));
    else
%         aa = 1;
%         varpiM    = pai;
%         rau       = mean(varpiM,2)*ones(1,T);
%         rau = rau;
    end
    
 	if EM == 1
        tao_w = mean(mean(abs(Y-z_0).^2+tao_z));
        g = sum(sum(pai.*(abs(gamma).^2+veta)))/sum(sum(pai))*ones(N,T);
        mu = sum(sum(gamma.*pai))/sum(sum(pai)) *ones(N,T);
        pa = mean(pai_mnt_l(:,1).*pai_nt_r(:,1).*pai_tn_l(:,2) ./ (pai_mnt_l(:,1).*pai_nt_r(:,1).*pai_tn_l(:,2)+(1-pai_mnt_l(:,1)).*(1-pai_nt_r(:,1)).*(1-pai_tn_l(:,2))));
        p_temp1 = pai_nt_r(:,1:(T-1)).*pai_mnt_l(:,1:(T-1)) ./ (pai_nt_r(:,1:(T-1)).*pai_mnt_l(:,1:(T-1))+(1-pai_nt_r(:,1:(T-1))).*(1-pai_mnt_l(:,1:(T-1))));
        p_temp2 = zeros(N,T-1);
        p_temp2(:,1:(T-2)) = pai_tn_l(:,3:T).*pai_mnt_l(:,2:(T-1)) ./ (pai_tn_l(:,3:T).*pai_mnt_l(:,2:(T-1)) + (1-pai_tn_l(:,3:T)).*(1-pai_mnt_l(:,2:(T-1))));
        p_temp2(:,T-1) = pai_mnt_l(:,T);
        PP = P(2,2).*p_temp1.*p_temp2+P(1,2).*p_temp1.*(1-p_temp2)+P(2,1).*(1-p_temp1).*p_temp2+P(1,1).*(1-p_temp1).*(1-p_temp2);
        p_s_temp = pai_mnt_l(:,1:(T-1)).*pai_nt_r(:,1:(T-1)).*pai_tn_l(:,2:T) ./ ((1-pai_mnt_l(:,1:(T-1))).*(1-pai_nt_r(:,1:(T-1))).*(1-pai_tn_l(:,2:T)) + pai_mnt_l(:,1:(T-1)).*pai_nt_r(:,1:(T-1)).*pai_tn_l(:,2:T));
        Es = p_s_temp;
        Ess = P(2,2).*p_temp1.*p_temp2./PP;
        P(1,2) = (sum(sum(Es))-sum(sum(Ess)))./sum(sum(Es));
        P(2,2) = 1-P(1,2);
        P(2,1) = pa*P(1,2)/(1-pa);
        P(1,1) = 1 - P(2,1);
    end 

end

X_est = x_hat;

end

