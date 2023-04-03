function [ nmse_simu, mse_simu, mdp_simu, fap_simu, err_simu, se ] = TCMC_sample( L, N, T, K, A, p, P, SNR, iter, thre )
%MS_MIMO Summary of this function goes here
%   Detailed explanation goes here

mdp_simu  = zeros(length(thre),5);
fap_simu  = zeros(length(thre),5);
err_simu  = zeros(length(thre),5);
nmse_simu = zeros(iter,5);
mse_simu  = zeros(iter,5);
se        = zeros(iter,5);
% nmse_simu_gamp = zeros(iter,1);
% nmse_simu_ASP  = 0;

T1            = T;
ind_T1        = zeros(N,T1);
Nrand1        = randperm(N);
Ka1           = Nrand1(1:K);
Rand_NT       = rand(N,T1);
LKa1          = length(Ka1);
ind_T1(Ka1,1) = ones(LKa1,1);
for iT = 2 : T
    for iN = 1 : N
        if ind_T1(iN,iT-1) == 0
            if Rand_NT(iN,iT) <= P(2,1)
                ind_T1(iN,iT) = 1;
            end
        else
            if Rand_NT(iN,iT) <= P(2,2)
                ind_T1(iN,iT) = 1;
            end
        end
    end
end
sum_u   = sum(sum(ind_T1));
Ch      = 1/sqrt(2)*(randn(N,T)+1j*randn(N,T));
Ch      = Ch.*ind_T1;
ind_T   = ind_T1;

noise_var = sqrt(10^(-SNR/10))/sqrt(L);
Noise     = noise_var/sqrt(2)*(randn(L,T)+1j*randn(L,T));
Y         = A*Ch;
Y         = Y + Noise;

EM = 0;

% % Simulation for GAMP, HyGAMP-DCS algorithm

% sp = 20;
% eta1 = 0.7;
% eta2 = 0.9;
% [X_est_ASP, nmse_ASP, nmse_ASP_T]                      = PIA_ASP_SMV(Y,A,noise_var^2,Ch,sp,eta1,eta2);

% [X_est_GAMP,nmse_GAMP,mse_GAMP]                        = GAMP_CS(Y, A, P, ones(N,T), K/N, noise_var^2, Ch);
% [X_est_SAMP,nmse_SAMP,mse_SAMP,nmse_SAMPIT,mse_SAMPIT] = S_AMP(Y, A, P, ones(N,T), K/N, noise_var^2, Ch, iter);
[X_est_DCSNE, nmse_DCSNE, mse_DCSNE, se_DCSNE]        = HyGAMP_DCS_SMV( Y, A, P, p, ones(N,T), noise_var^2, Ch, iter, 1, 0 );
[X_est_DCSE, nmse_DCSE, mse_DCSE, se_DCSE]            = HyGAMP_DCS_SMV( Y, A, P, p, ones(N,T), noise_var^2, Ch, iter, 1, 1 );

% % SE to select SNR0
% nmse_DCSE = zeros(iter,6);
% mse_DCSE  = zeros(iter,6);
% se_DCSE   = zeros(iter,6);
% 
% [X_est_DCSE, nmse_DCSE(:,1), mse_DCSE(:,1), se_DCSE(:,1)] = HyGAMP_DCS_SMV( Y, A, P, p, ones(N,T), noise_var^2, Ch, iter, 1, 1, -35 );
% [X_est_DCSE, nmse_DCSE(:,2), mse_DCSE(:,2), se_DCSE(:,2)] = HyGAMP_DCS_SMV( Y, A, P, p, ones(N,T), noise_var^2, Ch, iter, 1, 1, -25 );
% [X_est_DCSE, nmse_DCSE(:,3), mse_DCSE(:,3), se_DCSE(:,3)] = HyGAMP_DCS_SMV( Y, A, P, p, ones(N,T), noise_var^2, Ch, iter, 1, 1, -15 );
% [X_est_DCSE, nmse_DCSE(:,4), mse_DCSE(:,4), se_DCSE(:,4)] = HyGAMP_DCS_SMV( Y, A, P, p, ones(N,T), noise_var^2, Ch, iter, 1, 1, -5 );
% [X_est_DCSE, nmse_DCSE(:,5), mse_DCSE(:,5), se_DCSE(:,5)] = HyGAMP_DCS_SMV( Y, A, P, p, ones(N,T), noise_var^2, Ch, iter, 1, 1, 5 );
% [X_est_DCSE, nmse_DCSE(:,6), mse_DCSE(:,6), se_DCSE(:,6)] = HyGAMP_DCS_SMV( Y, A, P, p, ones(N,T), noise_var^2, Ch, iter, 1, 1, 15 );

% nmse_simu(:,1) = nmse_ASP*ones(iter,1);
% nmse_simu(:,2) = nmse_GAMP;
% nmse_simu(:,3) = nmse_SAMP*ones(iter,1);
nmse_simu(:,4) = nmse_DCSNE;
nmse_simu(:,5) = nmse_DCSE;

% mse_simu(:,1) = zeros(iter,1);
% mse_simu(:,2) = mse_GAMP;
% mse_simu(:,3) = mse_SAMP*ones(iter,1);
mse_simu(:,4) = mse_DCSNE;
mse_simu(:,5) = mse_DCSE;


se(:,4) = se_DCSNE;
se(:,5) = se_DCSE;

for iA = 1 :5
    switch iA
%         case 1
%             X_est_imp = abs(X_est_ASP).^2;
%         case 2
%             X_est_imp = abs(X_est_GAMP).^2;
%         case 3
%             X_est_imp = abs(X_est_SAMP).^2;
        case 4
            X_est_imp = abs(X_est_DCSNE).^2;
        case 5
            X_est_imp = abs(X_est_DCSE).^2;
    end

    for ith = 1 : length(thre)
        acd = zeros(N,T);
%         acd_ub = zeros(N,T);
        temp_Xt = X_est_imp - thre(ith);
        acd_index = find(temp_Xt >= 0);
        acd(acd_index) = ones(length(acd_index),1);

        mdp_simu(ith,iA) = 1 - sum(sum(and(ind_T, acd))) / sum_u;
        fap_simu(ith,iA) = sum(sum(and(1-ind_T, acd))) / (N*T-sum_u);
        err_simu(ith,iA) = ( (sum_u-sum(sum(and(ind_T,acd)))) + sum(sum(and(1-ind_T,acd))) ) / (N*T);

    end
    
end


end

