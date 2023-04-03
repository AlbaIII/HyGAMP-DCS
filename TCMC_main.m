clc;
clear;
tic;

L    = 200;
N    = 1000;
M    = 1;
K    = floor(N*0.1);
p    = K/N;
T    = 4;

Pa     = K/N;
P      = zeros(2,2);
P(2,2) = 0.75;
P(1,2) = 1 - P(2,2);
P(2,1) = (1-P(2,2))*Pa/(1-Pa);
P(1,1) = 1 - P(2,1);

SNR = -10;

thre = [linspace(1e-11, 1e-9, 10), linspace(1e-9, 1e-7, 100), linspace(1e-7, 1e-5, 200), linspace(1e-5, 1e-3, 200), linspace(1e-3, 1e-2, 200), linspace(1e-2, 1e-1, 200), linspace(1e-1, 10, 200)];

% Complex Gaussian pilot
A = (randn(L,N)+1j*randn(L,N))/sqrt(2);
A = A./(ones(L,1)*sqrt(sum(abs(A).^2,1)));

num_simu  = 20;
iter      = 100;
nmse_simu = zeros(iter, 5, num_simu);
mse_simu  = zeros(iter, 5, num_simu);
mdp_simu  = zeros(length(thre), 5, num_simu);
fap_simu  = zeros(length(thre), 5, num_simu);
err_simu  = zeros(length(thre), 5, num_simu);
se_simu   = zeros(iter, 5, num_simu);


tau_w = (10^(-SNR/10))/L;

for i = 1 : num_simu
    [nmse_simu(:,:,i), mse_simu(:,:,i), mdp_simu(:,:,i), fap_simu(:,:,i), err_simu(:,:,i), se_simu(:,:,i)] = TCMC_sample(L, N, T, K, A, p, P, SNR, iter, thre);
    disp(i); 
end    
    
% mypar = parpool('local',3);
% parfor i = 1 : num_simu
%     % Complex Gaussian pilot
% %     A = (randn(L,N)+1j*randn(L,N))/sqrt(2);
% %     A = A./(ones(L,1)*sqrt(sum(abs(A).^2,1)));
%     [nmse_simu(:,:,i), mse_simu(:,:,i), mdp_simu(:,:,i), fap_simu(:,:,i), err_simu(:,:,i), se_simu(:,:,i)] = TCMC_sample(L, N, T, K, A, p, P, SNR, iter, thre);
%     disp(i); 
% %     if mod(i,10) == 0
% %         aa = 1;
% %     end
% end
% delete(mypar);

nmse_mean = mean(nmse_simu,3);
mse_mean  = mean(mse_simu,3);
mdp_mean = mean(mdp_simu,3);
fap_mean  = mean(fap_simu,3);
err_mean = mean(err_simu,3);
se_mean  = mean(se_simu,3);

nmse_EM = squeeze(mse_simu(:,5,:));
se_EM   = squeeze(se_simu(:,5,:));
toc;

% save('New_Result_on_TARE_NMSE_vs_p11_1010','nmse_mean','mse_mean','mdp_mean','fap_mean','err_mean','se_mean','L','N','K','p');