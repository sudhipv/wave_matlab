
%%%% Experiment for convergence of PC 




%%%% Checking whether the input approximation by PCE is enough for log
%%%% normal random variable
clearvars

mu_g = 0;
sigma_g = 0.1;

% number of samples
n = 10000;
batch = 100;

batch_LNRV = zeros(n,batch);

for bb = 1: batch


xi_pce = randn(1,n)';


psi(:,1) = ones(1,n);
psi(:,2) = xi_pce;
psi(:,3) = xi_pce.^2-1;
psi(:,4) = xi_pce.^3-3.*xi_pce;
psi(:,5) = xi_pce.^4-6*xi_pce.^2+3;
psi(:,6) = xi_pce.^5-10*xi_pce.^3+15*xi_pce;
psi(:,7) = xi_pce.^6-15*xi_pce.^4+45*xi_pce.^2-15;
psi(:,8) = xi_pce.^7-21*xi_pce.^5+105*xi_pce.^3-105*xi_pce;
psi(:,9) = xi_pce.^8-28*xi_pce.^6+210*xi_pce.^4-420*xi_pce.^2+105;
psi(:,10) = xi_pce.^9-36*xi_pce.^7+378*xi_pce.^5-1260*xi_pce.^3+945*xi_pce;
psi(:,11) = xi_pce.^10-45*xi_pce.^8+630*xi_pce.^6-3150*xi_pce.^4+4725*xi_pce.^2-945;


mu_l = exp(mu_g + sigma_g^2/2);

LN_RV = mu_l .* (1.*psi(:,1) + (sigma_g^1/1).*psi(:,2)+(sigma_g^2/2).*psi(:,3)+(sigma_g^3/6).*psi(:,4)+...
        (sigma_g^4/24).*psi(:,5));
%     +(sigma_g^5/120).*psi(:,6));
%     + (sigma_g^6/720).*psi(:,7)+ (sigma_g^7/5040).*psi(:,8) ...
%         + (sigma_g^8/40320).*psi(:,9) + (sigma_g^9/362880).*psi(:,10) + (sigma_g^10/3628800).*psi(:,11));
%    
    

batch_LNRV(:,bb) = LN_RV;

end

batch_sd = sqrt(var(batch_LNRV));











