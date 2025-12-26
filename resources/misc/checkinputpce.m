

%%%% Checking whether the input approximation by PCE is enough for log
%%%% normal random variable
clearvars

mu_g = 0;
sigma_g = 0.3;

% MCS
n_MC = 2000000;
xi_MC = randn(1,n_MC);
LN_MCS = exp(mu_g+sigma_g.*xi_MC);
[f_MC,x_MC] = ksdensity(LN_MCS);


% Log PDF matlab
x_matlab = 0:0.001:10;
y = lognpdf(x_matlab,mu_g,sigma_g);


% PCE
n = 50000;
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
        (sigma_g^4/24).*psi(:,5)+(sigma_g^5/120).*psi(:,6)+ (sigma_g^6/720).*psi(:,7)+ (sigma_g^7/5040).*psi(:,8)); ...
%         + (sigma_g^8/40320).*psi(:,9) + (sigma_g^9/362880).*psi(:,10) + (sigma_g^10/3628800).*psi(:,11));
%    
    
    
[f_pce,x_pce] = ksdensity(LN_RV);




% openfig('MCandAnalytical.fig')
% hold on
% plot(x_pce,f_pce)
% hold on

figure(5)
plot(x_matlab,y,'linewidth',2)
hold on
plot(x_MC,f_MC,'linewidth',2)
hold on
plot(x_pce,f_pce,'linewidth',2)
grid on;
axis([-1 4 min(y) 1.5]);
% xlabel('X axis','fontSize',14);
% ylabel('PDF','fontSize',14);   
legend                   
  h=gca; 
  get(h,'FontSize');
  set(h,'FontSize',14);



% % % % % mupdf_MC = mean(LN_MCS)
% % % % % sigmapdf_MC = sqrt(var(LN_MCS))
% % % % % skew_MC = skewness(LN_MCS)
% % % % % kur_MC = kurtosis(LN_MCS)
% % % % % 
% % % % % 
% % % % % [m_mat,v_mat] = lognstat(mu_g,sigma_g)
% % % % % 
% % % % % v_mat = sqrt(v_mat)


% mupdf_pce = mean(LN_RV)
% sigmapdf_pce = sqrt(var(LN_RV))
% skew_pce = skewness(LN_RV)
% kur_pce = kurtosis(LN_RV)

% % % % figure(2)
% % % % plot(x_pce,f_pce)
% % % % hold on
% % % % plot(x_matlab,y,'k')

