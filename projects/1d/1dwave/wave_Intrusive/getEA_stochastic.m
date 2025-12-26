

function EA_ss = getEA_stochastic(i,ipce,EA_0,l,x,mu_g,sigma_g,norm,mindex,dim)


%%% Log Normal Random Variable- General Way

  L_0 = exp(mu_g + 0.5 * sigma_g^2);

    
  m_index = mindex.m_index_out;
    
    numerator = 1;
    for ii = 1:dim
        numerator = numerator * sigma_g^(m_index(ipce,ii));
    end
        EA_ss = L_0*(numerator/sqrt(norm.norm_squared(ipce)));


%      EA_plot(i) = EA_ss;
     
%%% Log Normal Random Variable - Analytical     
     
% mu_l = exp(mu_g + sigma_g^2/2);
% kappa(1) = mu_l;
% kappa(2) = mu_l*sigma_g;
% kappa(3) = mu_l*sigma_g^2/sqrt(2);   
% kappa(4) = mu_l*sigma_g^3/sqrt(6);
% kappa(5) = mu_l*sigma_g^4/sqrt(24);
% kappa(6) = mu_l*sigma_g^5/sqrt(120);
% kappa(7) = mu_l*sigma_g^6/sqrt(720);
% kappa(8) = mu_l*sigma_g^7/sqrt(5040);
% kappa(9) = mu_l*sigma_g^8/sqrt(40320);
% kappa(10) = mu_l*sigma_g^9/sqrt(362880);
% kappa(11) = mu_l*sigma_g^10/sqrt(3628800);
% 
%      EA_ss = kappa(ipce);
     
%%%%%%% Random process

% % % % switch ipce
% % % %     
% % % %     case 1 
% % % %         EA_ss = EA_0;
% % % %     case 2
% % % %         
% % % %         EA_ss = EA_0*abs(cos(4*pi*x_s));
% % % % 
% % % %    case 3
% % % %         EA_ss = EA_0*abs(cos(8*pi*x_s));

        
        
        
end