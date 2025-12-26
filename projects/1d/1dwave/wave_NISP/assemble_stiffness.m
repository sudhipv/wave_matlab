

function K = assemble_stiffness(ne,nN,K,xi,mu_g,sigma_g,x,l,ipce,psi,ss,n_inpce,expnsn)

EA_ss = 0;

switch expnsn
    
    case 1
        %%% Direct Sampling
EA_ss = exp(mu_g+ sigma_g*xi);


 K = K * EA_ss;
        
    case 2
%%% PCE expansion


mu_l = exp(mu_g + sigma_g^2/2);
kappa(1) = mu_l;
kappa(2) = mu_l*sigma_g;
kappa(3) = mu_l*sigma_g^2/sqrt(2);   
kappa(4) = mu_l*sigma_g^3/sqrt(6);
kappa(5) = mu_l*sigma_g^4/sqrt(24);
kappa(6) = mu_l*sigma_g^5/sqrt(120);
kappa(7) = mu_l*sigma_g^6/sqrt(720);
kappa(8) = mu_l*sigma_g^7/sqrt(5040);
kappa(9) = mu_l*sigma_g^8/sqrt(40320);
kappa(10) = mu_l*sigma_g^9/sqrt(362880);
kappa(11) = mu_l*sigma_g^10/sqrt(3628800);

%%% Non Normalized
% % mu_l = exp(mu_g + sigma_g^2/2);
% % kappa(1) = mu_l;
% % kappa(2) = mu_l*sigma_g;
% % kappa(3) = mu_l*sigma_g^2/2;   
% % kappa(4) = mu_l*sigma_g^3/6;
% % kappa(5) = mu_l*sigma_g^4/24;
% % kappa(6) = mu_l*sigma_g^5/120;
% % kappa(7) = mu_l*sigma_g^6/720;
% % kappa(8) = mu_l*sigma_g^7/5040;
% % kappa(9) = mu_l*sigma_g^8/40320;
% % kappa(10) = mu_l*sigma_g^9/362880;
% % kappa(11) = mu_l*sigma_g^10/3628800;



for pp = 1:n_inpce

EA_ss = EA_ss + kappa(pp)*psi(ss,pp);

end



%%%% Random Process %%%%



% % % for i = 1:ne
% % %     
% % %     %%% Center of each element is taken as ea varying point
% % %    x_s = x(i) + l*0.5;
% % %    alpha1 = 2;
% % %    alpha2 = .5;
% % %   %%% Multiplied by each random amplitude to form stiffness matrix
% % %    EA_ss = 1 + xi*cos(x(i)/l) + (xi^2 - 1)*cos(x(i)/l);
% % %    
% % %    easss(i) = EA_ss*5;
% % % 
% % %     K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + k_e*EA_ss;
% % % 
% % % 
% % % end

% %  K(1,1) = K(1,1) + 1* 10^(7);

end

K = K * EA_ss;

end



