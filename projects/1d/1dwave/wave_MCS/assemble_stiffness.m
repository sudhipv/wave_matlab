

function K = assemble_stiffness(ne,nN,K,xi,mu_g,sigma_g,x,l)



EA_ss = exp(mu_g + sigma_g*xi);

K = K * EA_ss;


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



