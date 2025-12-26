
%%%%% Function to Update %%%%%%%%

function [u_n,p_n,a_n] = Update(u_0,p_0,a_0,deltaT,beta,gamma,a_new)

a_beta = (1-2*beta)*a_0 + 2*beta*a_new;

a_gamma = (1-gamma)*a_0+gamma*a_new;

p_n = p_0 + deltaT* a_gamma;

u_n = u_0 + deltaT*p_n + (deltaT^2/2)*a_beta;

a_n = a_new;


end
