% Function for pre-conditioning


function [Z_o] = GetPreconditionedZLumped(numpartition,node_gamma,r_gamma,R)


Z_o = zeros;

for i=1:numpartition
    
n_gamma = unique(node_gamma(any(node_gamma(:,1,i),2),:,i));

Aname = sprintf('Kt_%d.mat',i);
fileA = matfile(Aname);
A = fileA.K_T;


Z_o = Z_o + R(1:size(n_gamma,1),:,i)'*inv(A(n_gamma,n_gamma))*R(1:size(n_gamma,1),:,i)*r_gamma;

end


end
    
