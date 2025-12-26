% Function for pre-conditioning


function [Z_o] = GetPreconditionedZLumped(numpartition,node_gamma,r_gamma,R,D)


Z_o = zeros;

for i=1:numpartition
    
n_gamma = unique(node_gamma(any(node_gamma(:,1,i),2),:,i));

sc = size(n_gamma,1);

Aname = sprintf('Kt_%d.mat',i);
fileA = matfile(Aname);
A = fileA.K_T;


Z_o = Z_o + R(1:sc,:,i)'*D(1:sc,1:sc,i)'*inv(A(n_gamma,n_gamma))*D(1:sc,1:sc,i)*R(1:sc,:,i)*r_gamma;

end


end
    
