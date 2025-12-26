function [M_NN] = GetPreconditionedNeumann(npart,node_gamma,node_interior,r_gamma,R,D,NBcount)

M_NN = 0;
for i=1:npart
    
n_gamma = unique(node_gamma(any(node_gamma(:,1,i),2),:,i));
n_interior = unique(node_interior(any(node_interior(:,1,i),2),:,i));

% Retrieving saved sub-domain matrices

Aname = sprintf('Kt_%d.mat',i);
fileA = matfile(Aname);
A = fileA.K_T;

S_local = A(n_gamma,n_gamma) - (A(n_gamma,n_interior) * inv(A(n_interior,n_interior)) * A(n_interior,n_gamma));


if(NBcount ==2)
    strcond = 'Condition number of Local Schur complement for subdomain %d is %d \n';
    fprintf(strcond,i,cond(S_local))
end

% % % for k=1:size(S_local,1)
% % %    S_local(k,k) = S_local(k,k) + 0.001; 
% % % end

M_NN = M_NN + R(1:size(n_gamma,1),:,i)'* D(1:size(n_gamma,1),1:size(n_gamma,1),i)'*inv(S_local)*D(1:size(n_gamma,1),1:size(n_gamma,1),i)*R(1:size(n_gamma,1),:,i)*r_gamma;

end

end