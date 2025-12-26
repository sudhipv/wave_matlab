function [M_NN] = GetPreconditionedNeumann(npart,p_f,e_f,t_f,node_gamma,node_interior,r_gamma,R,D)

M_NN = 0;
for i=1:npart
    
p = p_f(any(p_f(:,1,i),2),2:3,i)';
e = e_f(:,:,i)';
t = t_f(any(t_f(:,1:3,i),2),:,i)';

n_gamma = unique(node_gamma(any(node_gamma(:,1,i),2),:,i));
n_interior = unique(node_interior(any(node_interior(:,1,i),2),:,i));

% Retrieving saved sub-domain matrices

Aname = sprintf('A_%d.mat',i);
fileA = matfile(Aname);
A = fileA.A;

S_local = A(n_gamma,n_gamma) - (A(n_gamma,n_interior) * inv(A(n_interior,n_interior)) * A(n_interior,n_gamma));

M_NN = M_NN + R(1:size(n_gamma,1),:,i)'* D(1:size(n_gamma,1),1:size(n_gamma,1),i)'*inv(S_local)*D(1:size(n_gamma,1),1:size(n_gamma,1),i)*R(1:size(n_gamma,1),:,i)*r_gamma;

end

end