% Function for matrix vector product


function [Q] = GetMatrixVectorProduct(numpartition,node_gamma,node_interior,R,P)


Q_global = zeros;

for i=1:numpartition
    
n_gamma = unique(node_gamma(any(node_gamma(:,1,i),2),:,i));
n_interior = unique(node_interior(any(node_interior(:,1,i),2),:,i));


Aname = sprintf('Mn_%d.mat',i);
fileA = matfile(Aname);
A = fileA.M_n;

p_s = R(1:size(n_gamma,1),:,i)*P;
u1 = A(n_interior,n_gamma) * p_s;
u2 = A(n_interior,n_interior)\u1;
u3 = A(n_gamma,n_interior)*u2;
u4 = A(n_gamma,n_gamma)*p_s;
Q_local = u4 - u3;
Q_global = Q_global + R(1:size(n_gamma,1),:,i)'*Q_local;

end

Q = Q_global;
end
    
