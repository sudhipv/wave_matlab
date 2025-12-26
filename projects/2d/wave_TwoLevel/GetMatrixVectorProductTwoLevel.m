% Function for matrix vector product


function [Q] = GetMatrixVectorProductTwoLevel(numpartition,node_corner_part,node_remaining_part,node_interior,Bc,P)

Q_global = zeros;

for i=1:numpartition
    
n_c = unique(node_corner_part(any(node_corner_part(:,1,i),2),:,i));
n_r = unique(node_remaining_part(any(node_remaining_part(:,1,i),2),:,i));
n_i = unique(node_interior(any(node_interior(:,1,i),2),:,i));

sr = size(n_r,1);
sc = size(n_c,1);

Aname = sprintf('Kt_%d.mat',i);
fileA = matfile(Aname);
A = fileA.K_T;

s_rr = A(n_r,n_r) - (A(n_r,n_i) * inv(A(n_i,n_i)) * A(n_i,n_r));

s_rc = A(n_r,n_c) - (A(n_r,n_i) * inv(A(n_i,n_i)) * A(n_i,n_c));

s_cr = A(n_c,n_r) - (A(n_c,n_i) * inv(A(n_i,n_i)) * A(n_i,n_r));

s_cc = A(n_c,n_c) - (A(n_c,n_i) * inv(A(n_i,n_i)) * A(n_i,n_c));

p_s = Bc(1:sc,:,i)*P;
u1 = s_rc * p_s;
u2 = s_rr\u1;
u3 = s_cr*u2;
u4 = s_cc*p_s;
Q_local = u4 - u3;
Q_global = Q_global + Bc(1:sc,:,i)'*Q_local;

end

Q = Q_global;
end
    
