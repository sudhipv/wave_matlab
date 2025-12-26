
% Function to compute Coarse grid residue

function Z_coarse = GetPreconditionedTwolevelNeumann(npart,node_corner_part,node_remaining_part,node_corner,node_gamma,node_interior,r_gamma_o,Rc,Rr,Bc,R,D,tolc)

d_global = 0;
z_global = 0;

for i=1:npart

n_c = unique(node_corner_part(any(node_corner_part(:,1,i),2),:,i));
n_r = unique(node_remaining_part(any(node_remaining_part(:,1,i),2),:,i));
n_i = unique(node_interior(any(node_interior(:,1,i),2),:,i));
n_gamma = unique(node_gamma(any(node_gamma(:,1,i),2),:,i));

sg = size(n_gamma,1);
sr = size(n_r,1);
sc = size(n_c,1);
% Retrieving saved sub-domain matrices

Aname = sprintf('Mn_%d.mat',i);
fileA = matfile(Aname);
A = fileA.M_n;

s_rr = A(n_r,n_r) - (A(n_r,n_i) * inv(A(n_i,n_i)) * A(n_i,n_r));

s_cr = A(n_c,n_r) - (A(n_c,n_i) * inv(A(n_i,n_i)) * A(n_i,n_r));

s_rc = A(n_r,n_c) - (A(n_r,n_i) * inv(A(n_i,n_i)) * A(n_i,n_c));

r_gamma_local = D(1:sg,1:sg,i) * R(1:sg,:,i)* r_gamma_o;

f_r(:,:,i) = Rr(:,1:sg,i)*r_gamma_local;

f_c(:,:,i) = Rc(:,1:sg,i)*r_gamma_local;

v1 = s_rr\f_r(1:sr,:,i);

d_local = f_c(1:sc,:,i) - (s_cr * v1);

d_global = d_global + Bc(1:sc,:,i)'*d_local;

end

%%%%%%%%%%%%%%%%%%%

% Initial assumption
U_o = zeros(length(node_corner),1);
r_gamma_o = d_global;
z_c = DoPCGMLumpedTwoLevel(npart,U_o,r_gamma_o,node_corner_part,node_remaining_part,node_interior,Bc,tolc);

%%%%%%%%%%%%%%%%%%%

for i=1:npart
     
n_c = unique(node_corner_part(any(node_corner_part(:,1,i),2),:,i));
n_r = unique(node_remaining_part(any(node_remaining_part(:,1,i),2),:,i));
n_i = unique(node_interior(any(node_interior(:,1,i),2),:,i));
n_gamma = unique(node_gamma(any(node_gamma(:,1,i),2),:,i));

sg = size(n_gamma,1);
sr = size(n_r,1);
sc = size(n_c,1);

% Retrieving saved sub-domain matrices

Aname = sprintf('Mn_%d.mat',i);
fileA = matfile(Aname);
A = fileA.M_n;

s_rr = A(n_r,n_r) - (A(n_r,n_i) * inv(A(n_i,n_i)) * A(n_i,n_r));

s_rc = A(n_r,n_c) - (A(n_r,n_i) * inv(A(n_i,n_i)) * A(n_i,n_c));


z_c_local = Bc(1:sc,:,i) * z_c;

v2 = f_r(1:sr,:,i) - s_rc * z_c_local;

z_r_local = s_rr\v2;

z_local = Rr(1:sr,1:sg,i)'*z_r_local + Rc(1:sc,1:sg,i)'*z_c_local;

z_global = z_global +  R(1:sg,:,i)' * D(1:sg,1:sg,i) * z_local;

end

Z_coarse = z_global;
end