
% Function for PCGM solver using Lumped Preconditioner  

function [U,iter,nres] = DoPCGMNeumann(npart,U_o,d_global,node_interface_part,node_interior_part,R,D,tol,count,NBcount)

r_gamma_o = d_global;

Z_o = GetPreconditionedNeumann(npart,node_interface_part,node_interior_part,r_gamma_o,R,D,NBcount);

P_o = Z_o;

rho_o = r_gamma_o' * Z_o;

n = length(d_global);

U(:,:,1) = U_o;
P(:,:,1) = P_o;
rho(1) = rho_o;
Z(:,:,1) = Z_o;
r_gamma(:,:,1) = r_gamma_o;

iter = zeros(n,1);
nres = zeros(n,1);

for j = 1:n
    
%    fprintf('PCGM Iteration %d\n',j);
   iter(j) = j;
   
   Q(:,:,j) = GetMatrixVectorProduct(npart,node_interface_part,node_interior_part,R,P(:,:,j)); 
    
   rho_temp(j) = P(:,:,j)' * Q(:,:,j);
   
   alpha(j) = rho(j)/rho_temp(j);
   
   U(:,:,j+1) = U(:,:,j)+alpha(j)*P(:,:,j);
   
   r_gamma(:,:,j+1) = r_gamma(:,:,j) - alpha(j)*Q(:,:,j);
   
  Z(:,:,j+1) = GetPreconditionedNeumann(npart,node_interface_part,node_interior_part,r_gamma(:,:,j+1),R,D,NBcount);
  
  rho(j+1) = r_gamma(:,:,j+1)' * Z(:,:,j+1);
  
  beta(j) = rho(j+1)/rho(j);
    
  P(:,:,j+1) = Z(:,:,j+1) + beta(j)*P(:,:,j);  
  
  nres(j) = norm(r_gamma(:,:,j+1))/norm(r_gamma(:,:,1));
  
  if(norm(r_gamma(:,:,j+1))/norm(r_gamma(:,:,1)) < tol)
    break;
  end
    
end

U =  U(:,:,j+1);
end