
% Function for PCGM solver using Lumped Preconditioner  

function [U] = DoPCGMLumpedTwoLevel(npart,U_o,d_global,node_corner_part,node_remaining_part,node_interior,Bc,tolc)

r_gamma_o = d_global;

z_o = GetPreconditionedZLumped(npart,node_corner_part,r_gamma_o,Bc);

P_o = z_o;

rho_o = r_gamma_o' * z_o;

n = length(d_global);

U(:,:,1) = U_o;
P(:,:,1) = P_o;
rho(1) = rho_o;
z(:,:,1) = z_o;
r_gamma(:,:,1) = r_gamma_o;

for j = 1:n
    
   Q(:,:,j) = GetMatrixVectorProductTwoLevel(npart,node_corner_part,node_remaining_part,node_interior,Bc,P(:,:,j)); 
    
   rho_temp(j) = P(:,:,j)' * Q(:,:,j);
   
   alpha(j) = rho(j)/rho_temp(j);
   
   U(:,:,j+1) = U(:,:,j)+alpha(j)*P(:,:,j);
   
   r_gamma(:,:,j+1) = r_gamma(:,:,j) - alpha(j)*Q(:,:,j);
   
  z(:,:,j+1) = GetPreconditionedZLumped(npart,node_corner_part,r_gamma(:,:,j+1),Bc);
  
  rho(j+1) = r_gamma(:,:,j+1)' * z(:,:,j+1);
  
  beta(j) = rho(j+1)/rho(j);
    
  P(:,:,j+1) = z(:,:,j+1) + beta(j)*P(:,:,j);  
  
  if(norm(r_gamma(:,:,j+1))/norm(r_gamma(:,:,1)) < tolc)
    break;
  end
    
end

U =  U(:,:,j+1);
end