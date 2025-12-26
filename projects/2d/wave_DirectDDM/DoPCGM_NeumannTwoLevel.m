

function [U_gamma,iter,nres] = DoPCGM_NeumannTwoLevel(npart,node_c,node_r,...
node_corner,node_gamma,node_interior,G_global,U_gamma_o,r_gamma_o,Rc,Rr,Bc,R,D,tol,tolc)



 Z_o = GetPreconditionedTwolevelNeumann(npart,node_c,node_r,...
node_corner,node_gamma,node_interior,r_gamma_o,Rc,Rr,Bc,R,D,tolc);

P_o = Z_o;

rho_o = r_gamma_o' * Z_o;

n = length(G_global);

U_gamma(:,:,1) = U_gamma_o;
P(:,:,1) = P_o;
rho(1) = rho_o;
Z(:,:,1) = Z_o;
r_gamma(:,:,1) = r_gamma_o;
iter = zeros(n,1);
nres = zeros(n,1);

for j = 1:n
    
   sprintf('PCGM Iteration %d',j);
   iter(j) = j;
   
   Q(:,:,j) = GetMatrixVectorProduct(npart,node_gamma,node_interior,R,P(:,:,j)); 
    
   rho_temp(j) = P(:,:,j)' * Q(:,:,j);
   
   alpha(j) = rho(j)/rho_temp(j);
   
   U_gamma(:,:,j+1) = U_gamma(:,:,j)+alpha(j)*P(:,:,j);
   
   r_gamma(:,:,j+1) = r_gamma(:,:,j) - alpha(j)*Q(:,:,j);
   
   Z(:,:,j+1) = GetPreconditionedTwolevelNeumann(npart,node_c,node_r,...
         node_corner,node_gamma,node_interior,r_gamma(:,:,j+1),Rc,Rr,Bc,R,D,tolc);
  
  rho(j+1) = r_gamma(:,:,j+1)' * Z(:,:,j+1);
  
  beta(j) = rho(j+1)/rho(j);
    
  P(:,:,j+1) = Z(:,:,j+1) + beta(j)*P(:,:,j);  
  
  nres(j) = norm(r_gamma(:,:,j+1))/norm(r_gamma(:,:,1));
  
  if(norm(r_gamma(:,:,j+1))/norm(r_gamma(:,:,1)) < tol)
    break;
  end
  
        
end

U_gamma =  U_gamma(:,:,j+1);

end