

function [U_gamma,iter,nres] = DoPCGM_NeumannTwoLevel(npart,U_gamma_o,r_gamma_o,tol,tolc)



 Z_o = GetPreconditionedTwolevelNeumann(npart,r_gamma_o,tolc);

P_o = Z_o;

rho_o = r_gamma_o' * Z_o;

n = length(r_gamma_o);

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
   
   Q(:,:,j) = GetMatrixVectorProduct(npart,P(:,:,j)); 
    
   rho_temp(j) = P(:,:,j)' * Q(:,:,j);
   
   alpha(j) = rho(j)/rho_temp(j);
   
   U_gamma(:,:,j+1) = U_gamma(:,:,j)+alpha(j)*P(:,:,j);
   
   r_gamma(:,:,j+1) = r_gamma(:,:,j) - alpha(j)*Q(:,:,j);
   
   Z(:,:,j+1) = GetPreconditionedTwolevelNeumann(npart,r_gamma(:,:,j+1),tolc);
  
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