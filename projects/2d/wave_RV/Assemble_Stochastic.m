
%%%%%%%%% Function to Assemble Stochastic Matrix, Ajk %%%%%%%%

function Assemble_Stochastic(s_spectral,Cijk,npart,node_gamma,node_interior,num_interface,R,beta,deltaT)


%%% Assembling the Domain Decomposed Stochastic Sub matrices %%%%
%%%%% Assembling the interface sub matrix will allow us to solve for the
%%%%% interface problem and then parallely solving for all the remaining
%%%%% sub domain level matrices %%%%%%%%%%%%

n_outpce = s_spectral.n_outpce;
n_inpce = s_spectral.n_inpce;

% % % % Expectation of psi_k <b psi_k>/<psi_k^2>
% % % exp_pce = zeros(1,n_outpce);
% % % exp_pce(1,1) = 1;

mkdir AbSto

for np = 1:npart

    
  n_gamma = unique(node_gamma(any(node_gamma(:,1,np),2),:,np));
  n_interior = unique(node_interior(any(node_interior(:,1,np),2),:,np));   
  
  
  cd ./AbDet
  Aname = sprintf('A_%d.mat',np);
  Mname = sprintf('M_%d.mat',np);
  Amat = load(Aname);
  Mmat = load(Mname);
  A = Amat.A;
  M = Mmat.M;
  cd ../
   
  
  nAgg = size(A(n_gamma,n_gamma,1),1);
  
  nAii = size(A(n_interior,n_interior,1),1);
  
  
  %%% Stiffness/Diffusion Coefficient
size_gamma = nAgg*n_outpce;
size_intr = nAii*n_outpce;

Ajk_gamma = zeros(size_gamma,size_gamma);
Ajk_intr = zeros(size_intr,size_intr);
Ajk_igamma = zeros(size_intr,size_gamma);
Ajk_gintr = zeros(size_gamma,size_intr);


Mjk_gamma = zeros(size_gamma,size_gamma);
Mjk_intr = zeros(size_intr,size_intr);
Mjk_igamma = zeros(size_intr,size_gamma);

% Force Vector
% % % bk_gamma = zeros(size_gamma,1);
% % % bk_intr = zeros(size_intr,1);
  

% Restriction Matrix - stochastic
R_row = nAgg*n_outpce;
R_col = num_interface*n_outpce;
R_sto = zeros(R_row,R_col);

    
% Assembly block struture
for k = 1: n_outpce
    
        
      
    for j = 1:n_outpce
        
        for i= 1: n_inpce
       
    
     %% A_gamma_gamma
    Ajk_gamma((k-1)*nAgg+1:(k-1)*nAgg+nAgg, (j-1)*nAgg+1:(j-1)*nAgg+nAgg)...
        = Ajk_gamma((k-1)*nAgg+1:(k-1)*nAgg+nAgg, (j-1)*nAgg+1:(j-1)*nAgg+nAgg) +...
                         Cijk.moment_PsiNorm(i,j,k) * A(n_gamma,n_gamma,i);
                                    
     %% A_interior_interior                
   Ajk_intr((k-1)*nAii+1:(k-1)*nAii+nAii, (j-1)*nAii+1:(j-1)*nAii+nAii)...
        = Ajk_intr((k-1)*nAii+1:(k-1)*nAii+nAii, (j-1)*nAii+1:(j-1)*nAii+nAii) +...
                         Cijk.moment_PsiNorm(i,j,k) * A(n_interior,n_interior,i);
                     
      %% A_interior_gamma                
   Ajk_igamma((k-1)*nAii+1:(k-1)*nAii+nAii, (j-1)*nAgg+1:(j-1)*nAgg+nAgg)...
        = Ajk_igamma((k-1)*nAii+1:(k-1)*nAii+nAii, (j-1)*nAgg+1:(j-1)*nAgg+nAgg) +...
                         Cijk.moment_PsiNorm(i,j,k) * A(n_interior,n_gamma,i);               
    
      %% A_gamma_interior                
   Ajk_gintr((k-1)*nAgg+1:(k-1)*nAgg+nAgg, (j-1)*nAii+1:(j-1)*nAii+nAii)...
        = Ajk_gintr((k-1)*nAgg+1:(k-1)*nAgg+nAgg, (j-1)*nAii+1:(j-1)*nAii+nAii) +...
                         Cijk.moment_PsiNorm(i,j,k) * A(n_gamma,n_interior,i);   
    
        end
        
    end
    
    %%% Force vector is same for all PCE terms
% % %    bk_gamma((k-1)*nAgg+1:(k-1)*nAgg+nAgg) = b(n_gamma,1)* exp_pce(k);
% % %    bk_intr((k-1)*nAii+1:(k-1)*nAii+nAii) = b(n_interior,1)* exp_pce(k);
% % %    
   
   %%% Restriction Matrix assembly
   
   R_sto((k-1)*nAgg+1:k*nAgg,(k-1)*num_interface+1:k*num_interface) = R(1:nAgg,:,np);
   
   
   Mjk_gamma((k-1)*nAgg+1:k*nAgg,(k-1)*nAgg+1:k*nAgg) = ...
                                  M(n_gamma,n_gamma);
%                               *sqrt(s_spectral.norm.norm_squared(k));
    
   Mjk_intr((k-1)*nAii+1:k*nAii,(k-1)*nAii+1:k*nAii) = ...
                                  M(n_interior,n_interior);
%                                                        
                              
   Mjk_igamma((k-1)*nAii+1:k*nAii,(k-1)*nAgg+1:k*nAgg) = ...
                                  M(n_interior,n_gamma);
                                                        
                                    
                              
end


Mjk_gintr = Mjk_igamma';

%%%%%% Finding out K_T


K_T_gamma = Mjk_gamma/(beta*deltaT^2) + Ajk_gamma;

K_T_intr = Mjk_intr/(beta*deltaT^2) + Ajk_intr;

K_T_igamma = Mjk_igamma/(beta*deltaT^2) + Ajk_igamma;

K_T_gintr = Mjk_gintr/(beta*deltaT^2) + Ajk_gintr;

cd ./AbSto
foldername = sprintf('part_%d',np);
mkdir(foldername);
cd (foldername)

save('K_T_gg','K_T_gamma')
save('K_T_ii','K_T_intr')
save('K_T_ig','K_T_igamma')
save('K_T_gi','K_T_gintr')

save('Mjk_gg','Mjk_gamma')
save('Mjk_ii','Mjk_intr')
save('Mjk_ig','Mjk_igamma')
save('Mjk_gi','Mjk_gintr')

% save('Bk_g','bk_gamma')
% save('Bk_i','bk_intr')
save('R_pce','R_sto')

cd ../../

end


end