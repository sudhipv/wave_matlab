
%%%%%%%%% Function to Assemble Stochastic Matrix, Ajk %%%%%%%%

function Assemble_Stochastic(s_spectral,Cijk,npart,node_gamma,node_interior,node_c,node_r,num_interface,...
    node_corner,R,Rc,Rr,Bc,D,beta,deltaT)


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
  n_cor = unique(node_c(any(node_c(:,1,np),2),:,np));
  n_rem = unique(node_r(any(node_r(:,1,np),2),:,np));
  
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
  
  nAcc = size(n_cor,1); % subdomain corner size
  
  nArr = size(n_rem,1); % subdomain remaining nodes size
  
  
  %%% Stiffness/Diffusion Coefficient
size_gamma = nAgg*n_outpce;
size_intr = nAii*n_outpce;
size_rem = nArr*n_outpce;
size_cor = nAcc*n_outpce;

Ajk_gamma = zeros(size_gamma,size_gamma);
Ajk_intr = zeros(size_intr,size_intr);
Ajk_igamma = zeros(size_intr,size_gamma);
Ajk_gintr = zeros(size_gamma,size_intr);
Ajk_cc = zeros(size_cor,size_cor);
Ajk_rr = zeros(size_rem,size_rem);
Ajk_ci = zeros(size_cor,size_intr);
Ajk_ri = zeros(size_rem,size_intr);
Ajk_cr = zeros(size_cor,size_rem);


Mjk_gamma = zeros(size_gamma,size_gamma);
Mjk_intr = zeros(size_intr,size_intr);
Mjk_igamma = zeros(size_intr,size_gamma);
Mjk_cc = zeros(size_cor,size_cor);
Mjk_rr = zeros(size_rem,size_rem);
Mjk_ci = zeros(size_cor,size_intr);
Mjk_ri = zeros(size_rem,size_intr);
Mjk_cr = zeros(size_cor,size_rem);


% Force Vector
% % % bk_gamma = zeros(size_gamma,1);
% % % bk_intr = zeros(size_intr,1);
  

% Restriction Matrix - stochastic
R_row = nAgg*n_outpce;
R_col = num_interface*n_outpce;
R_sto = zeros(R_row,R_col);


% Rc_sto
Rc_sto = zeros(size_cor,size_gamma); 

%Rr_sto
Rr_sto = zeros(size_rem,size_gamma); 


%Bc_sto
nAcc_global = size(node_corner,1); % Global corner nodes size
size_corner_global = nAcc_global*n_outpce;
Bc_sto = zeros(size_cor,size_corner_global); 


% Scaling Matrix - stochastic
D_row = nAgg*n_outpce;
D_sto = zeros(D_row); 
    
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
    
     
          %% A_remaining_remaining  
       
    Ajk_rr((k-1)*nArr+1:(k-1)*nArr+nArr, (j-1)*nArr+1:(j-1)*nArr+nArr)...
        = Ajk_rr((k-1)*nArr+1:(k-1)*nArr+nArr, (j-1)*nArr+1:(j-1)*nArr+nArr) +...
                         Cijk.moment_PsiNorm(i,j,k) * A(n_rem,n_rem,i);
        
        
        %% A_corner_corner
    Ajk_cc((k-1)*nAcc+1:(k-1)*nAcc+nAcc, (j-1)*nAcc+1:(j-1)*nAcc+nAcc)...
        = Ajk_cc((k-1)*nAcc+1:(k-1)*nAcc+nAcc, (j-1)*nAcc+1:(j-1)*nAcc+nAcc) +...
                         Cijk.moment_PsiNorm(i,j,k) * A(n_cor,n_cor,i);    
       
        
        %% A_corner_interior
        
     Ajk_ci((k-1)*nAcc+1:(k-1)*nAcc+nAcc, (j-1)*nAii+1:(j-1)*nAii+nAii)...
        = Ajk_ci((k-1)*nAcc+1:(k-1)*nAcc+nAcc, (j-1)*nAii+1:(j-1)*nAii+nAii) +...
                         Cijk.moment_PsiNorm(i,j,k) * A(n_cor,n_interior,i);       
        
        
        %% A_remaining_interior
        
      Ajk_ri((k-1)*nArr+1:(k-1)*nArr+nArr, (j-1)*nAii+1:(j-1)*nAii+nAii)...
        = Ajk_ri((k-1)*nArr+1:(k-1)*nArr+nArr, (j-1)*nAii+1:(j-1)*nAii+nAii) +...
                         Cijk.moment_PsiNorm(i,j,k) * A(n_rem,n_interior,i);
               
        
        %% A_corner_remaining
        
      Ajk_cr((k-1)*nAcc+1:(k-1)*nAcc+nAcc, (j-1)*nArr+1:(j-1)*nArr+nArr)...
        = Ajk_cr((k-1)*nAcc+1:(k-1)*nAcc+nAcc, (j-1)*nArr+1:(j-1)*nArr+nArr) +...
                         Cijk.moment_PsiNorm(i,j,k) * A(n_cor,n_rem,i);       
        
        
        
        end
        
    end
    
    %%% Force vector is same for all PCE terms
% % %    bk_gamma((k-1)*nAgg+1:(k-1)*nAgg+nAgg) = b(n_gamma,1)* exp_pce(k);
% % %    bk_intr((k-1)*nAii+1:(k-1)*nAii+nAii) = b(n_interior,1)* exp_pce(k);
% % %    
   
   %%% Restriction Matrix assembly
   
   R_sto((k-1)*nAgg+1:k*nAgg,(k-1)*num_interface+1:k*num_interface) = R(1:nAgg,:,np);
   
   Rc_sto((k-1)*nAcc+1:k*nAcc,(k-1)*nAgg+1:k*nAgg) = Rc(1:nAcc,1:nAgg,np);
   
   Rr_sto((k-1)*nArr+1:k*nArr,(k-1)*nAgg+1:k*nAgg) = Rr(1:nArr,1:nAgg,np);
   
   Bc_sto((k-1)*nAcc+1:k*nAcc,(k-1)*nAcc_global+1:k*nAcc_global) = Bc(1:nAcc,1:nAcc_global,np);
   
   %%% Scaling Matrix for Lumped Preconditioner
   
   D_sto((k-1)*nAgg+1:k*nAgg,(k-1)*nAgg+1:k*nAgg) = D(1:nAgg,1:nAgg,np);
   
   
   Mjk_gamma((k-1)*nAgg+1:k*nAgg,(k-1)*nAgg+1:k*nAgg) = ...
                                  M(n_gamma,n_gamma);
%                               *sqrt(s_spectral.norm.norm_squared(k));
    
   Mjk_intr((k-1)*nAii+1:k*nAii,(k-1)*nAii+1:k*nAii) = ...
                                  M(n_interior,n_interior);
%                                                        
                              
   Mjk_igamma((k-1)*nAii+1:k*nAii,(k-1)*nAgg+1:k*nAgg) = ...
                                  M(n_interior,n_gamma);
                              
   
   Mjk_rr((k-1)*nArr+1:k*nArr,(k-1)*nArr+1:k*nArr) = ...
                                  M(n_rem,n_rem);
                                                        
   Mjk_cc((k-1)*nAcc+1:k*nAcc,(k-1)*nAcc+1:k*nAcc) = ...
                                  M(n_cor,n_cor);   
                              
   Mjk_ci((k-1)*nAcc+1:k*nAcc,(k-1)*nAii+1:k*nAii) = ...
                                  M(n_cor,n_interior);   
                              
                              
   Mjk_ri((k-1)*nArr+1:k*nArr,(k-1)*nAii+1:k*nAii) = ...
                                  M(n_rem,n_interior);                           
                              
      
    Mjk_cr((k-1)*nAcc+1:k*nAcc,(k-1)*nArr+1:k*nArr) = ...
                                  M(n_cor,n_rem);                               
                              
                              
end


Mjk_gintr = Mjk_igamma';
Mjk_rc = Mjk_cr';

%%%%%% Finding out K_T


K_T_gamma = Mjk_gamma/(beta*deltaT^2) + Ajk_gamma;

K_T_intr = Mjk_intr/(beta*deltaT^2) + Ajk_intr;

K_T_igamma = Mjk_igamma/(beta*deltaT^2) + Ajk_igamma;

K_T_gintr = Mjk_gintr/(beta*deltaT^2) + Ajk_gintr;


K_T_rr = Mjk_rr/(beta*deltaT^2) + Ajk_rr;

K_T_cc = Mjk_cc/(beta*deltaT^2) + Ajk_cc;

K_T_ri = Mjk_ri/(beta*deltaT^2) + Ajk_ri;

K_T_ci = Mjk_ci/(beta*deltaT^2) + Ajk_ci;

K_T_cr = Mjk_cr/(beta*deltaT^2) + Ajk_cr;


cd ./AbSto
foldername = sprintf('part_%d',np);
mkdir(foldername);
cd (foldername)

save('K_T_gg','K_T_gamma')
save('K_T_ii','K_T_intr')
save('K_T_ig','K_T_igamma')
save('K_T_gi','K_T_gintr')

save('K_T_cc','K_T_cc')
save('K_T_rr','K_T_rr')
save('K_T_ci','K_T_ci')
save('K_T_ri','K_T_ri')
save('K_T_cr','K_T_cr')


save('Mjk_gg','Mjk_gamma')
save('Mjk_ii','Mjk_intr')
save('Mjk_ig','Mjk_igamma')
save('Mjk_gi','Mjk_gintr')

save('Mjk_cc','Mjk_cc')
save('Mjk_rr','Mjk_rr')
save('Mjk_ci','Mjk_ci')
save('Mjk_ri','Mjk_ri')
save('Mjk_cr','Mjk_cr')




% save('Bk_g','bk_gamma')
% save('Bk_i','bk_intr')
save('R_pce','R_sto')
save('Rc_pce','Rc_sto')
save('Rr_pce','Rr_sto')
save('Bc_pce','Bc_sto')
save('D_pce','D_sto')

cd ../../

end


end