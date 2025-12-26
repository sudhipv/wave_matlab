


function [U, iter, nres] = solve_NeumannTwoLevel(U_gamma_o,r_gamma_o,node_gamma,node_interior,...
         npart, s_spectral,GLnode,node_interfacewhole,num_interior,node_total,tol,tolc,NBcount)
                      
% %   Neumann -Two Level                     
                      


[U_gamma,iter,nres] = DoPCGM_NeumannTwoLevel(npart,U_gamma_o,r_gamma_o,tol,tolc);


U_gamma_time(:,NBcount) = U_gamma;


node_global = nonzeros(unique(GLnode));

n_outpce = s_spectral.n_outpce;


% Constructing the A and b Matrices after getting interface solution. Can
% be avoided by saving them and reloading at this stage.
cd ./Absto

for j=1:npart
  
 foldername = sprintf('part_%d',j);
    cd (foldername)
    

Amat = load('K_T_ii.mat','K_T_intr');
 Ajk_ii = Amat.K_T_intr;
Amat = load('K_T_ig.mat','K_T_igamma');
 Ajk_ig = Amat.K_T_igamma;
Amat = load('Bk_i.mat','b_t_interior');
Bjk_i = Amat.b_t_interior;   

Rmat = load('R_pce.mat','R_sto');
R_sto = Rmat.R_sto;    
     

n_gamma = unique(node_gamma(any(node_gamma(:,1,j),2),:,j));
n_interior = unique(node_interior(any(node_interior(:,1,j),2),:,j));

U_interior  = Ajk_ii\(Bjk_i - (Ajk_ig * R_sto * U_gamma));

%%% Interior to Global Nodes by PCE

for ni = 1:size(n_interior,1)
 node_in = GLnode(n_interior(ni),j);
 U(node_in:node_total:node_total*n_outpce,1) = U_interior(ni:size(n_interior,1):...
                                               size(n_interior,1)*n_outpce);
end

cd ../

end


%%% Whole solution found %%%%%


%%%% Arranging the nodes


% Interface to Global Nodes by PCE
num_interface = node_total - num_interior;

for nr = 1:size(node_interfacewhole,1)
    
     node_g = node_global(node_interfacewhole(nr));
     U(node_g:node_total:node_total*n_outpce,1) = U_gamma(nr:num_interface:...
                                                            num_interface*n_outpce);
     
end

cd ../                     

end
