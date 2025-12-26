


function U = solve_DirectDDM(node_gamma,node_interior,npart, s_spectral,GLnode,...
                         node_interfacewhole,num_inter,node_total)
%%%% Function for solving Direct DDM 

%..............%

S_global = 0;

G_global = 0;

U = zeros(node_total*s_spectral.n_outpce,1);

cd ./AbSto


for i=1:npart
  
    foldername = sprintf('part_%d',i);
    cd (foldername)
    
 Amat = load('K_T_gg.mat');
 Ajk_gg = Amat.K_T_gamma;
Amat = load('K_T_ii.mat','K_T_intr');
 Ajk_ii = Amat.K_T_intr;
Amat = load('K_T_ig.mat','K_T_igamma');
 Ajk_ig = Amat.K_T_igamma;
Amat = load('K_T_gi.mat','K_T_gintr');
 Ajk_gi = Amat.K_T_gintr;
Amat = load('Bk_g.mat','b_t_interface');
 Bjk_g = Amat.b_t_interface;
Amat = load('Bk_i.mat','b_t_interior');
Bjk_i = Amat.b_t_interior;   

Rmat = load('R_pce.mat','R_sto');
R_sto = Rmat.R_sto;


    

n_gamma = unique(node_gamma(any(node_gamma(:,1,i),2),:,i));
n_interior = unique(node_interior(any(node_interior(:,1,i),2),:,i));


S_local = Ajk_gg - (Ajk_gi *(Ajk_ii\Ajk_ig));

G_local = Bjk_g - (Ajk_gi *(Ajk_ii\Bjk_i));
      

S_global = S_global + R_sto' * S_local * R_sto;

G_global = G_global + R_sto' * G_local;

cd ../

end

U_gamma = S_global\G_global;

node_global = nonzeros(unique(GLnode));

n_outpce = s_spectral.n_outpce;


% Constructing the A and b Matrices after getting interface solution. Can
% be avoided by saving them and reloading at this stage.


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
num_interface = node_total - num_inter;

for nr = 1:size(node_interfacewhole,1)
    
     node_g = node_global(node_interfacewhole(nr));
     U(node_g:node_total:node_total*n_outpce,1) = U_gamma(nr:num_interface:...
                                                            num_interface*n_outpce);
     
end

cd ../

end
