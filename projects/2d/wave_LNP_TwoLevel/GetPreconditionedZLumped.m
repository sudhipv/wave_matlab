% Function for pre-conditioning


function [Z_o] = GetPreconditionedZLumped(numpartition,r_gamma)

cd ./AbSto

Z_o = zeros;

for i=1:numpartition
    
foldername = sprintf('part_%d',i);
    cd (foldername)
    
 Amat = load('K_T_cc.mat');
 Ajk_cc = Amat.K_T_cc;

Rmat = load('Bc_pce.mat','Bc_sto');
Bc_sto = Rmat.Bc_sto;


%% Without scaling -Extended Lumped
r = Bc_sto*r_gamma;
z = Ajk_cc\r;
Z_o = Z_o + Bc_sto'* z;

cd ../
end

cd ../
end
    
