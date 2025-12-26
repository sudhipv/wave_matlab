% Function for matrix vector product


function [Q] = GetMatrixVectorProduct(numpartition,P)

cd ./AbSto

Q_global = zeros;

for i=1:numpartition
    
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


Rmat = load('R_pce.mat','R_sto');
R_sto = Rmat.R_sto;

p_s = R_sto*P;
u1 = Ajk_ig * p_s;
u2 = Ajk_ii\u1;
u3 = Ajk_gi*u2;
u4 = Ajk_gg*p_s;
Q_local = u4 - u3;
Q_global = Q_global + R_sto'*Q_local;

cd ../

end

Q = Q_global;

cd ../
end
    
