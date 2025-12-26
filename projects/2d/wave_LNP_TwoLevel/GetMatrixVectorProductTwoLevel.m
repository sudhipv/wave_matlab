% Function for matrix vector product


function [Q] = GetMatrixVectorProductTwoLevel(numpartition,P)

Q_global = zeros;

cd ./AbSto

for i=1:numpartition
    

foldername = sprintf('part_%d',i);
    cd (foldername)
    
 Amat = load('K_T_ii.mat','K_T_intr');
 Ajk_ii = Amat.K_T_intr;   
    
Amat = load('K_T_cc.mat');
Ajk_cc = Amat.K_T_cc;
Amat = load('K_T_rr.mat','K_T_rr');
 Ajk_rr = Amat.K_T_rr;
Amat = load('K_T_ci.mat','K_T_ci');
 Ajk_ci = Amat.K_T_ci;
Amat = load('K_T_ri.mat','K_T_ri');
 Ajk_ri = Amat.K_T_ri;
Amat = load('K_T_cr.mat','K_T_cr');
 Ajk_cr = Amat.K_T_cr;  
 
 Ajk_ir = Ajk_ri';
 Ajk_ic = Ajk_ci';
 Ajk_rc = Ajk_cr';

Rmat = load('Bc_pce.mat','Bc_sto');
Bc_sto = Rmat.Bc_sto;

s_rr = Ajk_rr - (Ajk_ri * (Ajk_ii\Ajk_ir));

s_cr = Ajk_cr - (Ajk_ci * (Ajk_ii\Ajk_ir));

s_rc = Ajk_rc - (Ajk_ri * (Ajk_ii\Ajk_ic));

s_cc = Ajk_cc - (Ajk_ci * (Ajk_ii\Ajk_ic));

p_s = Bc_sto*P;
u1 = s_rc * p_s;
u2 = s_rr\u1;
u3 = s_cr*u2;
u4 = s_cc*p_s;
Q_local = u4 - u3;
Q_global = Q_global + Bc_sto'*Q_local;

cd ../
end

Q = Q_global;
cd ../
end
    
