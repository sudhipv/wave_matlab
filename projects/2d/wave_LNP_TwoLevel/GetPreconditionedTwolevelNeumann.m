
% Function to compute Coarse grid residue



function Z_coarse = GetPreconditionedTwolevelNeumann(npart,r_gamma_o,tolc)

d_global = 0;
z_global = 0;

cd ./AbSto

for i=1:npart
    

% % % n_c = unique(node_corner_part(any(node_corner_part(:,1,i),2),:,i));
% % % n_r = unique(node_remaining_part(any(node_remaining_part(:,1,i),2),:,i));
% % % n_i = unique(node_interior(any(node_interior(:,1,i),2),:,i));
% % % n_gamma = unique(node_gamma(any(node_gamma(:,1,i),2),:,i));
% % % 
% % % sg = size(n_gamma,1);
% % % sr = size(n_r,1);
% % % sc = size(n_c,1);


% Retrieving saved sub-domain matrices

  foldername = sprintf('part_%d',i);
    cd (foldername)
    
 Amat = load('K_T_ii.mat','K_T_intr');
 Ajk_ii = Amat.K_T_intr;   
    
% %  Amat = load('Ajk_cc.mat');
% %  Ajk_cc = Amat.Ajk_cc;
Amat = load('K_T_rr.mat','K_T_rr');
 Ajk_rr = Amat.K_T_rr;
Amat = load('K_T_ci.mat','K_T_ci');
 Ajk_ci = Amat.K_T_ci;
Amat = load('K_T_ri.mat','K_T_ri');
 Ajk_ri = Amat.K_T_ri;
Amat = load('K_T_cr.mat','K_T_cr');
 Ajk_cr = Amat.K_T_cr;  
 
 Ajk_ir = Ajk_ri';
% %  Ajk_ic = Ajk_ci';
% %  Ajk_rc = Ajk_cr';

Rmat = load('R_pce.mat','R_sto');
R_sto = Rmat.R_sto;

Rmat = load('Rr_pce.mat','Rr_sto');
Rr_sto = Rmat.Rr_sto;

Rmat = load('Rc_pce.mat','Rc_sto');
Rc_sto = Rmat.Rc_sto;

Rmat = load('Bc_pce.mat','Bc_sto');
Bc_sto = Rmat.Bc_sto;

Dmat = load('D_pce.mat','D_sto');
D_sto = Dmat.D_sto;


s_rr = Ajk_rr - (Ajk_ri * (Ajk_ii\Ajk_ir));

s_cr = Ajk_cr - (Ajk_ci * (Ajk_ii\Ajk_ir));

% s_rc = Ajk_rc - (Ajk_ri * (Ajk_ii\Ajk_ic));


r_gamma_local = D_sto * R_sto * r_gamma_o;

f_r = Rr_sto * r_gamma_local;

f_c = Rc_sto * r_gamma_local;

v1 = s_rr\f_r;

d_local = f_c - (s_cr * v1);

d_global = d_global + Bc_sto'*d_local;


save('f_r','f_r');

cd ../

end

cd ../
%%%%%%%%%%%%%%%%%%%

% Initial assumption
U_o = zeros(length(d_global),1);
r_gamma_o = d_global;
z_c = DoPCGMLumpedTwoLevel(npart,U_o,r_gamma_o,tolc);

%%%%%%%%%%%%%%%%%%%

cd ./AbSto

for i=1:npart
     
% Retrieving saved sub-domain matrices

foldername = sprintf('part_%d',i);
    cd (foldername)
    
 Amat = load('K_T_ii.mat','K_T_intr');
 Ajk_ii = Amat.K_T_intr;   
    
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

Rmat = load('R_pce.mat','R_sto');
R_sto = Rmat.R_sto;

Rmat = load('Rr_pce.mat','Rr_sto');
Rr_sto = Rmat.Rr_sto;

Rmat = load('Rc_pce.mat','Rc_sto');
Rc_sto = Rmat.Rc_sto;

Rmat = load('Bc_pce.mat','Bc_sto');
Bc_sto = Rmat.Bc_sto;

Dmat = load('D_pce.mat','D_sto');
D_sto = Dmat.D_sto;

fmat = load('f_r','f_r');
f_r = fmat.f_r;


s_rr = Ajk_rr - (Ajk_ri * (Ajk_ii\Ajk_ir));

s_rc = Ajk_rc - (Ajk_ri * (Ajk_ii\Ajk_ic));


z_c_local = Bc_sto * z_c;

v2 = f_r - s_rc * z_c_local;

z_r_local = s_rr\v2;

z_local = Rr_sto'*z_r_local + Rc_sto'*z_c_local;

z_global = z_global +  R_sto' * D_sto * z_local;

cd ../

end

Z_coarse = z_global;
cd ../
end