



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%% SSFEM - Domain Decomposition - Acoustic Wave Propagation Problem - Random Variable %

% Sudhi Sharma P V %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%GmeshExtractor
% Given the file name, 
% 1)(p_f,e_f,t_f)- extracted points, edges, and triangles for different domains. 
% 2) GLnode - Global nodes by partition
% 3) npart - Number of partitions
% 4) point_exact, triangle_exact - Used for plotting exact solutions through PDE Surf
% 5) node_gamma, node_interior - Boundary and interior nodes by partition -"RENAMED"
% 6) node_interfacewhole, node_interiorwhole - Full Interface nodes for
% domain...interior nodes for the whole domain by partition- "NOT RENAMED"

clearvars;

delete('./AbDet/*')
if exist('AbSto')
rmdir AbSto s
end


[p_f,e_f,t_f,GLnode,npart,point_exact,triangle_exact,node_gamma,node_interior,...
    node_interfacewhole,node_interiorwhole,R]...
= ExtractGmshPart_puff('./square.msh');


% Check for interface and interior decomposing

num_inter = size(nonzeros(unique(node_interiorwhole)),1);
num_interface = size(node_interfacewhole,1);

node_g1 = sort(vertcat(node_interfacewhole,nonzeros(unique(node_interiorwhole))));
node_g2 = nonzeros(unique(GLnode));

node_total = size(node_g1,1);

if(any(node_g1-node_g2))
    disp("error")
end

%%%%....Stochastic Parameters......%%%%%%%%%

ord_in  = 7;
ord_out = 10;
dim     = 1;

%%%% First index order of expansion and second dimension
%%% Can use highest order for generating the multiindex and norms
str5    = '0100001';
str6    = '0100001';
% Inherent Gaussian Properties
mu_g = 0;
sigma_g = 0.1;

% Correlation Length
b = 1;
% Half the domain length
a = 0.5;

% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + dim)/(factorial(ord_in)*factorial(dim))
% Number of Spectral Terms in expansion of putput random coefficient
n_outpce = factorial(ord_out + dim)/(factorial(ord_out)*factorial(dim))


% num_spectral : determines the kind of expansion employed for input
% parameter
%%%% 1 ........ Gaussian Random Variable
%%%% 2 ........ Log Normal Random Variable
%%%% 3 ........ Gaussian Process % Not yet active
%%%% 4 ........ LogNormal Process

num_spectral = 4;

if(num_spectral >2)
    GetLambdaSymmetric_2D(b,a,dim);
    file = load('lambda_2D.mat');
else
    file = 0;
end


filename_mi = strcat('./../misc/cijk_ord_dim/mindex',str6);
filename_norm = strcat('./../misc/cijk_ord_dim/norm_squared',str5);

s_spectral = struct;
s_spectral.num_spectral = num_spectral;
s_spectral.n_inpce = n_inpce;
s_spectral.n_outpce = n_outpce;
s_spectral.dim = dim;
s_spectral.ord_in = ord_in;
s_spectral.ord_out = ord_out;
s_spectral.mu_g = mu_g;
s_spectral.sigma_g = sigma_g;
s_spectral.mindex = load(filename_mi);
s_spectral.norm = load(filename_norm);


%.............................%

%%%%%%%%%%%%%%%%%%%%%%
% New-Mark Beta Scheme

T = 1;

gamma = 0.5;

beta = 0.25;

deltaT = 0.01;

count = ceil(T/deltaT);


%%%%%%% Assembling and saving all the subdomain level matrices for each PC
%%%%%%% term %%%%%%%%%%%%

s_gamma_max = 1;
s_interior_max = 1;

for i=1:npart

p = p_f(any(p_f(:,1,i),2),2:3,i)';
e = e_f(:,:,i)';
t = t_f(any(t_f(:,1:3,i),2),:,i)';


 n_gamma = unique(node_gamma(any(node_gamma(:,1,i),2),:,i));
 n_interior = unique(node_interior(any(node_interior(:,1,i),2),:,i));
 
 if(size(n_gamma,1) > s_gamma_max)
     s_gamma_max = size(n_gamma,1);
 end
 
 if(size(n_interior,1) > s_interior_max)
     s_interior_max = size(n_interior,1);
 end
    
A = zeros(size(p,2),size(p,2),n_inpce);
K_T = zeros(size(p,2),size(p,2),n_inpce);

% Assemble Mass matrix
M = AssembleMatrix_ss(p, e, t, 2, 'Acoustic_ss', [], 0,0,0,0);

for ipce = 1:n_inpce


% Assemble Stiffness matrix
A(:,:,ipce) = AssembleMatrix_ss(p, e, t, 1, 'Acoustic_ss', [], 0,s_spectral,ipce,file);


%  K_T(:,:,ipce) = (M*sqrt(s_spectral.norm.norm_squared(ipce))/(beta*deltaT^2))+ A(:,:,ipce); 


end


Aname = sprintf('A_%d.mat',i);
Mname = sprintf('M_%d.mat',i);
% bname = sprintf('b_0%d.mat',i);
cd ./AbDet
save(Aname,'A');
save(Mname,'M');
% save(bname,'b');
cd ../

end

%%% Multiplication Tensor
filename = strcat('./../misc/cijk_ord_dim/Cijk',str5);
Cijk = load(filename);
Cijk.moment_PsiNorm;


Assemble_Stochastic(s_spectral,Cijk,npart,node_gamma,node_interior,num_interface,R,beta,deltaT);


%%%%%%%%%%%%%%%%%%%%%%  TIME LOOP STARTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Space coordinates %%%%

% uini = zeros(length(node_g2),1);
p_n = zeros(length(node_g2)*n_outpce,1);
a_n = zeros(length(node_g2)*n_outpce,1);


x = point_exact(1,:)';

y = point_exact(2,:)';

m = 2;
n = 1;

%%% Initial Conditions %%%%%%%%%%
uini = sin(m.*pi.*x).*sin(n.*pi.*y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


U_wave = zeros(length(node_g2)*n_outpce,count+1); %%% Total Solution at all time step

U = zeros(size(nonzeros(unique(GLnode)),1)*n_outpce,1); %%%% Solution at one timestep

u_n = zeros(length(node_g2)*n_outpce,1);
u_n(1:length(node_g2),:) = uini;

U_gamma_time = zeros(length(node_interfacewhole)*n_outpce,count+1);

U_gamma = zeros(length(node_interfacewhole)*n_outpce,1);

time = zeros(count+1,1);

NBcount =1;

U_wave(:,NBcount) = u_n;

exp_pce = zeros(1,n_outpce);
exp_pce(1,1) = 1;


while time<=T
    
time(NBcount+1)=time(NBcount)+deltaT;

NBcount = NBcount+1;

S_global = 0;

G_global = 0;


for i=1:npart


p = p_f(any(p_f(:,1,i),2),2:3,i)';
e = e_f(:,:,i)';
t = t_f(any(t_f(:,1:3,i),2),:,i)';

node_ini = p_f(any(p_f(:,1,i),2),1,i)';

n_gamma = unique(node_gamma(any(node_gamma(:,1,i),2),:,i));
n_interior = unique(node_interior(any(node_interior(:,1,i),2),:,i));

ng = length(n_gamma);
ni = length(n_interior);

% % % % Assemble vector
b = AssembleVector(p, e, t, 4, 'Acoustic_ss', [], time(NBcount));



% Retrieving saved sub-domain matrices
cd ./AbSto
foldername = sprintf('part_%d',i);
cd (foldername)

Mmat = load('Mjk_gg.mat','Mjk_gamma');
 Mjk_gg = Mmat.Mjk_gamma;
Mmat = load('Mjk_ii.mat','Mjk_intr');
 Mjk_ii = Mmat.Mjk_intr;

Mmat = load('Mjk_ig.mat','Mjk_igamma');
 Mjk_ig = Mmat.Mjk_igamma;
Mmat = load('Mjk_gi.mat','Mjk_gintr');
 Mjk_gi = Mmat.Mjk_gintr;
 
Rmat = load('R_pce.mat','R_sto');
R_sto = Rmat.R_sto;

% fileA = matfile(Aname);
% fileb = matfile(Bname);
% fileM = matfile(Mname);
% fileC = matfile(Cname);

% A = fileA.K_T;
% b = fileb.b;
% M_s = fileM.M;
% C = fileC.C;
cd ../../

nNsub = length(node_ini);

un_sub = zeros(nNsub*n_outpce,1);
pn_sub = zeros(nNsub*n_outpce,1);
an_sub = zeros(nNsub*n_outpce,1);

un_interface = zeros(num_interface*n_outpce,1);
pn_interface = zeros(num_interface*n_outpce,1);
an_interface = zeros(num_interface*n_outpce,1);

un_interior = zeros(ni*n_outpce,1);
pn_interior = zeros(ni*n_outpce,1);
an_interior = zeros(ni*n_outpce,1);

for nk = 1:n_outpce
    
un_sub((nk-1)*nNsub+1:nk*nNsub) = u_n(node_ini+(nk-1)*length(node_g2),:);
pn_sub((nk-1)*nNsub+1:nk*nNsub) = p_n(node_ini+(nk-1)*length(node_g2),:);
an_sub((nk-1)*nNsub+1:nk*nNsub) = a_n(node_ini+(nk-1)*length(node_g2),:);

un_interface((nk-1)*num_interface+1:nk*num_interface) = u_n(node_interfacewhole+(nk-1)*length(node_g2),:);
pn_interface((nk-1)*num_interface+1:nk*num_interface) = p_n(node_interfacewhole+(nk-1)*length(node_g2),:);
an_interface((nk-1)*num_interface+1:nk*num_interface) = a_n(node_interfacewhole+(nk-1)*length(node_g2),:);

un_interior((nk-1)*ni+1:nk*ni) = un_sub(n_interior+(nk-1)*nNsub,:);
pn_interior((nk-1)*ni+1:nk*ni) = pn_sub(n_interior+(nk-1)*nNsub,:);
an_interior((nk-1)*ni+1:nk*ni) = an_sub(n_interior+(nk-1)*nNsub,:);


end


%%%%%% Assembling the Linear System to be solved for Implicit Method - for
%%%%%% Each Subdomain %%%%%% 

fun_Un_An_interior = ((un_interior + deltaT*pn_interior)/(deltaT^2 * beta)) ...
                       + ((1-2*beta)*an_interior)/(2*beta);


un_sub_g = R_sto*un_interface;
pn_sub_g = R_sto*pn_interface;
an_sub_g = R_sto*an_interface;

fun_Un_An_interface = ((un_sub_g + deltaT*pn_sub_g)/(deltaT^2 * beta)) + ((1-2*beta)*an_sub_g)/(2*beta);


M_n_interior = Mjk_ii * fun_Un_An_interior + Mjk_ig*fun_Un_An_interface;


M_n_interface = Mjk_gi * fun_Un_An_interior + Mjk_gg*fun_Un_An_interface;

%%% Force vector is same for all PCE terms
bk_gamma = zeros(ng*n_outpce,1);
bk_intr = zeros(ni*n_outpce,1);

for kk = 1:n_outpce
    
   bk_gamma((kk-1)*ng+1:(kk-1)*ng+ng,1) = b(n_gamma) * exp_pce(kk);
   bk_intr((kk-1)*ni+1:(kk-1)*ni+ni,1) = b(n_interior) * exp_pce(kk);
% % % %%%%%   
end

b_t_interior = bk_intr + M_n_interior ;
b_t_interface = bk_gamma + M_n_interface ;


cd ./AbSto
foldername = sprintf('part_%d',i);
cd (foldername)

save('Bk_i','b_t_interior')
save('Bk_g','b_t_interface')

cd ../../

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


U = solve_DirectDDM(node_gamma,node_interior,npart, s_spectral,GLnode,...
                    node_interfacewhole,num_inter,node_total);


%%%%%%%%%%%%%%%%% DDM solver ends %%%%%%%%%%%%%%%%%%%%%%%%%%

u_new = U;

U_wave(:,NBcount) = u_new;

[u_n,p_n,a_n] = Update_acceleration(u_n,p_n,a_n,deltaT,beta,gamma,u_new);

end


                
% Plot solution

nt = length(time);
np = length(node_g1);

%%% Mean Solution
for j =1:nt
    
    
     s1 = pdesurf(point_exact,triangle_exact,U_wave(1:np,j));
     axis([0 1 0 1 -1 1]);
     zlim([-1 1])
     hold on
     colorbar
     pause(0.05);
     
         
    if(j<count)
         delete(s1)
    end

end


%%%%%%% Wave plot of PC coefficients %%%%%%%%%%

t_pos = 2;
indexpc = 2;

for j =1:nt
    
    
     s2 = pdesurf(point_exact,triangle_exact,U_wave(np*(indexpc-1)+1:np*(indexpc),j));
     axis([0 1 0 1 -1 1]);
     zlim([-1 1])
     hold on
     colorbar
     pause(0.05);
     
         
    if(j<count)
         delete(s2)
    end

end



%%%%%%%%%%% PC coefficients at partcular time step %%%%%%%%

index_coeff = size(U_wave,1)/n_outpce;

for j = 1:n_outpce
figure(j); clf
pdesurf(point_exact,triangle_exact,U_wave((j-1)*index_coeff+1:j*index_coeff,t_pos))
shading faceted
if(j==1)
title('Mean vector of solution expansion')
else
    titlename = sprintf('%d th PC coefficient of solution using %d RV with %drd order expansion',j,dim,ord_out);
    title(titlename)
end

end



% Calculation of standard deviation

U_std = zeros(size(point_exact,2),count+1);

for j = 2:n_outpce
U_std(:,:) = U_std(:,:) + U_wave((j-1)*node_total+1:(j-1)*node_total+node_total,:).^2;
end

U_std = sqrt(U_std);


% % % % % % % for j =1:nt
% % % % % % %     
% % % % % % %     
% % % % % % %      s3 = pdesurf(point_exact,triangle_exact,U_std(:,j));
% % % % % % %      axis([0 1 0 1 -1 1]);
% % % % % % %      zlim([-1 1])
% % % % % % %      hold on
% % % % % % %      colorbar
% % % % % % %      pause(0.05);
% % % % % % %      
% % % % % % %          
% % % % % % %     if(j<count)
% % % % % % %          delete(s3)
% % % % % % %     end
% % % % % % % 
% % % % % % % end




%%%% Mean at a point

pos = 400;

figure(34)
plot(time,U_wave(pos,:))
grid on


figure(35)
plot(time,U_std(pos,:))
grid on










