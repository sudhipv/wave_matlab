



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%% SSFEM - Acoustic Wave Propagation Problem - Random Variable %

% Sudhi Sharma P V %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;


% Run ExtractGmsh.m

p = load('points.txt')';
e =load('edges.txt')';
t =load('triangles.txt')';


nN = size(p,2);

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
sigma_g = 0.1;%%% Change folder name at end to save file

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

num_spectral = 2;

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

deltaT = 0.005;

count = ceil(T/deltaT);

    
A = zeros(nN,nN,n_inpce);

% Assemble Mass matrix
M = AssembleMatrix_ss(p, e, t, 2, 'Acoustic_ss', [], 0,0,0,0);

for ipce = 1:n_inpce


% Assemble Stiffness matrix
A(:,:,ipce) = AssembleMatrix_ss(p, e, t, 1, 'Acoustic_ss', [], 0,s_spectral,ipce,file);


%  K_T(:,:,ipce) = (M*sqrt(s_spectral.norm.norm_squared(ipce))/(beta*deltaT^2))+ A(:,:,ipce); 


end

%%% For Gaussian Pulse
C = 0.5 .* M + 0.01*A(:,:,1);

%%% Intrusive SSFEM solve %%%%
A_jk = zeros(size(A,1)*n_outpce);
%%% Multiplication Tensor
filename = strcat('./../misc/cijk_ord_dim/Cijk',str5);
Cijk = load(filename);
Cijk.moment_PsiNorm;
% Cijk.moment_ND;
n_A = size(A,1);
n_M = size(M,1);

 
% Assembly block struture
for k = 1: n_outpce
    for j = 1:n_outpce
        
        for i= 1: n_inpce
     
    A_jk((k-1)*n_A+1:(k-1)*n_A+n_A, (j-1)*n_A+1:(j-1)*n_A+n_A)...
        = A_jk((k-1)*n_A+1:(k-1)*n_A+n_A, (j-1)*n_A+1:(j-1)*n_A+n_A) + Cijk.moment_PsiNorm(i,j,k) * A(:,:,i);
    
        end
        
    end
    
    M_jk((k-1)*n_M+1:k*n_M,(k-1)*n_M+1:k*n_M) = M(:,:);
    C_jk((k-1)*n_M+1:k*n_M,(k-1)*n_M+1:k*n_M) = C(:,:);
    
end


K_T = M_jk/(beta*deltaT^2) + (C_jk*gamma/(beta*deltaT))+ A_jk;

%%%%%%%%%%%%%%%%%%%%%%  TIME LOOP STARTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Space coordinates %%%%

x = p(1,:)';

y = p(2,:)';

m = 2;
n = 1;

%%% Initial Conditions %%%%%%%%%%
% uini = sin(m.*pi.*x).*sin(n.*pi.*y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


U_wave = zeros(nN*n_outpce,count+1); %%% Total Solution at all time step

U = zeros(nN*n_outpce,1); %%%% Solution at one timestep

u_n = zeros(nN*n_outpce,1);
p_n = zeros(nN*n_outpce,1);
a_n = zeros(nN*n_outpce,1);


% u_n(1:nN) = uini;

time = zeros(count+1,1);

NBcount =1;

U_wave(:,NBcount) = u_n;

exp_pce = zeros(1,n_outpce);
exp_pce(1,1) = 1;


while time<=T
    
time(NBcount+1)=time(NBcount)+deltaT;

NBcount = NBcount+1;


fun_Un_An = ((u_n + deltaT*p_n)/(deltaT^2 * beta)) + ((1-2*beta)*a_n)/(2*beta);

M_n = M_jk * fun_Un_An;

% % % % Assemble vector
b = AssembleVector(p, e, t, 4, 'Acoustic_ss', [], NBcount);


C_n = C_jk*(p_n + deltaT*(1-gamma)*a_n - (deltaT * gamma * fun_Un_An));


%%% Force vector is same for all PCE terms
b_jk = zeros(nN*n_outpce,1);

for kk = 1:n_outpce
    
   b_jk((kk-1)*nN+1:(kk-1)*nN+nN,1) = b * exp_pce(kk);

end


b_t = b_jk + M_n - C_n;


u_new = K_T\b_t;


U_wave(:,NBcount) = u_new;

[u_n,p_n,a_n] = Update_acceleration(u_n,p_n,a_n,deltaT,beta,gamma,u_new);

end

               
% Plot solution

nt = length(time);

%%% Mean Solution
for j =1:nt
    
    
     s1 = pdesurf(p,t,U_wave(1:nN,j));
     axis([0 1 0 1 -0.5 0.5]);
     view(-45,30)
%      zlim([-1 1])
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
    
    
     s2 = pdesurf(p,t,U_wave(nN*(indexpc-1)+1:nN*(indexpc),j));
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


for j = 1:n_outpce
figure(j); clf
pdesurf(p,t,U_wave((j-1)*nN+1:j*nN,t_pos))
shading faceted
if(j==1)
title('Mean vector of solution expansion')
else
    titlename = sprintf('%d th PC coefficient of solution using %d RV with %drd order expansion',j,dim,ord_out);
    title(titlename)
end

end



% Calculation of standard deviation

U_std = zeros(size(p,2),count+1);

for j = 2:n_outpce
U_std(:,:) = U_std(:,:) + U_wave((j-1)*nN+1:(j-1)*nN+nN,:).^2;
end

U_std = sqrt(U_std);


for j =1:nt
    
    
     s3 = pdesurf(p,t,U_std(:,j));
     axis([0 1 0 1 -1 1]);
     zlim([-1 1])
     hold on
     colorbar
     pause(0.05);
     
         
    if(j<count)
         delete(s3)
    end

end



%%%% Mean at a point

pos = 246;

figure(24)
plot(time,U_wave(pos,:))
grid on


figure(25)
plot(time,U_std(pos,:))
grid on



% % cd ./results/7_10_sd_3_nodamp
% % save('U_intrusive.mat','U_wave')
% % save('U_intr_std.mat','U_std')
% % save('p.mat','p')
% % save('t.mat','t')
% % save('time.mat','time')





