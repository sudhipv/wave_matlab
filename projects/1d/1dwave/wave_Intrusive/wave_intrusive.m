



%%% 1D mesh %%%
clc
clear
%%% Element Size
h = 0.01;

%%% Meshing
x = 0:h:1;

% Number of points
nN = length(x);

% Number of elements
ne = nN-1;


%%%%%%%%%%%%%%

T = 1;

gamma = 0.5;

beta = 0.25;

deltaT = 0.002;

count = ceil(T/deltaT);

t = 0:deltaT:T;


%%%%....Stochastic Parameters......%%%%%%%%%

ord_in  = 10;
ord_out = 15;
dim     = 1;

%%%% First index order of expansion and second dimension
%%% Can use highest order for generating the multiindex and norms
str_orddim    = '0150001';
str_mindex = '0150001';
norm = strcat('./../../misc/cijk_ord_dim/norm_squared',str_orddim);
mindex = strcat('./../../misc/cijk_ord_dim/mindex',str_mindex);
norm = load(norm)
mindex = load(mindex)

% Inherent Gaussian Properties
mu_g = 0;
sigma_g = 0.2; %%%% change folder name to save results at bottom

% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + dim)/(factorial(ord_in)*factorial(dim))
% Number of Spectral Terms in expansion of putput random coefficient
n_outpce = factorial(ord_out + dim)/(factorial(ord_out)*factorial(dim))



%%%%%%Assemble Mass Matrix%%%%%%
   
rho = 1;
area = 1;
l = h;
% l = p(5,1)- p(4,1)

m =    [2   1;
       1   2];
   
m_e = (1/6)*(rho * area * l)* m;


M = zeros(nN,nN);

for i = 1:ne

         M(i:i+1,i:i+1) = M(i:i+1,i:i+1) + m_e;


end


%%%%%%Assemble Stiffness Matrix%%%%%%

EA = 1;
EA_0 = 5;
k_e = (EA_0/l)* [ 1 -1;
                -1 1];

          
            
K = zeros(nN,nN,n_inpce);

K_damp = zeros(nN,nN);

%%% For finding the damping term
for i = 1:ne

        K_damp(i:i+1,i:i+1) = K_damp(i:i+1,i:i+1) + k_e;

end

 K_damp(1,1) = K_damp(1,1) + 1* 10^(15);
%%% Damping matrix

m_alpha = 0.01;
k_alpha = 0.001;
 
C = m_alpha *M + k_alpha*K_damp;

for ipce = 1:n_inpce  
    

% for i = 1:ne
    
       %%% For generalizing with random process
% %         EA_ss = getEA_stochastic(i,ipce,EA_0,l,x,mu_g,sigma_g,norm,mindex,dim);
% % 
% %         K(i:i+1,i:i+1,ipce) = K(i:i+1,i:i+1,ipce) + k_e*EA_ss;


% end

%%% Boundary Condition - Fixed at left end (node 1) adding a big penalty
%%% term

%  K(1,1,ipce) = K(1,1,ipce) + 1* 10^(7);


EA_ss = getEA_stochastic(i,ipce,EA_0,l,x,mu_g,sigma_g,norm,mindex,dim);


K(:,:,ipce) = K_damp(:,:) * EA_ss;

 
end


%%% Force Vector

% May have to check for nyquist frequency condition
% deltaT < bt * 2*PI/(w max- frequency of force) 

wmax = 100;
f = 1000*sin(wmax.*t); 
B = zeros(nN*n_outpce,count+1);
% 
B(nN,:) = f;

%%% Force vector have all other PC coefficients 0 since <B
%%% psi_k>/sqrt(<psi_k^2> is 1 for first term and 0 for all others


t_nyqst = 0.5 * 2* pi/(wmax);
fprintf("delta T nyquist %2.2f\n",t_nyqst)
t_cfl = h/(sqrt(EA_0/rho));
fprintf("delta T CFL %4.4f\n",t_cfl)

fprintf("delta T used %4.3f\n",deltaT)

%%%%% Stochastic Assembly

 % Assembly block struture
 
 %%% Intrusive SSFEM solve %%%%
K_jk = zeros(size(K,1)*n_outpce);
%%% Multiplication Tensor
filename = strcat('./../../misc/cijk_ord_dim/Cijk',str_orddim);
Cijk = load(filename);
Cijk.moment_PsiNorm;
n_A = size(K,1);
n_M = size(M,1);

for k = 1: n_outpce
    for j = 1:n_outpce
        
        for i= 1: n_inpce
     
    K_jk((k-1)*n_A+1:(k-1)*n_A+n_A, (j-1)*n_A+1:(j-1)*n_A+n_A)...
        = K_jk((k-1)*n_A+1:(k-1)*n_A+n_A, (j-1)*n_A+1:(j-1)*n_A+n_A) + Cijk.moment_PsiNorm(i,j,k) * K(:,:,i);
    
        end
        
    end
    
    M_jk((k-1)*n_M+1:(k-1)*n_M+n_M, (k-1)*n_M+1:(k-1)*n_M+n_M) = M(:,:);
%     sqrt(norm.norm_squared(k)) * 
   
   C_jk((k-1)*n_M+1:(k-1)*n_M+n_M, (k-1)*n_M+1:(k-1)*n_M+n_M) = C(:,:);

end

%%%%%%%%

K_T = M_jk/(beta*deltaT^2) + C_jk*gamma/(beta*deltaT)+ K_jk;


 
%%%% Initial Conditions %%%%%%%%%%
uini = zeros(nN*n_outpce,1);

p_n = zeros(nN*n_outpce,1);

a_n = zeros(nN*n_outpce,1);

u_n = uini;

U = zeros(nN*n_outpce,count+1);

time = zeros(count+1,1);

NBcount = 1;

U(:,NBcount) = u_n;
% f = 10;
% U(:,NBcount) = sin(2*pi*f.*time);


while NBcount<=count
    
time(NBcount+1)=time(NBcount)+deltaT;

NBcount = NBcount+1;

fun_Un_An = ((u_n + deltaT*p_n)/(deltaT^2 * beta)) + ((1-2*beta)*a_n)/(2*beta);

M_n = M_jk * fun_Un_An;

C_n = C_jk*gamma*deltaT*fun_Un_An - (C_jk*p_n) - (C_jk*(1-gamma)*a_n*deltaT);


F_T = B(:,NBcount) + M_n + C_n;
 
u_new = K_T\F_T;           

U(:,NBcount) = u_new;

[u_n,p_n,a_n] = Update_acceleration(u_n,p_n,a_n,deltaT,beta,gamma,u_new);


end

%%%%%%%%%%%%%%%%%%%%%%%
nt = length(time);

index_coeff = size(U,1)/n_outpce;
% Calculation of standard deviation

U_std = zeros(nN,count+1);

for kk = 1:nt

for j = 2:n_outpce
U_std(:,kk) = U_std(:,kk) + U((j-1)*index_coeff+1:(j-1)*index_coeff+index_coeff,kk).^2;
end

end

U_std = sqrt(U_std);


% Plot the results
plot_times = [100 150 250 350 500];


%%% Mean 

for i = 1:length(plot_times)
    
  figure(i)
  k = plot_times(i);
  plot(x,U(1:nN,k),'linewidth',2);
  grid on;
  axis([min(x) max(x) -25 25]);
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);              
  titlestring = ['TIME STEP = ',num2str(k), ' TIME = ',num2str(time(k)), 'second'];
  title(titlestring ,'fontsize',14);                            
  h=gca; 
  get(h,'FontSize');
  set(h,'FontSize',14);
  fh = figure(i);
  set(fh, 'color', 'white'); 
  
  
end

%%%% Standard Deviation

for i = 1:length(plot_times)
    
  figure(i+10)
  k = plot_times(i);
  plot(x,U_std(:,k),'linewidth',2);
  grid on;
  axis([min(x) max(x) 0 25]);
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);              
  titlestring = ['TIME STEP = ',num2str(k), ' TIME = ',num2str(time(k)), 'second'];
  title(titlestring ,'fontsize',14);                            
  h=gca; 
  get(h,'FontSize');
  set(h,'FontSize',14);
  fh = figure(i);
  set(fh, 'color', 'white'); 
  
  
end


%%%%% PCE coefficient

pc_num = 6;


for i = 1:length(plot_times)
    
  figure(i+20)
  k = plot_times(i);
  plot(x,U((pc_num-1)*nN+1:(pc_num-1)*nN+nN,k),'linewidth',2);
  grid on;
  axis([min(x) max(x) -3 3]);
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);              
  titlestring = ['TIME STEP = ',num2str(k), ' TIME = ',num2str(time(k)), 'second'];
  title(titlestring ,'fontsize',14);                            
  h=gca; 
  get(h,'FontSize');
  set(h,'FontSize',14);
  fh = figure(i);
  set(fh, 'color', 'white'); 
  
  
end

cd ./Results/sd_2
save('U_intrusive.mat','U')
save('U_std_intrusive.mat','U_std')
save('x_intrusive.mat','x')






%-------------------------------------------------------------------------%
% Movie for the travelling wave

for j = 1:nt              
  plot(x,U(1:nN,j),'linewidth',2);
  grid on;
  axis([min(x) max(x) -25 25]);
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);              
  titlestring = ['TIME STEP = ',num2str(j), ' TIME = ',num2str(time(j)), 'second'];
  title(titlestring ,'fontsize',14);                            
  h=gca; 
  get(h,'FontSize');
  set(h,'FontSize',14);
  fh = figure(length(plot_times)+1);
  set(fh, 'color', 'white'); 
  G=getframe;
            
end

movie(G,1,10)
%-------------------------------------------------------------------------%

%%% Code to save variables to be used in Comparing with Intrusive and NISP






