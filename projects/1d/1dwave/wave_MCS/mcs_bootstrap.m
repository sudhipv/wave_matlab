

tic

%%% 1D mesh %%%

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



%%%%%%%% MCS Parameters %%%%%%%%%%%%%%%%%%%%%%%

ns_orig = 2000;
xi_org = randn(1,ns_orig);

btstrap = 1000;






mu_g = 0;
sigma_g = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

         M(i:i+1,i:i+1) =M(i:i+1,i:i+1) + m_e;


end


%%%%%%Assemble Stiffness Matrix%%%%%%

EA = 5;
k_e = (EA/l)* [ 1 -1;
                -1 1];

K = zeros(nN,nN);

for i = 1:ne

        K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + k_e;

end

K(1,1) = K(1,1) + 1* 10^(15);

m_alpha = 0.01;
k_alpha = 0.001;
 
C = m_alpha *M + k_alpha *K;

%%% Force Vector

% May have to check for nyquist frequency condition
% deltaT < bt * 2*PI/(w max- frequency of force) 

wmax = 100;
f = 1000*sin(wmax.*t); 
B = zeros(nN,count+1);
% 
B(nN,:) = f;


t_nyqst = 0.5 * 2* pi/(wmax);
fprintf("delta T nyquist %2.2f\n",t_nyqst)
t_cfl = h/(sqrt(EA/rho));
fprintf("delta T CFL %4.4f\n",t_cfl)

fprintf("delta T used %4.3f\n",deltaT)


%%% Boundary Condition - Fixed at left end (node 1) adding a big penalty
%%% term

%  K(1,1) = K(1,1) + 1* 10^(7);


boot_mean = zeros(nN,count+1,btstrap);
boot_sd = zeros(nN,count+1,btstrap);

for bb = 1: btstrap
    
    
ns = 100;

xi = datasample(xi_org,100);

%%%% MCS Sampling %%%%%%

u_MCS = zeros(nN,count+1,ns);


for ss = 1:ns
    
%%%% Initial Conditions %%%%%%%%%%
uini = zeros(nN,1);

p_n = zeros(nN,1);

a_n = zeros(nN,1);

u_n = uini;

U = zeros(nN,count+1);

time = zeros(count+1,1);

NBcount = 1;

U(:,NBcount) = u_n;


while NBcount<=count
    
    
time(NBcount+1)=time(NBcount)+deltaT;

NBcount = NBcount+1;

fun_Un_An = ((u_n + deltaT*p_n)/(deltaT^2 * beta)) + ((1-2*beta)*a_n)/(2*beta);

M_n = M * fun_Un_An;

C_n = C*gamma*deltaT*fun_Un_An - (C*p_n) - (C*(1-gamma)*a_n*deltaT);

F_T = B(:,NBcount) + M_n + C_n;


K_ss = assemble_stiffness(ne,nN,K,xi(1,ss),mu_g,sigma_g,x,l);
            
   
K_T = M/(beta*deltaT^2) + C*gamma/(beta*deltaT) + K_ss ;


u_new = K_T\F_T;  
    

U(:,NBcount) = u_new;

[u_n,p_n,a_n] = Update_acceleration(u_n,p_n,a_n,deltaT,beta,gamma,u_new);


end


%%% 1 realization of MCS sampling 

u_MCS(:,:,ss) = U;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% % % Mean Solution
U_wave = sum(u_MCS,3)/ns;


% % % Standard Deviation of Solution
sdum = zeros(nN,count+1);

 for j = 1:ns
   
     sdum = sdum + (u_MCS(:,:,j) - U_wave).^2;
     
     
 end
 
 U_sd = sdum/(ns-1);
 
 U_sd = sqrt(U_sd);
 
 
 boot_mean(:,:,bb) = U_wave;
 boot_sd(:,:,bb) = U_sd;

end




toc






