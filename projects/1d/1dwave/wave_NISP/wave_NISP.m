



%%% 1D mesh %%%
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

ord_in  = 7;
ord_out = 10;
dim     = 1;

%%%% First index order of expansion and second dimension
%%% Can use highest order for generating the multiindex and norms
str_orddim    = '0100001';
str_mindex = '0100001';
norm = strcat('./../../misc/cijk_ord_dim/norm_squared',str_orddim);
mindex = strcat('./../../misc/cijk_ord_dim/mindex',str_mindex);
norm = load(norm)
mindex = load(mindex)


% Number of Spectral Terms in expansion of input random coefficient
n_inpce = factorial(ord_in + dim)/(factorial(ord_in)*factorial(dim))
% Number of Spectral Terms in expansion of putput random coefficient
n_outpce = factorial(ord_out + dim)/(factorial(ord_out)*factorial(dim))


%%%%%%%% MCS Parameters %%%%%%%%%%%%%%%%%%%%%%%

ns = 1000;
xi = randn(1,ns);

mu_g = 0;
sigma_g = 0.1;
%%%%%%%%% NISP parameter %%%%%%%%%%%%%%%

psi(:,1) = ones(1,ns);
psi(:,2) = xi;
psi(:,3) = xi.^2-1;
psi(:,4) = xi.^3-3.*xi;
psi(:,5) = xi.^4-6*xi.^2+3;
psi(:,6) = xi.^5-10*xi.^3+15*xi;
psi(:,7) = xi.^6-15*xi.^4+45*xi.^2-15;
psi(:,8) = xi.^7-21*xi.^5+105*xi.^3-105*xi;
psi(:,9) = xi.^8-28*xi.^6+210*xi.^4-420*xi.^2+105;
psi(:,10) = xi.^9-36*xi.^7+378*xi.^5-1260*xi.^3+945*xi;
psi(:,11) = xi.^10-45*xi.^8+630*xi.^6-3150*xi.^4+4725*xi.^2-945;

expnsn = 2; %%% Case value to have direct and PC expansion
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

%%% Damping matrix

m_alpha = 0.01;
k_alpha = 0.001;
 
C = m_alpha *M + k_alpha*K;


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
%  M(1,1) = M(1,1) + 1* 10^(7);

 
%  K_T = M/(beta*deltaT^2) + K;
 


u_NISP = zeros(nN,count+1,n_outpce);



for ipce = 1:n_outpce
    
    
    U_MCS = zeros(nN,count+1,ns);
    
for ss = 1:ns
    

 %%%% Initial Conditions %%%%%%%%%%
uini = zeros(nN,1);

p_n = zeros(nN,1);

a_n = zeros(nN,1);

u_n = uini;

time = zeros(count+1,1);

NBcount = 1;

U_MCS(:,NBcount,ss) = u_n;
    
    
    
while NBcount<=count
    
time(NBcount+1)=time(NBcount)+deltaT;

NBcount = NBcount+1;

fun_Un_An = ((u_n + deltaT*p_n)/(deltaT^2 * beta)) + ((1-2*beta)*a_n)/(2*beta);

M_n = M * fun_Un_An;

C_n = C*gamma*deltaT*fun_Un_An - (C*p_n) - (C*(1-gamma)*a_n*deltaT);


F_T = B(:,NBcount) + M_n + C_n;


K_ss = assemble_stiffness(ne,nN,K,xi(1,ss),mu_g,sigma_g,x,l,ipce,psi,ss,n_inpce,expnsn);
    
    
 K_T = M/(beta*deltaT^2) + C*gamma/(beta*deltaT) + K_ss;
    
    
 u_new = K_T\F_T;
    

U_MCS(:,NBcount,ss) = u_new;

[u_n,p_n,a_n] = Update_acceleration(u_n,p_n,a_n,deltaT,beta,gamma,u_new);


end


U_MCS(:,:,ss) = U_MCS(:,:,ss) .* psi(ss,ipce);


end

  % % % % Numerator
u_numerator = sum(U_MCS,3)/ns;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Normalized Coefficient

u_NISP(:,:,ipce) = u_numerator/sqrt(norm.norm_squared(ipce)); 


end


U = u_NISP(:,:,1);

%%%%%%%%%%%%%%%%%%%%%%%
nt = length(time);

% Plot the results
plot_times = [100 150 250 350 500];

% Mean Solution
for i = 1:length(plot_times)
    
  figure(i)
  k = plot_times(i);
  plot(x,U(:,k),'linewidth',2);
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


% Calculation of standard deviation
U_std = zeros(nN,count+1);


% for tt = 1:nt   
    
for ipce = 2:n_outpce
    
U_std(:,:) = U_std(:,:) + u_NISP(:,:,ipce).^2;

end

% end

U_std = sqrt(U_std);

% Sd Solution
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
  fh = figure(i+10);
  set(fh, 'color', 'white'); 
  
  
end



%%%%% PCE coefficient

pc_num = 6;


for i = 1:length(plot_times)
    
  figure(i+20)
  k = plot_times(i);
  plot(x,u_NISP(:,k,pc_num),'linewidth',2);
  grid on;
  axis([min(x) max(x) -2 2]);
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);              
  titlestring = ['TIME STEP = ',num2str(k), ' TIME = ',num2str(time(k)), 'second'];
  title(titlestring ,'fontsize',14);                            
  h=gca; 
  get(h,'FontSize');
  set(h,'FontSize',14);
  fh = figure(i+20);
  set(fh, 'color', 'white'); 
  
  
end




%-------------------------------------------------------------------------%
% Movie for the travelling wave

for j = 1:nt              
  plot(x,U(:,j),'linewidth',2);
  grid on;
  axis([min(x) max(x) -5 5]);
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
cd ./results/
save('U_NISP.mat','u_NISP')
save('U_std_NISP.mat','U_std')


