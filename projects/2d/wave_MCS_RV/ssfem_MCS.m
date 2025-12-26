

%%% NON INTRUSIVE SSFEM
clearvars
%rng('default')
% Load the mesh
% square
%  square_refined

% Generated from 'ExtractGmsh.m' 
p = load('points.txt')';
e =load('edges.txt')';
t =load('triangles.txt')';
%  square_refined

c = 1;

% Assemble Stiffness matrix
A = AssembleMatrix(p, e, t, 1, 'Acoustic', [], 0);


% Assemble Mass matrix
M = AssembleMatrix(p, e, t, 2, 'Acoustic', [], 0);

%%%%%%%%%%%%%%%%%%%%%%
% New-Mark Beta Scheme

T = 1;

gamma = 0.5;

beta = 0.25;

deltaT = 0.01;

count = ceil(T/deltaT);

%%% Space coordinates %%%%

x = p(1,:)';

y = p(2,:)';

m = 2;
n = 1;


%%%%%%%% MCS Parameters %%%%%%%%%%%%%%%%%%%%%%%

ns = 3000;
xi = randn(1,ns);


psi(:,1) = ones(1,ns);
psi(:,2) = xi;
psi(:,3) = xi.^2-1;
psi(:,4) = xi.^3-3.*xi;

mu_g = 0;
sigma_g = 0.3; %%% Change folder name at end to save file
mu_l = exp(mu_g + sigma_g^2/2);
kappa(1) = mu_l;
kappa(2) = mu_l*sigma_g;
kappa(3) = mu_l*sigma_g^2/2;   
kappa(4) = mu_l*sigma_g^3/6;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% MCS Sampling %%%%%%

u_MCS = zeros(length(p(1,:)),count+1,ns);
b_time = zeros(length(p(1,:)),count+1);


for ss = 1:ns

%%%% Initial Conditions %%%%%%%%%%
uini = sin(m.*pi.*x).*sin(n.*pi.*y);

p_n = zeros(length(p(1,:)),1);

a_n = zeros(length(p(1,:)),1);

u_n = uini;

U = zeros(length(p(1,:)),count+1);

time = zeros(count+1,1);

i =1;

U(:,i) = u_n;


while time<=T
    
time(i+1)=time(i)+deltaT;

i = i+1;

fun_Un_An = ((u_n + deltaT*p_n)/(deltaT^2 * beta)) + ((1-2*beta)*a_n)/(2*beta);

M_n = M * fun_Un_An;

if(ss==1)
% % % % Assemble vector
b = AssembleVector(p, e, t, 4, 'Acoustic', [], time(i));
b_time(:,i) = b;
else
    b = b_time(:,i);
end
% C_n = C*(p_n + deltaT*(1-gamma)*a_n - (deltaT * gamma * fun_Un_An));

b_n = b + M_n;

%    A_MCS = (A.*kappa(1).*psi(i,1)+A.*kappa(2).*psi(i,2)+A.*kappa(3).*psi(i,3)+A.*kappa(4).*psi(i,4));
 
 A_MCS = A.*(exp(mu_g+sigma_g*xi(1,ss)));

%    A_MCS = A.*((mu_g+sigma_g*xi(1,ss)));
   
 K_transient = (M/(beta*deltaT^2)) + A_MCS;

 u_new = K_transient\b_n;     
      
U(:,i) = u_new;

[u_n,p_n,a_n] = Update_acceleration(u_n,p_n,a_n,deltaT,beta,gamma,u_new);

end

%%% 1 realization of MCS sampling 

u_MCS(:,:,ss) = U;


end


%%%%%%%%%%%%%%%%%%%%%%%
nt = length(time);


% % % Mean Solution
U_wave = sum(u_MCS,3)/ns;


% % % Standard Deviation of Solution
sdum = zeros(length(p(1,:)),count+1);

 for j = 1:ns
   
     sdum = sdum + (u_MCS(:,:,j) - U_wave).^2;
     
     
 end
 
 U_sd = sdum/(ns-1);
 
 U_sd = sqrt(U_sd);



%%% Mean Solution
for j =1:nt
    
    
     s1 = pdesurf(p,t,U_wave(:,j));
     axis([0 1 0 1 -1 1]);
     zlim([-1 1])
     hold on
     colorbar
     pause(0.05);
     
         
    if(j<count)
         delete(s1)
    end

end

%%%% Stadard Deviation of Solution


for j =1:nt
    
    
     s2 = pdesurf(p,t,U_sd(:,j));
     axis([0 1 0 1 -1 1]);
     zlim([-1 1])
     hold on
     colorbar
     pause(0.05);
     
         
    if(j<count)
         delete(s2)
    end

end



%%%% Mean at a point

pos = 92;

figure(24)
plot(time,U_wave(200,:))
grid on


figure(25)
plot(time,U_sd(200,:))
grid on

cd ./results/3000_sd_3_nodamp
save('U_MCS.mat','U_wave')
save('U_MCS_sd.mat','U_sd')








