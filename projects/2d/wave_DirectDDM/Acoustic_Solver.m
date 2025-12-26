

% Edited from Poissor Solver : 2003 Johan Hoffman and Anders Logg.
% Licensed under the GNU GPL Version 2.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Domain Decomposition PCGM- % Acoustic Wave propagation solver %
% Two level Preconditioner
% Sudhi Sharma P V %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%GmeshExtractor
% Given the file name, 
% 1)(p_f,e_f,t_f)- extracted points, edges, and triangles for different domains. 
% 2) GLnode - Global nodes by partition
% 3) npart - Number of partitions
% 4) point_exact, triangle_exact - Used for plotting exact solutions through PDE Surf
% 5) node_gamma, node_interior - Boundary and interior nodes by partition -"RENAMED"
% 6) node_interfacewhole, node_remainingwhole - Full Interface nodes for domain...remaining
% nodes for the whole domain by partition- "NOT RENAMED"
% 7) node_corner,node_remaining -  Coarse Grid nodes "NOT RENAMED"
% 8) node_c,node_r - Coarse Grid nodes "RENAMED"
% 9) R_c,R_r - Restriction matrix from interface subdomain to subdomain corner and remaining
% nodes
% 10) Bc -  Restriction operator from corner global to corner sub domain
% 11) node_gamma_OG,node_interior_OG,node_corner_OG,node_remaining_OG -
% Original node numbers of interior, interface, corner and remaining by
% parts


clear;
delete *.mat  % Command deletes all files in folder with .mat extension

%%%%%%%% Mesh extraction from Gmsh %%%%%%%%%%%%%%%


[p_f,e_f,t_f,GLnode,npart,point_exact,triangle_exact,node_gamma,node_interior,...
    node_interfacewhole,node_interiorwhole,R,D]...
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



S_global = 0;

G_global = 0;

tol = 1*10^(-5);
tolc = 1*10^(-6);

%%%%%%%%%%%%%%%%%%%%%%
% New-Mark Beta Scheme

c = 1; %%% Wave Velocity

T = 1;

gamma = 0.5;

beta = 0.25;

deltaT = 0.01;

count = ceil(T/deltaT);

%%% Space coordinates %%%%

uini = zeros(length(node_g2),1);
p_n = zeros(length(node_g2),1);
a_n = zeros(length(node_g2),1);


x = point_exact(1,:)';

y = point_exact(2,:)';

m = 2;
n = 1;

%%% Initial Conditions %%%%%%%%%%
uini = sin(m.*pi.*x).*sin(n.*pi.*y);

% uini = cos(-x-y);

% uini = 1*exp(-50*((x-0.5).^2+(y-0.5).^2));


%   uini = zeros(1,length(p(1,:)));


% p_n = -sqrt(2)*sin(-x-y);

%  p_n(:,i) = zeros(length(p(1,:)),1);

% a_n(:,i) = zeros(length(p(1,:)),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


U_wave = zeros(length(node_g2),count+1); %%% Total Solution at all time step

U = zeros(size(nonzeros(unique(GLnode)),1),1); %%%% Solution at one timestep




for i=1:npart


p = p_f(any(p_f(:,1,i),2),2:3,i)';
e = e_f(:,:,i)';
t = t_f(any(t_f(:,1:3,i),2),:,i)';

node_ini = p_f(any(p_f(:,1,i),2),1,i);

n_gamma = unique(node_gamma(any(node_gamma(:,1,i),2),:,i));
n_interior = unique(node_interior(any(node_interior(:,1,i),2),:,i));

% Assemble Stiffness matrix
A = AssembleMatrix(p, e, t, c, 1, 'Acoustic', [], 0);

% Assemble Mass matrix
M = AssembleMatrix(p, e, t, c, 2, 'Acoustic', [], 0);

% C = zeros(size(M,1),size(M,2));
C = 0.0 * M;

% M = 0.001*M;

% % % % Assemble vector
% b = AssembleVector(p, e, t, c, 4, 'Acoustic', [], 0);


Aname = sprintf('A_%d.mat',i);
% Bname = sprintf('B_%d.mat',i);
Mname = sprintf('M_%d.mat',i);
Cname = sprintf('C_%d.mat',i);
save(Aname,'A');
% save(Bname,'b');
save(Mname,'M');
save(Cname,'C');

end


u_n = uini;

U_gamma_time = zeros(length(node_interfacewhole),count+1);

U_gamma = zeros(length(node_interfacewhole),1);

time = zeros(count+1,1);

NBcount =1;

U_wave(:,NBcount) = u_n;

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


% % % % Assemble vector
b = AssembleVector(p, e, t, c, 4, 'Acoustic', [], time(NBcount));



% Retrieving saved sub-domain matrices

Aname = sprintf('A_%d.mat',i);
% Bname = sprintf('B_%d.mat',i);
Mname = sprintf('M_%d.mat',i);
Cname = sprintf('C_%d.mat',i);

fileA = matfile(Aname);
% fileb = matfile(Bname);
fileM = matfile(Mname);
fileC = matfile(Cname);

A = fileA.A;
% b = fileb.b;
M = fileM.M;
C = fileC.C;


un_sub = u_n(node_ini);
pn_sub = p_n(node_ini);
an_sub = a_n(node_ini);

un_interface = u_n(node_interfacewhole);
pn_interface = p_n(node_interfacewhole);
an_interface = a_n(node_interfacewhole);


un_interior = un_sub(n_interior);
pn_interior = pn_sub(n_interior);
an_interior = an_sub(n_interior);

%%%%%% Assembling the Linear System to be solved for Implicit Method - for
%%%%%% Each Subdomain %%%%%% 

fun_Un_An_interior = ((un_interior + deltaT*pn_interior)/(deltaT^2 * beta)) ...
                       + ((1-2*beta)*an_interior)/(2*beta);


un_sub_g = R(1:size(n_gamma,1),:,i)*un_interface;
pn_sub_g = R(1:size(n_gamma,1),:,i)*pn_interface;
an_sub_g = R(1:size(n_gamma,1),:,i)*an_interface;

fun_Un_An_interface = ((un_sub_g + deltaT*pn_sub_g)/(deltaT^2 * beta)) + ((1-2*beta)*an_sub_g)/(2*beta);



K_transient = (M/(beta*deltaT^2))+ A; 
% + (C*gamma/(beta*deltaT))

M_n_interior = M(n_interior,n_interior) * fun_Un_An_interior + M(n_interior,n_gamma)*fun_Un_An_interface;

% M_n_schur = M(n_gamma,n_gamma) - (M(n_gamma,n_interior)* (M(n_interior,n_interior)\M(n_interior,n_gamma)));

% M_n_interface = M_n_schur * fun_Un_An_interface + (M(n_gamma,n_interior)*(M(n_interior,n_interior)\M_n_interior));


M_n_interface = M(n_gamma,n_interior) * fun_Un_An_interior + M(n_gamma,n_gamma)*fun_Un_An_interface;


% C_n = C*(p_n(node_ini) + deltaT*(1-gamma)*a_n(node_ini) - (deltaT * gamma * fun_Un_An));

b_t_interior = b(n_interior) + M_n_interior ;
b_t_interface = b(n_gamma) + M_n_interface ;

% - C_n;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% K_transient u_n+1 = b_n


%%% u_new = K_transient\b_n will be solved by Two level DDM solver 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_T = K_transient;

% fprintf('Condition number of original matrix')
% cond(M_n)

S_local = K_T(n_gamma,n_gamma) - (K_T(n_gamma,n_interior) * ((K_T(n_interior,n_interior))\ K_T(n_interior,n_gamma)));

G_local = b_t_interface - (K_T(n_gamma,n_interior) * (K_T(n_interior,n_interior) \ b_t_interior));

S_global = S_global + R(1:size(n_gamma,1),:,i)' * S_local * R(1:size(n_gamma,1),:,i);

G_global = G_global + R(1:size(n_gamma,1),:,i)' * G_local;

Ktransientname = sprintf('Kt_%d.mat',i);
Biname = sprintf('Bt_interior_%d.mat',i);
Bgname = sprintf('Bt_interface_%d.mat',i);
save(Ktransientname,'K_T');
save(Biname,'b_t_interior');
save(Bgname,'b_t_interface');


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% Starting Interface Solution %%%%%%%%%%%%%%%%%%%


% U_gamma_o = U_gamma_time(:,NBcount-1);
% r_gamma_o = G_global;

U_gamma = S_global\G_global;

%  [U_gamma,iter,nres] = DoPCGM_NeumannTwoLevel(npart,node_c,node_r,...
%     node_corner,node_gamma,node_interior,G_global,U_gamma_o,r_gamma_o,Rc,Rr,Bc,R,D,tol,tolc);


U_gamma_time(:,NBcount) = U_gamma;

node_global = nonzeros(unique(GLnode));

for nr = 1:size(node_interfacewhole,1)
 U(node_global(node_interfacewhole(nr),:)) = U_gamma(nr);
end


% Constructing the A and b Matrices after getting interface solution. Can
% be avoided by saving them and reloading at this stage.

for j=1:npart
  
    
p = p_f(any(p_f(:,1,j),2),2:3,j)';
e = e_f(:,:,j)';
t = t_f(any(t_f(:,1:3,j),2),:,j)';

n_gamma = unique(node_gamma(any(node_gamma(:,1,j),2),:,j));
n_interior = unique(node_interior(any(node_interior(:,1,j),2),:,j));

% Retrieving saved sub-domain matrices

Ktransientname = sprintf('Kt_%d.mat',j);
Biname = sprintf('Bt_interior_%d.mat',j);
Bgname = sprintf('Bt_interface_%d.mat',j);
fileA = matfile(Ktransientname);
filebi = matfile(Biname);
filebg = matfile(Bgname);


K_T = fileA.K_T;
b_t_interior = filebi.b_t_interior;
b_t_interface = filebg.b_t_interface;


U_interior  = K_T(n_interior,n_interior)\(b_t_interior - (K_T(n_interior,n_gamma) * R(1:size(n_gamma,1),:,j) * U_gamma));

for ni = 1:size(n_interior,1)
 U(GLnode(n_interior(ni),:,j)) = U_interior(ni);
end


end


%%%%%%%%%%%%%%%%% DDM solver ends %%%%%%%%%%%%%%%%%%%%%%%%%%

u_new = U;

U_wave(:,NBcount) = u_new;

[u_n,p_n,a_n] = Update_acceleration(u_n,p_n,a_n,deltaT,beta,gamma,u_new);

end


fprintf("Condition number of Schur Complement matrix")
cond(S_global)
%%%%%%%%%%%%%%%%%%%%%%%
% % % Analytic solution

nt = length(time);

SOL = zeros(length(x),nt); 
S = uini';

for i=1:nt
    SOL(:,i) = S.*cos(c*pi*(sqrt(m^2+n^2)*time(i)));
end

%Absolute error
E = abs(SOL-U_wave);
max_E = max(max(E))


for j =1:nt
    
    
     s1 = pdesurf(point_exact,triangle_exact,U_wave(:,j));
     axis([0 1 0 1 -1 1]);
     zlim([-1 1])
     hold on
     colorbar
     pause(0.05);
     
         
    if(j<count)
         delete(s1)
    end

end



%%% Draw the time vaiation at a particular point

pos = 12;
% figure(20)
% plot(time,U_wave(pos,:));


figure(10)
plot(time,U_wave(pos,:),time,SOL(pos,:));

E_norm = zeros(nt,1);
for kk = 1:nt
E_norm(kk) = norm(U_wave(:,kk)-SOL(:,kk));
end
figure(12)
plot(time,E_norm);

E_norm(kk-1)













