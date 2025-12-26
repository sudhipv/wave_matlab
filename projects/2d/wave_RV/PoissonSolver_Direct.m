%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Direct Domain Decomposition - Poisson Problem %

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

clear;
delete *.mat  % Command deletes all files in folder with .mat extension

[p_f,e_f,t_f,GLnode,npart,point_exact,triangle_exact,node_gamma,node_interior,...
    node_interfacewhole,node_remainingwhole,R]...
= ExtractGmshPart_puff('mymeshpart10.msh');

% Check for interface and interior decomposing

node_g1 = sort(vertcat(node_interfacewhole,nonzeros(unique(node_remainingwhole))));
node_g2 = nonzeros(unique(GLnode));

if(any(node_g1-node_g2))
    disp("error")
end

%..............%

S_global = 0;

G_global = 0;

U = zeros(size(nonzeros(unique(GLnode)),1),1);

for i=1:npart


p = p_f(any(p_f(:,1,i),2),2:3,i)';
e = e_f(:,:,i)';
t = t_f(any(t_f(:,1:3,i),2),:,i)';

n_gamma = unique(node_gamma(any(node_gamma(:,1,i),2),:,i));
n_interior = unique(node_interior(any(node_interior(:,1,i),2),:,i));

% Assemble matrix
A = AssembleMatrix(p, e, t, 'Poisson', [], 0);

% Assemble vector
b = AssembleVector(p, e, t, 'Poisson', [], 0);

S_local = A(n_gamma,n_gamma) - (A(n_gamma,n_interior) * inv(A(n_interior,n_interior)) * A(n_interior,n_gamma));

G_local = b(n_gamma) - (A(n_gamma,n_interior) * inv(A(n_interior,n_interior)) * b(n_interior));

S_global = S_global + R(1:size(n_gamma,1),:,i)' * S_local * R(1:size(n_gamma,1),:,i);

G_global = G_global + R(1:size(n_gamma,1),:,i)' * G_local;

Aname = sprintf('A_%d.mat',i);
bname = sprintf('b_%d.mat',i);
save(Aname,'A');
save(bname,'b');

end

U_gamma = S_global\G_global;

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

Aname = sprintf('A_%d.mat',j);
bname = sprintf('b_%d.mat',j);
fileA = matfile(Aname);
fileb = matfile(bname);
A = fileA.A;
b = fileb.b;


U_interior  = A(n_interior,n_interior)\(b(n_interior) - (A(n_interior,n_gamma) * R(1:size(n_gamma,1),:,j) * U_gamma));

for ni = 1:size(n_interior,1)
 U(GLnode(n_interior(ni),:,j)) = U_interior(ni);
end


end


% Compute exact solution
u = zeros(size(U));
for i = 1:size(point_exact,2)
  x = point_exact(:,i);
  u(i) = sin(pi*x(1)) * sin(2*pi*x(2));
end

% Compute the error
error = U - u;
enorm = max(abs(error));
disp(['Maximum norm error: ' num2str(enorm)])

% Plot solution, exact solution, and error
figure(1); clf
pdesurf(point_exact,triangle_exact,U)
shading faceted
title('Computed solution')

figure(2); clf
pdesurf(point_exact,triangle_exact,u);
shading faceted
title('Exact solution')

figure(3); clf
pdesurf(point_exact,triangle_exact,error)
shading faceted
title('Error')
