% Edited from Poissor Solver : 2003 Johan Hoffman and Anders Logg.
% Licensed under the GNU GPL Version 2.


% Acoustic Wave propagation solver -JOURNAL NEWMARK BETA APPROACH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M d^2u(t)/dt^2 + K u(t) = f

% M = 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Load the mesh
% square
clear all
% square_refined

p = load('points.txt')';
e =load('edges.txt')';
t =load('triangles.txt')';

c = 1;

% Assemble Stiffness matrix
A = AssembleMatrix(p, e, t, c, 1, 'Acoustic', [], 0);


% Assemble Mass matrix
M = AssembleMatrix(p, e, t, c, 2, 'Acoustic', [], 0);

% % % Assemble Damping Matrix
% % 
% % C = AssembleMatrix(p, e, t, c, 3, 'Acoustic', [], 0);



 C = zeros(size(M,1),size(M,2));

% C = 0.0*M;


% Assemble vector
b = AssembleVector(p, e, t, c, 4, 'Acoustic', [], 0);


%%%%%%%%%%%%%%%%%%%%%%
% New-Mark Beta Scheme

T = 1;

gamma = 3/2;

beta = 1;


%%%%

% CFL_mult = (1/sqrt(2));

% deltaT = (CFL_mult/c) * 0.05;

deltaT = 0.01;

count = ceil(T/deltaT);

%%% Space coordinates %%%%

x = p(1,:);

y = p(2,:);

m = 2;
n = 1;

%%%% Initial Conditions %%%%%%%%%%
uini = sin(m.*pi.*x).*sin(n.*pi.*y);

% u_ini = 0*uini;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


U = zeros(length(p(1,:)),count+3);


time = zeros(count+3,1);

i =1;

U(:,i) = zeros;

i = i+1;

U(:,2) = uini;

while time<=T
    
time(i+1)=time(i)+deltaT;

%%%% New mark beta - 2 step %%%%%%%%%%%%%%%

A_n = A*( ((1/2)- 2*beta + gamma).*U(:,i) + ((1/2)+ beta - gamma).*U(:,i-1));

C_n = (C/deltaT)*((1-2*gamma)*U(:,i) + (gamma-1)*U(:,i-1));

M_n = (M/deltaT^2)*(U(:,i-1) - 2*U(:,i));

if(i == 2)
    b_n = beta*b + ((1/2)- 2*beta + gamma)*b ;
else
    b_n = beta*b + ((1/2)- 2*beta + gamma)*b + ((1/2)+ beta - gamma)*b;
end


K_T = (M/deltaT^2) + (C*gamma/deltaT) + beta*A;

F_T = b_n - M_n - C_n - A_n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%Houbolt %%%%%%%%%%%%%%%%%%%%%%

% % % % % M_n = (M/deltaT^2)*(-5*U(:,i) + 4*U(:,i-1)- U(:,i-2));
% % % % % 
% % % % % C_n =  (C/6*deltaT)*(-18*U(:,i)+9*U(:,i-1) - 2*U(:,i-2));
% % % % % 
% % % % % b_n = b-M_n-C_n;
% % % % % 
% % % % % K_T = (M*2/deltaT^2) + (C*11/6*deltaT) + A;
% % % % % 
% % % % % 
% % % % % F_T = b_n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_new = K_T\F_T;

i = i+1;

% [u_n,p_n,a_n] = Update_acceleration(u_n,p_n,a_n,deltaT,beta,gamma,u_new);

U(:,i) = u_new;


 end

cnd_num = cond(K_T);

%%%%%%%%%%%%%%%%%%%%%%%
% % % Analytic solution

nt = length(time);
SOL = zeros(length(x),nt); 
S = uini';

for i=1:nt
    SOL(:,i) = S.*cos(c*pi*(sqrt(m^2+n^2)*time(i)));
end

%Absolute error
E = abs(SOL-U);
max_E = max(max(E))


for j =1:nt
    
    
     s1 = pdesurf(p,t,U(:,j));
     axis([0 1 0 1 -1 1]);
     zlim([-1 1])
     hold on
     colorbar
     pause(0.05);
     
         
    if(j<count)
         delete(s1)
    end




     
% % % % %      subplot(2,1,1)  
% % % % %      s2 = pdesurf(p,t,U(:,j)); 
% % % % %      pause(0.01);
% %   

% % % %      subplot(2,1,2)  
% % % %      s3 = pdesurf(p,t,E(:,j));
% % % %      title(sprintf('Absolute error at t = %1.2f',time(j)),'Fontsize',11);
% % % %      xlabel('x','Fontsize',11); ylabel('y','Fontsize',11);
% % % %      zlabel(sprintf('E(x,y,t = %1.2f)',time(j)),'Fontsize',11);
% % % %      pause(0.01);



end



%%% Draw the time vaiation at a particular point

pos = 81;

figure(10)
plot(time,U(pos,:),time,SOL(pos,:));

E_norm = zeros(nt,1);
for kk = 1:nt
E_norm(kk) = norm(U(:,kk)-SOL(:,kk));
end
figure(12)
plot(time,E_norm);

E_norm(kk-1)

% movie(Mov)













