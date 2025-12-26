

% Edited from Poissor Solver : 2003 Johan Hoffman and Anders Logg.
% Licensed under the GNU GPL Version 2.


% Acoustic Wave propagation solver

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M d^2u(t)/dt^2 + K u(t) = f

% M = 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

% Load the mesh
% square
clear
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
% M = 0.001*M;


% C = zeros(size(M,1),size(M,2));

%%% For Gaussian Pulse
C = 0.5 .* M + 0.01*A;

% C = 0.5 .* M + 0.001*A;

%%%%%%%%%%%%%%%%%%%%%%
% New-Mark Beta Scheme

T = 1;

gamma = 0.5;

beta = 0.25;

deltaT = 0.005;

count = ceil(T/deltaT);

%%% Space coordinates %%%%

x = p(1,:)';

y = p(2,:)';

m = 2;
n = 1;

%%%% Initial Conditions %%%%%%%%%%
% Initial Condition for verified
% solution
% uini = sin(m.*pi.*x).*sin(n.*pi.*y); 


%%% Initail Condition for Forcing term
uini = zeros(length(p(1,:)),1);

p_n = zeros(length(p(1,:)),1);

a_n = zeros(length(p(1,:)),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_n = uini;


U = zeros(length(p(1,:)),count+1);

time = zeros(count+1,1);

i =1;

U(:,i) = u_n;

while time<=T
    
time(i+1)=time(i)+deltaT;

i = i+1;

% Assemble vector
b = AssembleVector(p, e, t, c, 4, 'Acoustic', [], i);


fun_Un_An = ((u_n + deltaT*p_n)/(deltaT^2 * beta)) + ((1-2*beta)*a_n)/(2*beta);

K_transient = (M/(beta*deltaT^2)) + (C*gamma/(beta*deltaT))+ A;

M_n = M * fun_Un_An;

C_n = C*(p_n + deltaT*(1-gamma)*a_n - (deltaT * gamma * fun_Un_An));


b_n = b + M_n - C_n;

u_new = K_transient\b_n;


U(:,i) = u_new;

[u_n,p_n,a_n] = Update_acceleration(u_n,p_n,a_n,deltaT,beta,gamma,u_new);

end

cnd_num = cond(K_transient);

toc

nt = length(time);


for j =1:nt
    
    
     s1 = pdesurf(p,t,U(:,j));
     axis([0 1 0 1 -.5 .5]);
%      zlim([-1 1])
     view(-45,30)
     hold on
     colorbar
     pause(0.05);
     
         
    if(j<count)
         delete(s1)
    end


end



% % pos =81 ; % 81, % 247 % 157
% % 
% % figure(10)
% % plot(time,U(pos,:));

%%% Draw the time vaiation at a particular point

pos =246 ; % 81, % 247 % 157


%%%% To plot the error norm for verification

figure(11)
plot(time,U(pos,:));


% fprintf("2nd norm of error matrix", norm(U(:,:)-SOL(:,:)))




%%%%%%% Verification with Poisson %%%%%%%%%%%%%

% % % % u = zeros(length(U(:,15)),1);
% % % % for i = 1:size(p,2)
% % % %   x = p(:,i);
% % % %   u(i) = sin(pi*x(1)) * sin(2*pi*x(2));
% % % % end
% % % % 
% % % % % Compute the error
% % % % error = U(:,15) - u;
% % % % enorm = max(abs(error));
% % % % disp(['Maximum norm error: ' num2str(enorm)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








