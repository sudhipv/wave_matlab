% Edited from Poissor Solver : 2003 Johan Hoffman and Anders Logg.
% Licensed under the GNU GPL Version 2.


% Acoustic Wave propagation solver

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

% C = 0.0 .* M;

C = zeros(size(M,1),size(M,2));


% Assemble vector
b = AssembleVector(p, e, t, c, 4, 'Acoustic', [], 0);


%%%%%%%%%%%%%%%%%%%%%%
% New-Mark Beta Scheme

T = 3;

gamma = 0.5;

beta = 0.25;

% CFL_mult = (1/sqrt(2));

% deltaT = (CFL_mult/c) * 0.05;

deltaT = 0.001;

count = ceil(T/deltaT);

%%% Space coordinates %%%%

x = p(1,:);

y = p(2,:);

m = 2;
n = 1;

%%%% Initial Conditions %%%%%%%%%%
uini = sin(m.*pi.*x).*sin(n.*pi.*y);

% uini = cos(-x-y);

%  uini = 1*exp(-50*((x-0.5).^2+(y-0.5).^2));


%  uini = zeros(1,length(p(1,:)));


% p_n = -sqrt(2)*sin(-x-y);

 p_n = zeros(length(p(1,:)),1);

a_n = zeros(length(p(1,:)),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_n = uini';


U = zeros(length(p(1,:)),count+1);


time = zeros(count+1,1);

i =1;

U(:,i) = u_n;

while time<=T
    
time(i+1)=time(i)+deltaT;

i = i+1;

A_n = A*(u_n + deltaT.*p_n + (deltaT^2/2)*(1-2*beta).*a_n);

C_n = C*(p_n + deltaT*(1-gamma).*a_n);

b_n = b - (A_n + C_n);

M_n = M + beta*deltaT^2*A + C*gamma;

a_new = M_n\b_n;

[u_n,p_n,a_n] = Update(u_n,p_n,a_n,deltaT,beta,gamma,a_new);

U(:,i) = u_n;


end

cnd_num = cond(M_n);

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
    
    
% % %      s1 = pdesurf(p,t,U(:,j));
% % %      axis([0 1 0 1 -1 1]);
% % %      zlim([-1 1])
% % %      hold on
% % %      colorbar
% % %      pause(0.05);
% % %      
% % %          
% % %     if(j<count)
% % %          delete(s1)
% % %     end




     
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













