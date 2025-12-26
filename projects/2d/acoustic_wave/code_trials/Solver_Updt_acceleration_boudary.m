

% Edited from Poissor Solver : 2003 Johan Hoffman and Anders Logg.
% Licensed under the GNU GPL Version 2.


% Acoustic Wave propagation solver

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M d^2u(t)/dt^2 + K u(t) = f

% M = 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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



 C = zeros(size(M,1),size(M,2));

% C = 0.5 .* M;


% Assemble vector
% b = AssembleVector(p, e, t, c, 4, 'Acoustic', [], 0);


%%%%%%%%%%%%%%%%%%%%%%
% New-Mark Beta Scheme

T = 1;

gamma = 0.5;

beta = 0.25;

deltaT = 0.01;

count = ceil(T/deltaT);

%%% Space coordinates %%%%

x = p(1,:);

y = p(2,:);

m = 2;
n = 1;

%%%% Initial Conditions %%%%%%%%%%
 uini = sin(m.*pi.*x).*sin(n.*pi.*y);

% uini = cos(-x-y)';

%  uini = 1*exp(-50*((x-0.5).^2+(y-0.5).^2));


%  uini = zeros(1,length(p(1,:)));


% p_n = -sqrt(2)*sin(-x-y)';

 p_n = zeros(length(p(1,:)),1);

a_n = zeros(length(p(1,:)),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_n = uini'*0.0;


U = zeros(length(p(1,:)),count+1);

time = zeros(count+1,1);

i =1;

U(:,i) = u_n;

while time<=T
    
time(i+1)=time(i)+deltaT;

% Assemble vector
b = AssembleVector(p, e, t, c, 4, 'Acoustic', [], time(i+1));

i = i+1;

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

%%%%%%%%%%%%%%%%%%%%%%%
% % % Analytic solution

nt = length(time);
% % SOL = zeros(length(x),nt); 
% % S = uini';
% % 
% % for i=1:nt
% %     SOL(:,i) = cos(sqrt(2)*time(i) - x - y);
% % end
% % 
% % %Absolute error
% % E = abs(SOL-U);
% % max_E = max(max(E))
% % 

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

pos =113 ; % 81, % 247 % 157
 
figure(10)
plot(time,U(pos,:));


% % % figure(10)
% % % plot(time,U(pos,:),time,SOL(pos,:));
% % % 
% % % 
% % % figure(11)
% % % plot(time,E(pos,:));
% % % 
% % % 
% % % E_norm = zeros(nt,1);
% % % for kk = 1:nt
% % % E_norm(kk) = norm(U(:,kk)-SOL(:,kk));
% % % end
% % % figure(12)
% % % plot(time,E_norm);
% % % 
% % % E_norm(kk-1)

% movie(Mov)

% fprintf("2nd norm of error matrix", norm(U(:,:)-SOL(:,:)))











