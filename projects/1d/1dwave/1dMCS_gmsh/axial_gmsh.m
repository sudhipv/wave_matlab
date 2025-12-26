%%%%%% Axial Bar problem - Monte Carlo %%%%%%%

%%% EA(x,theta) = EA_o * (1+ xi *cos(4*pi*x/l)*alpha1 + (xi^2-1)
%%% *cos(8*pi*x/l)*alpha2)

%%% 1D mesh %%%

h = 0.25;


p = load('points.txt');
e = load('elements.txt');
nodes = load('nodes.txt');
% x = 0:h:1;
% 
nN = length(p);

%%%%%%%%%%%%%%

T = 1;

gamma = 0.5;

beta = 0.25;

deltaT = 0.1;

count = ceil(T/deltaT);

%%%%%%%%%%%%%%%%%
   
rho = 1;
area = 1;
l = h;
% l = p(5,1)- p(4,1)

m =    [2   1;
       1   2];
   
m_e = (1/6)*(rho * area * l)* m;

%%% Assemble mass matrix %%%

M = zeros(nN,nN);

for i = 1:length(e)
    
     
        pp = e(i,:);
        M(pp,pp) = M(pp,pp) + m;
    
    
end

EA = 1;
k_e = (EA/l)* [ 1 -1;
                -1 1];

K = zeros(nN,nN);

for i = 1:length(e)
                
K(pp,pp) = K(pp,pp) + k_e;

 end

K_T = M/(beta*deltaT^2) + K;







