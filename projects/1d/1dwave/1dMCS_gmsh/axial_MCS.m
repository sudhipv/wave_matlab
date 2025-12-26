

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

            
 
 %%%%%%%% MCS Parameters %%%%%%%%%%%%%%%%%%%%%%%

ns = 1000;
xi = randn(1,ns);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


EA = 1;
k_e = (EA/l)* [ 1 -1;
                -1 1];

% EA = EA * (1 + xi * cos (4*pi*p(

K = zeros(nN,nN);


%%%% Initial Conditions %%%%%%%%%%
uini = zeros(nN,1);

p_n = zeros(nN,1);

a_n = zeros(nN,1);

u_n = uini;

U = zeros(nN,count+1);
U_std = zeros(nN,count+1);

time = zeros(count+1,1);

NBcount =1;

U(:,NBcount) = u_n;
f = 10;
U(:,NBcount) = sin(2*pi*f.*time);


while time<=T
    
time(NBcount+1)=time(NBcount)+deltaT;

NBcount = NBcount+1;

fun_Un_An = ((u_n + deltaT*p_n)/(deltaT^2 * beta)) + ((1-2*beta)*a_n)/(2*beta);

M_n = M * fun_Un_An;

u_mcs = zeros(nN,ns);

%%%% MCS sampling 

for ss = 1:ns

    
        for i = 1:length(e)
                
               %%% Center of each element is taken as ea varying point
                x_s = abs(p(pp(1)) - p(pp(2)))*0.5;
                %%% Multiplied by each random amplitude to form stiffness matrix
                EA_ss = 1 + xi(1,ss)*cos(4*pi*x_s/l) + (xi(1,ss)^2 - 1)*cos(8*pi*x_s/l);
                K(pp,pp) = K(pp,pp) + k_e * (EA_ss);
            


        end

           
            K_T = M/(beta*deltaT^2) + K;


            u_mcs(:,ss) = K_T\F_T;  


end

%%% Mean Solution
u_new = sum(u_mcs,2)/ns;


%%% Standard Deviation of Solution
sdum = zeros(nN,1);

 for j = 1:ns
   
     sdum = sdum + (u_mcs(:,j) - u_new).^2;
     
     
 end
 
 U_sd = sdum/(ns-1);
 U_std(:,NBcount) = sqrt(U_sd);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U(:,NBcount) = u_new;

[u_n,p_n,a_n] = Update_acceleration(u_n,p_n,a_n,deltaT,beta,gamma,u_new);


end


%%%%%%%%%%%%%%%%%%%%%%%
nt = length(time);

% Plot the results
plot_times = [2 4 10];

for i = 1:3
    
  figure(i)
  k = plot_times(i);
  plot(nodes,U(:,k),'linewidth',2);
  grid on;
%   axis([min(nodes) max(nodes) inf inf]);
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);              
  titlestring = ['TIME STEP = ',num2str(k), ' TIME = ',num2str(time(k)), 'second'];
  title(titlestring ,'fontsize',14);                            
  h=gca; 
  get(h,'FontSize') 
  set(h,'FontSize',14);
  fh = figure(i);
  set(fh, 'color', 'white'); 
  
  
end
%-------------------------------------------------------------------------%
% Movie for the travelling wave

for j = 1:T              
  plot(nodes,U(:,j),'linewidth',2);
  grid on;
%   axis([min(x) max(x) -2 2]);
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);              
  titlestring = ['TIME STEP = ',num2str(j), ' TIME = ',num2str(t(j)), 'second'];
  title(titlestring ,'fontsize',14);                            
  h=gca; 
  get(h,'FontSize') 
  set(h,'FontSize',14);
  fh = figure(5);
  set(fh, 'color', 'white'); 
  F=getframe;
            
end

movie(F,T,1)
%-------------------------------------------------------------------------%

%%% Mean Solution
for j =1:nt
    
    
     s1 = plot(nodes,U(:,j));
     hold on
     colorbar
     pause(0.05);
     
         
    if(j<count)
         delete(s1)
    end

end

%%%% Stadard Deviation of Solution


for j =1:nt
    
    
     s2 = plot(p(:,1),U_std(:,j));
     hold on
     colorbar
     pause(0.05);
     
         
    if(j<count)
         delete(s2)
    end

end








