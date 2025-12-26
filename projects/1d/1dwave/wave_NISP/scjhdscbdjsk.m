



%%% 1D mesh %%%

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

deltaT = 0.001;

count = ceil(T/deltaT);

t = 0:deltaT:T;


 %%%%%%%% MCS Parameters %%%%%%%%%%%%%%%%%%%%%%%

ns = 5000;
xi = randn(1,ns);


sigma_g = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%% Boundary Condition
M(1,1) = M(1,1) + 1* 10^(7);

%%%%%%Assemble Stiffness Matrix%%%%%%

EA = 1;
EA_0 = 5;
k_e = (EA/l)* [ 1 -1;
                -1 1];

K = zeros(nN,nN);


%%% Force Vector

% May have to check for nyquist frequency condition
% deltaT < bt * 2*PI/(w max- frequency of force) 

wmax = 100;
f = 1000*sin(wmax.*t); 
B = zeros(nN,count+1);
% 
B(nN,:) = f;


t_nyqst = 0.5 * 2* pi/(wmax);
fprintf("delta T nyquist %4.4f\n",t_nyqst)
t_cfl = h/(sqrt(EA/rho));
fprintf("delta T CFL %4.4f\n",t_cfl)

fprintf("delta T used %4.4f\n",deltaT)

 
%%%% Initial Conditions %%%%%%%%%%
uini = zeros(nN,1);

p_n = zeros(nN,1);

a_n = zeros(nN,1);

u_n = uini;

U = zeros(nN,count+1);

time = zeros(count+1,1);

NBcount = 1;

U(:,NBcount) = u_n;
U_std = zeros(nN,count+1);


while NBcount<=count
    
time(NBcount+1)=time(NBcount)+deltaT;

NBcount = NBcount+1;

fun_Un_An = ((u_n + deltaT*p_n)/(deltaT^2 * beta)) + ((1-2*beta)*a_n)/(2*beta);

M_n = M * fun_Un_An;


F_T = B(:,NBcount) + M_n;


u_mcs = zeros(nN,ns);

%%%% MCS sampling 

% for ss = 1:ns

            for i = 1:ne

                  
                    %%% Center of each element is taken as ea varying point
                    x_s = abs(x(i) - x(i+1))*0.5;
                    %%% Multiplied by each random amplitude to form stiffness matrix
%                     EA_ss = 1 + xi(1,ss)*cos(4*pi*x_s/l) + (xi(1,ss)^2 - 1)*cos(8*pi*x_s/l);
                    
%                    EA_ss = exp(0 + sigma_g*xi(1,ss))+ EA_0;
                   EA_ss = 10;
                    
                    
                    K(i:i+1,i:i+1) = K(i:i+1,i:i+1) + k_e *(EA_ss);

            end

            
            %%% Boundary Condition - Fixed at left end (node 1) adding a big penalty
           %%% term

            K(1,1) = K(1,1) + 1* 10^(7);
            K_T = M/(beta*deltaT^2) + K;


%             u_mcs(:,ss) = K_T\F_T;  

u_new =K_T\F_T;
% end

%%% Mean Solution
% u_new = sum(u_mcs,2)/ns;


%%% Standard Deviation of Solution
% % % % sdum = zeros(nN,1);
% % % % 
% % % %  for j = 1:ns
% % % %    
% % % %      sdum = sdum + (u_mcs(:,j) - u_new).^2;
% % % %      
% % % %      
% % % %  end
% % % %  
% % % %  U_sd = sdum/(ns-1);
% % % %  U_std(:,NBcount) = sqrt(U_sd);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  

U(:,NBcount) = u_new;

[u_n,p_n,a_n] = Update_acceleration(u_n,p_n,a_n,deltaT,beta,gamma,u_new);


end


%%%%%%%%%%%%%%%%%%%%%%%
nt = length(time);

% Plot the results
plot_times = [10 50 100];


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
%-------------------------------------------------------------------------%
% Movie for the travelling wave

for j = 1:nt              
  plot(x,U(:,j),'linewidth',2);
  grid on;
  axis([min(x) max(x) -25 25]);
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

movie(G,T)
%-------------------------------------------------------------------------%



