



%%% 1D mesh %%%

%%% Element Size
h = 0.005;

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


%%%%%%Assemble Stiffness Matrix%%%%%%

EA = 5;
k_e = (EA/l)* [ 1 -1;
                -1 1];

K = zeros(nN,nN);

for i = 1:ne

        K(i:i+1,i:i+1) =K(i:i+1,i:i+1) + k_e;


end

%%% Force Vector

% May have to check for nyquist frequency condition
% deltaT < bt * 2*PI/(w max- frequency of force) 

wmax = 100;
f = 1000*sin(wmax.*t); 
B = zeros(nN,count+1);
% 
B(nN,:) = f;


t_nyqst = 0.5 * 2* pi/(wmax);
fprintf("delta T nyquist %2.2f\n",t_nyqst)
t_cfl = h/(sqrt(EA/rho));
fprintf("delta T CFL %4.4f\n",t_cfl)

fprintf("delta T used %4.3f\n",deltaT)


%%% Boundary Condition - Fixed at left end (node 1) adding a big penalty
%%% term

 K(1,1) = K(1,1) + 1* 10^(7);
 M(1,1) = M(1,1) + 1* 10^(7);

 
 K_T = M/(beta*deltaT^2) + K;
 
%%%% Initial Conditions %%%%%%%%%%
uini = zeros(nN,1);

p_n = zeros(nN,1);

a_n = zeros(nN,1);

u_n = uini;

U = zeros(nN,count+1);

time = zeros(count+1,1);

NBcount = 1;

U(:,NBcount) = u_n;
% f = 10;
% U(:,NBcount) = sin(2*pi*f.*time);


while NBcount<=count
    
time(NBcount+1)=time(NBcount)+deltaT;

NBcount = NBcount+1;

fun_Un_An = ((u_n + deltaT*p_n)/(deltaT^2 * beta)) + ((1-2*beta)*a_n)/(2*beta);

M_n = M * fun_Un_An;


F_T = B(:,NBcount) + M_n;
 
u_new = K_T\F_T;           

U(:,NBcount) = u_new;

[u_n,p_n,a_n] = Update_acceleration(u_n,p_n,a_n,deltaT,beta,gamma,u_new);


end


%%%%%%%%%%%%%%%%%%%%%%%
nt = length(time);

% Plot the results
plot_times = [10 50 100 150 250];


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

movie(G,1,10)
%-------------------------------------------------------------------------%



