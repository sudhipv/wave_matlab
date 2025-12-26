function integral = Acoustic_ss(u, v, w, du, dv, dw, dx, ds, x, d, t, eq, s_spectral,ipce,file)

% Variational formulation of Poisson's equation.
% Copyright (C) 2003 Johan Hoffman and Anders Logg.
% Licensed under the GNU GPL Version 2.

if eq == 1
  kappa = GenerateCoeff(s_spectral,x, d, t,file,ipce);  
  integral = kappa * du'*dv*dx + g(x,d,t)*u*v*ds;
  %integral = kappa(n_ipce) * du'*dv*dx;
  
elseif eq == 2
  integral =  u*v*dx + g(x,d,t)*u*v*ds;
else
  integral = f(x,d,t)*v*dx + (g(x,d,t)*gd(x,d,t) - gn(x,d,t))*v*ds;
end

% Function to generate the coefficients of Spectral Expansion for Random
% Diffusivity

function kappa = GenerateCoeff(s_spectral, x, d, t,file,ipce)


switch (s_spectral.num_spectral) 
    
    % Gaussian Random variable
    case 1 
     
        if(ipce == 1)
          kappa = s_spectral.mu_g;
        else
          kappa = s_spectral.sigma_g;
        end
        
    % Log Normal Random Variable
    case 2 
     
     mu_g = s_spectral.mu_g;
     sigma_g = s_spectral.sigma_g;
     dim = s_spectral.dim;
     
     L_0 = exp(mu_g + 0.5 * sigma_g^2);

%       L_0 = 1;

% Direct value to check with ajits code
%     L_0 = 0.04;


m_index_in = s_spectral.mindex.m_index_out;
norm = s_spectral.norm;
    

% for npcin = 1:s_spectral.n_inpce
    
    numerator = 1;
    for ii = 1:dim
        numerator = numerator * sigma_g^(m_index_in(ipce,ii));
    end
        kappa = L_0*(numerator/sqrt(norm.norm_squared(ipce)));
        
        
     % Gaussian process  
     
    case 3 
      
     % Log Normal Random process
    case 4    
        
        
        [kappa] = GetKLEterms_Normalized(x,d,t,file,s_spectral,ipce);
        
% % %     cfolder = pwd;
% % %         % % % point to UQTk directopry
% % %     uqtk_path = '/Users/sudhipv/documents/UQTk_v3.0.4/UQTk-install/bin';
% % %     cd(uqtk_path);
% % %         
% % %     [status,out] = system(['./gen_mi -p', num2str(s_spectral.ord_out), ' -q', num2str(s_spectral.dim) , '-x TO'])

        
end


%--- Conductivity (penalty factor) ---
function y = g(x, d, t)
y = 1e7;
%y = 0.001;

%--- Dirichlet boundary condition ----
function y = gd(x, d, t)
y = 0;

%--- Neumann boundary condition ---
function y = gn(x, d, t)
y = 0;

%--- Right-hand side, source term ---
function y = f(x, d, t)
%     y = 0;
    
 if(t ==2)
   y = 3000*exp(-1000*((x(1)-0.5).^2+(x(2)-0.5).^2));
 else
    y = 0;
 end
%  y = 2*pi^2*sin(pi*x(1))*sin(pi*x(2));
%  y = 5*pi^2*sin(pi*x(1))*sin(2*pi*x(2));








