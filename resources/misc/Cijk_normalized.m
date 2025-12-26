
clc
% Moments of PCE - Multiplication Tensor 
clearvars;

dim = 1;

in_ord = 10;

out_ord = 15;


% One dimensional polynomial order (maximum of input and output order)
if(in_ord > out_ord)
    ord = in_ord;
else
    ord = out_ord;
end

% Number of terms in input expansion (k index for MultiDimensional Moment )

n_pce_in = factorial(in_ord + dim)/(factorial(in_ord)*factorial(dim));

% Number of terms in output expansion (i and j index for MultiDimensional Moment)

n_pce_out = factorial(out_ord + dim)/(factorial(out_ord)*factorial(dim));


% point to UQTk directopry
uqtk_path = '/Users/sudhipv/documents/UQTk_v3.0.4/UQTk-install/bin';
cd(uqtk_path);
level = 2 * ord +1;

dim_one = 1;

[status,out] = system(['./generate_quad -d ', num2str(dim_one),' -p ', num2str(level), ' -g HG -x full']);

% load the quadrature points and weights
x =load('qdpts.dat'); % Need to scale each variable, qdpts = f(qdpts)
w =load('wghts.dat');

n_qd = size(x,1);

psi = zeros(n_qd,ord+1);
psi_N = zeros(n_qd,ord+1);

% 1-st order Hermite PC: it is always fixed to psi(1) = 1
psi(:,1) = 1;
psi_N(:,1) = 1;

if(ord > 1 || ord == 1 )
% 2nd order Hermite PC: it?s always fixed to psi(2) = x if (nord > 0)
psi(:,2) = x;
psi_N(:,2) = x;

% 3rd and more order Hermite PC?s are solved by recursive formula for i = 3 : nord+1
    for i = 3:ord+1
        psi(:,i) = x .* psi(:,i-1) - (i-2) * psi(:,i-2);
        psi_N(:,i) = psi(:,i)/sqrt(factorial(i-1));
    end
end

% 1D Moments of Hermite PC basis <psi_i psi_j psi_k>
moment_1D = zeros(ord + 1, ord + 1, ord + 1);

for k = 1 : ord + 1
    for j = 1 : ord + 1
        for i = 1 : ord + 1
             
                  mult = 0;
                  mult_N = 0;
                  
                  for m = 1:n_qd
                  mult = mult + psi(m,i)*psi(m,j)*psi(m,k)*w(m);
                  mult_N = mult_N + psi_N(m,i)*psi_N(m,j)*psi_N(m,k)*w(m);
                  end
                  
                  moment_1D(i,j,k) = mult; 
                  moment_1D_N(i,j,k) = mult_N;%   *((1/(factorial(i-1))^(1/dim))); 
             
        end
    end
end

moment_1D;
moment_1D_N;


% Creating Multi indices from UQTk

% % % point to UQTk directopry
uqtk_path = '/Users/sudhipv/documents/UQTk_v3.0.4/UQTk-install/bin';
cd(uqtk_path);

[status,out] = system(['./gen_mi -p', num2str(out_ord), ' -q', num2str(dim) , '-x TO'])

m_index_in = load('mindex.dat')


[status,out] = system(['./gen_mi -p', num2str(out_ord), ' -q', num2str(dim) , '-x TO'])

m_index_out = load('mindex.dat')


% Multi Dimensional Moments from 1D moments <Psi_i Psi_j Psi_k> = 
  
    
    for i = 1:n_pce_out
        
        for j = 1: n_pce_out
            
            for k = 1:n_pce_out
                
                product = 1;
                product_N = 1;
                norm_product = 1;
                
                    for d = 1: dim
                        
                       mi1 = m_index_in(i,d)+1;
                       mi2 = m_index_out(j,d)+1;
                       mi3 = m_index_out(k,d)+1;
                       
                       %% Moment Normalized
                       moment = product * moment_1D(mi1,mi2,mi3);
                       product = moment;
                       
                       %% Each Basis normalized
                       moment_N = product_N * moment_1D_N(mi1,mi2,mi3);
                       product_N = moment_N;
                       
%                        if(d < 2 || d ==2)
                           norm_mik = m_index_out(k,d);
                           norm_moment = norm_product * factorial(norm_mik);
                           norm_product = norm_moment;
%                        end

                    end
                    
                  norm_squared(k) = norm_moment;
                  
                  if(moment < 1*10^(-6))
                      moment = 0;
                  end
                  moment_NonNormalized(i,j,k) = moment;
                  moment_ND(i,j,k) = moment/norm_squared(k);
                  
                  if(product_N < 1*10^(-6))
                      product_N = 0;
                  end
                  moment_PsiNorm(i,j,k) = product_N;
                  
            end
        end
    end
   
   moment_NonNormalized;
   moment_ND ;
   norm_squared;
   
   
   %%% Change the folders and variable name accordingly for NORMLAIZED and
   %%% NON-NORMALIZED multiplication TENSOR
   
 cd /Users/sudhipv/documents/wave_matlab/misc/cijk_ord_dim/;

% cd ../../../ssfem_matlab/Intrusive/;


str1 = strcat('0',int2str(ord));

if dim <= 10
    str2 = strcat('000',int2str(dim));
else
    str2 = strcat('00',int2str(dim));   
end

str3 = strcat('Cijk',str1,str2);
str4 = strcat('norm_squared',str1,str2);
str5 = strcat('mindex',str1,str2);
        
 save(str3,'moment_PsiNorm');
% save(str3,'moment_NonNormalized');

save(str4,'norm_squared');

save(str5,'m_index_out');




% % % save('/Users/sudhipv/documents/UQTk_Matlab/norm_k.mat','norm_squared');


% Quadrature weights and points

% % point to UQTk directopry
% uqtk_path = '/Users/sudhipv/documents/UQTk_v3.0.4/UQTk-build/bin';
% cd(uqtk_path);


% Non-dimentionalized, 1-dimensional Hermite polynomials 
% 1-dimensional Hermite polynomials at quad-points