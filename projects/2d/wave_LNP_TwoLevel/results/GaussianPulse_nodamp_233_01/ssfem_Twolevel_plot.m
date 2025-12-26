

%%%% Results Plotting for Two Level SSFEM 

load ('./s_spectral.mat')
load('point_exact.mat')
load('s_spectral.mat')
load('triangle_exact.mat')
load('U_std.mat')
load('U_wave.mat')



% Plot solution

nt = size(U_wave,2);
np = length(point_exact);

%%% Mean Solution
for j =1:nt
    
     s1 = pdesurf(point_exact,triangle_exact,U_wave(1:np,j));
     axis([0 1 0 1 -1 1]);
%      zlim([-1 1])
     hold on
     colorbar
     pause(0.05);
     
         
    if(j<nt)
         delete(s1)
    end

end


%%%%%%% Wave plot of PC coefficients %%%%%%%%%%

t_pos = 35;
indexpc = 2;

for j =1:nt
    
    
     s2 = pdesurf(point_exact,triangle_exact,U_wave(np*(indexpc-1)+1:np*(indexpc),j));
     axis([0 1 0 1 -1 1]);
     zlim([-1 1])
     hold on
     colorbar
     pause(0.05);
     
         
    if(j<nt)
         delete(s2)
    end

end



%%%%%%%%%%% PC coefficients at partcular time step %%%%%%%%


for j = 1:s_spectral.n_outpce
figure(j); clf
pdesurf(point_exact,triangle_exact,U_wave((j-1)*np+1:j*np,t_pos))
shading faceted
if(j==1)
title('Mean vector of solution expansion')
else
    titlename = sprintf('%d th PC coefficient of solution using %d RV with %drd order expansion',j,s_spectral.dim,s_spectral.ord_out);
    title(titlename)
end

end



% Calculation of standard deviation

U_std = zeros(size(point_exact,2),count+1);

for j = 2:n_outpce
U_std(:,:) = U_std(:,:) + U_wave((j-1)*node_total+1:(j-1)*node_total+node_total,:).^2;
end

U_std = sqrt(U_std);


for j =1:nt
    
    
     s3 = pdesurf(point_exact,triangle_exact,U_std(:,j));
     axis([0 1 0 1 -1 1]);
     zlim([-1 1])
     hold on
     colorbar
     pause(0.05);
     
         
    if(j<nt)
         delete(s3)
    end

end




%%%% Mean at a point

pos = 81;

figure(34)
plot(time,U_wave(pos,:))
grid on


figure(35)
plot(time,U_std(pos,:))
grid on


