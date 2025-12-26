

%%%% Code to plot MCS INtrusive and NISP together

%%% 7 in put 10 output
U_int = load('../wave_Intrusive/Results/sd_2/U_intrusive.mat');

%%%% 2000 samples
U_MCS = load('../wave_MCS/results/sd_2/U_MCS.mat','U_wave');

U_sd_MCS = load('../wave_MCS/results/sd_2/U_sd_MCS.mat');

U_sd_int = load('../wave_Intrusive/Results/sd_2/U_std_intrusive.mat');

time_MCS = load('../wave_MCS/results/sd_2/time_MCS.mat');

time = load('time.mat');

%%% 7 in 10 out, 1000 samples
% U_NISP = load('../wave_NISP/results/U_NISP.mat');
% 
% U_sd_NISP = load('../wave_NISP/results/U_std_NISP.mat');

x = load('x_intrusive.mat');

U_int = U_int.U;
U_MCS = U_MCS.U_wave;
U_sd_MCS = U_sd_MCS.U_sd;
U_sd_int = U_sd_int.U_std;
time = time.time;

time_MCS = time_MCS.time;
% U_NISP = U_NISP.u_NISP;
% U_sd_NISP = U_sd_NISP.U_std;
x = x.x;


nN = length(x);

% Plot the results
plot_times = [100 150 200 350 501];
plot_MCS_time = [199,299,399,699,1001];

%%% Mean 

for i = 1:length(plot_times)
    
  figure(i)
  k = plot_times(i);
  kk = plot_MCS_time(i);
  plot(x,U_MCS(1:nN,kk),'linewidth',2);
  hold on
%   plot(x,U_NISP(1:nN,k,1),'linewidth',2);
%   hold on
  plot(x,U_int(1:nN,k),'linewidth',2);
  
  grid on;
  legend;
  axis([min(x) max(x) -10 10]);
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);              
  titlestring = ['MEAN :' '   TIME STEP = ',num2str(k), '   TIME = ',num2str(time(k)), 'second'];
  title(titlestring ,'fontsize',14);                            
  h=gca; 
  get(h,'FontSize');
  set(h,'FontSize',14);
  fh = figure(i);
  set(fh, 'color', 'white'); 
  
  
end

%%%% Standard Deviation

for i = 1:length(plot_times)
    
  figure(i+10)
  k = plot_times(i);
  kk = plot_MCS_time(i);
  plot(x,U_sd_MCS(1:nN,kk),'linewidth',2);
  hold on
%   plot(x,U_sd_NISP(1:nN,k,1),'linewidth',2);
%   hold on
  plot(x,U_sd_int(1:nN,k),'linewidth',2);
  grid on;
  axis([min(x) max(x) 0 5]);
  xlabel('X axis','fontSize',14);
  ylabel('Wave Amplitude','fontSize',14);              
  titlestring = ['SD :' '   TIME STEP = ',num2str(k), '   TIME = ',num2str(time(k)), 'second'];
  title(titlestring ,'fontsize',14);                            
  h=gca; 
  get(h,'FontSize');
  set(h,'FontSize',14);
  fh = figure(i);
  set(fh, 'color', 'white'); 
  
  
end


%%%%% PCE coefficient


% % % plot_pc = [2,4,6,8];
% % % 
% % % k = 501;
% % % 
% % % for i = 1:length(plot_pc)
% % %     
% % %   figure(i+20)
% % %   pc_num = plot_pc(i);
% % %   plot(x,U_NISP(1:nN,k,pc_num),'linewidth',2);
% % %   hold on
% % %   plot(x,U_int((pc_num-1)*nN+1:(pc_num-1)*nN+nN,k),'linewidth',2);
% % %   grid on;
% % %   axis([min(x) max(x) -2 2]);
% % %   xlabel('X axis','fontSize',14);
% % %   ylabel('Wave Amplitude','fontSize',14);              
% % %   titlestring = ['PC COEFFICIENT = ',num2str(pc_num), ' TIME = ',num2str(time(k)), 'second'];
% % %   title(titlestring ,'fontsize',14);                            
% % %   h=gca; 
% % %   get(h,'FontSize');
% % %   set(h,'FontSize',14);
% % %   fh = figure(i);
% % %   set(fh, 'color', 'white'); 
% % %   
% % %   
% % % end



%%%% COV variation plot

% % U_cov = U_sd_int./abs(U_int(1:nN,:));

% % for i = 1:length(plot_times)
% %     
% %   figure(i)
% %   k = plot_times(i);
% %   plot(x,U_cov(:,k),'linewidth',2);
% %   
% %   grid on;
% %   legend;
% %   axis([min(x) max(x) -100 100]);
% %   xlabel('X axis','fontSize',14);
% %   ylabel('CoV','fontSize',14);              
% %   titlestring = ['CoV:' '   TIME STEP = ',num2str(k), '   TIME = ',num2str(time(k)), 'second'];
% %   title(titlestring ,'fontsize',14);                            
% %   h=gca; 
% %   get(h,'FontSize');
% %   set(h,'FontSize',14);
% %   fh = figure(i);
% %   set(fh, 'color', 'white'); 
% %   
% %   
% % end
% % 






