

%%%% Code to check MCS convergence


U_MCS_250 = load('../../wave_MCS/results/convergence/250/U_MCS.mat','U_wave');
U_MCS_500 = load('../../wave_MCS/results/convergence/500/U_MCS.mat','U_wave');
U_MCS_1000 = load('../../wave_MCS/results/convergence/1000/U_MCS.mat','U_wave');
U_MCS_2000 = load('../../wave_MCS/results/convergence/2000/U_MCS.mat','U_wave');

time = load('../time.mat');

x = load('../x_intrusive.mat');

U_MCS_250 = U_MCS_250.U_wave;
U_MCS_500 = U_MCS_500.U_wave;
U_MCS_1000 = U_MCS_1000.U_wave;
U_MCS_2000 = U_MCS_2000.U_wave;

time = time.time;


x = x.x;


nN = length(x);

% Plot the results
plot_times = [100 150 200 350 501];

%%% Mean 

for i = 1:length(plot_times)
    
  figure(i)
  k = plot_times(i);
  plot(x,U_MCS_250(1:nN,k),'linewidth',2);
  hold on
  plot(x,U_MCS_500(1:nN,k),'linewidth',2);
  hold on
  plot(x,U_MCS_1000(1:nN,k),'linewidth',2);
  hold on
  plot(x,U_MCS_2000(1:nN,k),'linewidth',2);
 
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

% % % for i = 1:length(plot_times)
% % %     
% % %   figure(i+10)
% % %   k = plot_times(i);
% % %   kk = plot_MCS_time(i);
% % %   plot(x,U_sd_MCS(1:nN,kk),'linewidth',2);
% % %   hold on
% % % %   plot(x,U_sd_NISP(1:nN,k,1),'linewidth',2);
% % % %   hold on
% % %   plot(x,U_sd_int(1:nN,k),'linewidth',2);
% % %   grid on;
% % %   axis([min(x) max(x) 0 5]);
% % %   xlabel('X axis','fontSize',14);
% % %   ylabel('Wave Amplitude','fontSize',14);              
% % %   titlestring = ['SD :' '   TIME STEP = ',num2str(k), '   TIME = ',num2str(time(k)), 'second'];
% % %   title(titlestring ,'fontsize',14);                            
% % %   h=gca; 
% % %   get(h,'FontSize');
% % %   set(h,'FontSize',14);
% % %   fh = figure(i);
% % %   set(fh, 'color', 'white'); 
% % %   
% % %   
% % % end





