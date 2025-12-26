

%%%%% Code to compare 2D acoustic Wave Results


pathintrusive = '/Users/sudhipv/documents/wave_matlab/wave_RV/results/7_10_sd_3_nodamp/';
pathMCS = '/Users/sudhipv/documents/wave_matlab/wave_MCS_RV/results/3000_sd_3_nodamp/';

fullFileIntr = fullfile(pathintrusive, 'U_intrusive.mat');
fullFileMCS = fullfile(pathMCS, 'U_MCS.mat');

fullFileIntr_sd = fullfile(pathintrusive, 'U_intr_std.mat');
fullFileMCS_sd = fullfile(pathMCS, 'U_MCS_sd.mat');

fullFiletime = fullfile(pathintrusive, 'time.mat');
fullFilep = fullfile(pathintrusive, 'p.mat');
fullFilet = fullfile(pathintrusive, 't.mat');

U_int = load(fullFileIntr);
U_sd_int = load(fullFileIntr_sd);
U_MCS = load(fullFileMCS);
U_sd_MCS = load(fullFileMCS_sd);

time = load(fullFiletime);
p = load(fullFilep);
t = load(fullFilet);


%%% 7 in put 10 output
% % U_int = load('/Users/sudhipv/documents/wave_matlab/wave_RV/results/7_10_sd_1_nodamp/U_intrusive.mat');
% % 
% % %%%% 2000 samples
% % U_MCS = load('/Users/sudhipv/documents/wave_matlab/wave_MCS_RV/results/1000_sd_1_nodamp/U_MCS');
% % 
% % U_sd_MCS = load('/Users/sudhipv/documents/wave_matlab/wave_MCS_RV/results/1000_sd_1_nodamp/U_MCS_sd.mat');
% % 
% % U_sd_int = load('/Users/sudhipv/documents/wave_matlab/wave_RV/results/7_10_sd_1_nodamp/U_intr_std.mat');
% % 
% % time = load('/Users/sudhipv/documents/wave_matlab/wave_RV/results/7_10_sd_1_nodamp/time.mat');
% % 
% % p = load('/Users/sudhipv/documents/wave_matlab/wave_RV/results/7_10_sd_1_nodamp/p.mat');
% % 
% % t = load('/Users/sudhipv/documents/wave_matlab/wave_RV/results/7_10_sd_1_nodamp/t.mat');


time = time.time;
p = p.p;
t= t.t;
U_int = U_int.U_wave;
U_MCS = U_MCS.U_wave;
U_sd_MCS = U_sd_MCS.U_sd;
U_sd_int = U_sd_int.U_std;



%%%% Mean at a point

pos = 400;

figure(24)
plot(time,U_MCS(pos,:),time,U_int(pos,:),'linewidth',2)
grid on


figure(25)
plot(time,U_sd_MCS(pos,:),time,U_sd_int(pos,:),'linewidth',2)
grid on



