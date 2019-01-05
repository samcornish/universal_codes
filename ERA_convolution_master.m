% This is a semi-universal script for calculating historical reconstructions of
% variables for which impulse responses have been calculated, where the
% forcing is the first three PCs of sea-level pressure in a given model. In
% my case, this is HiGEM or MO GC2.


% The script allows the user to specify the target timeseries and the
% impulse responses, as well as naming conventions. Pre-determined,
% however, is the forcing that it will be convolved with. This script can
% therefore only be used for variables which have model-derived impulses to the
% first three PCs of SLP. 

%% Loading variables and intial processing - user should edit in this section only 
clear
close all

% addpath to the location of the variables
% ------------------------------------------------- % vv
% ERA forcing
addpath /home/ocean2/samc/HiGEMArcticCRFs/fw_content_atm_circ/analysis/variables/
% Target time-series and impulse responses
addpath /home/ocean_personal_data/samc/HiGEMArcticCRFs/sea_ice_atm_circ/variables/
% addpath to the location of the function eshade2.m 
addpath /home/ocean2/samc/HiGEMArcticCRFs/fw_content_atm_circ/analysis/functions/
% ------------------------------------------------- % ^^

% load variables
% ------------------------------------------------- % vv
% load ERA forcing: regressions on eofs of choice.
load ERA_20Ci_pressure_regressions_on_eofs.mat
P_ERA = [detrend(norm_regression_ERA_e1,'linear'), detrend(norm_regression_ERA_e2,'linear'), detrend(norm_regression_ERA_e3,'linear')]; % change names and detrending convention as appropriate
% load impulses
load ice_extent_BK_HiGEM_G.mat 
load ice_extent_BK_HiGEM_step_resp.mat grand_SE grand_SR % for sake of plotting
load ice_extent_BK_HiGEM_recon.mat  % for sake of plotting
load /home/ocean_personal_data/samc/GS_experiments/variables/ice_extent.mat ice_extent_BK_ds
% load comparison dataset - comment out/change as appropriate
% load Rabe_data.mat 
% ref_data = Rabe_data; % re-name as appropriate
% ------------------------------------------------- % ^^

%process reference data
% s_R = size(ref_data);
% m_R = mean(ref_data(:,2));
% ref_data(:,2) = (ref_data(:,2) - m_R).*10^13;

% Specify time duration of ERA data?
t = 1416;   % 118 years

% flip any step responses when plotting? Insert -1 if so.
f1 = 1; f2 = 1; f3 = 1;

target_ds = detrend(ice_extent_BK_ds,'linear');

% Determining the strings used to save files
% saving variables
% insert the directory name for saving the variables
% ------------------------------------------------- % vv
var_dir = '/home/ocean_personal_data/samc/HiGEMArcticCRFs/sea_ice_atm_circ/variables/';
% directory for figures:
fig_dir = '/home/ocean_personal_data/samc/HiGEMArcticCRFs/sea_ice_atm_circ/figures/';
% ------------------------------------------------- % ^^

% filename identifier - this should be a word or several that uniquely
% identifies the saved variables to the input variables. Usually this will
% be the name of the target timeseries and the model or record it is from
% ------------------------------------------------- % vv
save_id = 'ice_extent_BK_Higem';
title_id = 'HiGEM sea-ice in BK seas';    % for the sake of plotting (spaces are fine!)
target_id = 'Sea-ice extent';          % for y axis label (spaces are fine!)
target_units = 'm^2';
tseries_units = 'months';
% ------------------------------------------------- % ^^

std_G = G_grand_std;   % G_grand_std
mean_G = G_grand_mean;    % G_grand_mean
impulse_err = grand_impulse_err;
SR = grand_SR;
SE = grand_SE;
g = length(mean_G);

CR_sum = sum(grand_CR,2);
CR_sum_errs = (grand_CR_errs(:,1).^2 + grand_CR_errs(:,2).^2 + grand_CR_errs(:,3).^2).^0.5; %We sum the contributions of each PC to the error on the control run reconstruction

t_plus = g;
S = zeros(t+t_plus,3);
E = zeros(t+t_plus,3);

% initialising
pred_target_ERA = zeros(t+t_plus,3);
pred_err_ERA_res = zeros(t+t_plus,3);
pred_err_ERA_spread = zeros(t+t_plus,3);
S_ERA = zeros(t+t_plus,3);
E_ERA = zeros(t+t_plus,3);


%% ERA REANALYSIS

for k = 1:3
        

            clearvars X Y M sd i G Y step_res epsilon var G_var S2

            X = detrend(P_ERA(:,k),'constant'); % X is de-trended NCEP kth PC with dimensions 456 x 1
            sd = std(X);
            X = X/sd;    % The column vector X (PC)

             M = zeros(g,t+t_plus);     %  initialising a matrix of the forcing, 
            for i = 1:g
                M(i,i:t-1+i) = X;  % second dimension will become larger that t + tplus
            end
            M = M(:,1:t+t_plus);  % trim second dimension
            M = M(1:360,:);   % Now trimming to minimum cutoff lag - subject to change
        
            pred_target_ERA(:,k) = mean_G(:,k)'*M;   
            pred_err_ERA_res(:,k) = (diag(M'*G_var_save(:,:,k)*M)).^0.5;
            pred_err_ERA_spread(:,k) = std_G(:,k)'*M;
        

    %reshape
    S_ERA(:,k) = pred_target_ERA(:,k);
    E_ERA(:,k) = (pred_err_ERA_res(:,k).^2 + pred_err_ERA_spread(:,k).^2).^0.5;
   
    
end

%% Plots and saving

figure;
hold on
plot(1:length(SR(:,3)),f3.*SR(:,3),'y','LineWidth',2)
plot(1:length(SR(:,2)),f2.*SR(:,2),'r','LineWidth',2)
plot(1:length(SR(:,1)),f1.*SR(:,1),'b','LineWidth',2)
legend('PC3','PC2','PC1')
errorbar(1:length(SR(:,3)), f3.*SR(:,3), SE(:,3), 'y');
errorbar(1:length(SR(:,2)), f2.*SR(:,2), SE(:,2), 'r');
errorbar(1:length(SR(:,1)), f1.*SR(:,1), SE(:,1), 'b');

set(gca,'FontSize',18)
set(gca,'FontName','Arial')
grid on
title('Step Resp. to PCs','FontName','Arial','FontSize',18); xlabel('Months','FontName','Arial','FontSize',18); ylabel(strcat(target_id,', ',target_units),'FontName','Arial','FontSize',18); grid on; hold off

filename = strcat(fig_dir,save_id,'_step_resps.fig')
saveas(gcf,filename)

figure
eshade2(1:length(CR_sum), detrend(CR_sum)-mean(detrend(CR_sum)), 2*CR_sum_errs, 'c')
hold on
eshade2(1:length(CR_sum),detrend(CR_sum)-mean(detrend(CR_sum)), CR_sum_errs, 'b')
plot(detrend(CR_sum)-mean(detrend(CR_sum)), 'g', 'LineWidth', 2.6)
plot(detrend(target_ds)-mean(detrend(target_ds)),'r','LineWidth',2.8);hold on
xlim([420-1 3778])
set(gca,'FontSize',18)
set(gca,'FontName','Arial')
title('Reconstruction using PC1, PC2, PC3','FontName','Arial','FontSize',18)
xlabel('Months','FontName','Arial','FontSize',18); ylabel('FW flux [m^3/s]','FontName','Arial','FontSize',18); grid on; hold off

filename = strcat(fig_dir,save_id,'_recon_summed.fig')
saveas(gcf,filename)


y1 = S_ERA(:,1);
y2 = S_ERA(:,2);
y3 = S_ERA(:,3);
y4 = sum(S_ERA,2);

e1 = E_ERA(:,1);
e2 = E_ERA(:,2);
e3 = E_ERA(:,3);
e4 = (e1.^2 + e2.^2 + e3.^2).^0.5;

% time: 1900 onwards

T = zeros(1,t+t_plus);
for i = 1:t+t_plus
    T(i) = 1900 + (i-1)/12;
end

% values before 360 months are less reliable

% xl = (1900+30)*ones(100,1);
% yl = linspace(-0.5*10^13,0.5*10^13);
filename = strcat(var_dir,save_id,'_ERA_20Ci_recon.mat');
save(filename,'y1','y2','y3','y4','e1','e2','e3','e4'); 

PC1_36 = [movmean(P_ERA(:,1),36);zeros(t_plus,1)];
PC2_36 = [movmean(P_ERA(:,2),36);zeros(t_plus,1)];
PC3_36 = [movmean(P_ERA(:,3),36);zeros(t_plus,1)];
PC1_36(PC1_36 == 0) = NaN;
PC2_36(PC2_36 == 0) = NaN;
PC3_36(PC3_36 == 0) = NaN;
PC_sum = abs(PC1_36) + abs(PC2_36) + abs(PC3_36);

G1 = mean_G(:,1);
G2 = mean_G(:,2);
G3 = mean_G(:,3);

GE1 = impulse_err(:,1);
GE2 = impulse_err(:,2);
GE3 = impulse_err(:,3);

filename = strcat(var_dir,save_id,'_plot_impulse.mat');
save(filename,'G1','G2','G3','GE1','GE2','GE3'); 

figure
subplot(3,1,1)
errorbar(T, y1, e1, 'Marker', '.', 'MarkerSize', 24);
hold on
plot(T, y1, 'Color', [0 0.5 0], 'linewidth', 3)
% line(xl,yl, 'color', [0.5 0 1], 'linewidth', 1); plot([1900 2017],[0 0]);
xlim([1900 2041])
ylabel(strcat(target_id,', ',target_units))
title(strcat('Predicted ',title_id,' response to ERA-20Ci, eof 1'))
hold off
subplot(3,1,2)
plot(T,PC1_36); hold on; 
ylabel('PC strength, 36-month mov-mean')
plot([1900 2041],[0 0]); hold off;
xlim([1900 2041])
xlabel(tseries_units)
subplot(3,1,3)
errorbar(G1,GE1,'Marker','.','MarkerSize',18); hold on;
xlim([0 1692])
ylabel('impulse response')
xlabel('months')
hold off;
filename = strcat(fig_dir,save_id,'_ERA_20Ci_pred_eof1.fig');
saveas(gcf,filename)

figure
subplot(3,1,1)
errorbar(T, y2, e2, 'Marker', '.', 'MarkerSize', 24);
hold on
plot(T, y2, 'Color', [0 0.5 0], 'linewidth', 3)
% line(xl,yl, 'color', [0.5 0 1], 'linewidth', 1); plot([1900 2010],[0 0]);
xlim([1900 2041])
ylabel(strcat(target_id,', ',target_units))
title(strcat('Predicted ',title_id,' response to ERA-20Ci, eof 2'))
hold off
subplot(3,1,2)
plot(T,PC2_36); hold on; 
ylabel('PC strength, 36-month mov-mean')
plot([1900 2041],[0 0]); hold off;
xlim([1900 2041])
xlabel('time')
subplot(3,1,3)
errorbar(G2,GE2,'Marker','.','MarkerSize',18); hold on;
xlim([0 1692])
ylabel('impulse response')
xlabel(tseries_units)
hold off;
filename = strcat(fig_dir,save_id,'_ERA_20Ci_pred_eof2.fig');
saveas(gcf,filename)

figure
subplot(3,1,1)
errorbar(T, y3, e3, 'Marker', '.', 'MarkerSize', 24);
hold on
plot(T, y3, 'Color', [0 0.5 0], 'linewidth', 3)
% line(xl,yl, 'color', [0.5 0 1], 'linewidth', 1); plot([1900 2010],[0 0]);
xlim([1900 2041])
ylabel('fw flux, m^3/s')
title(strcat('Predicted ',title_id,' response to ERA-20Ci, eof 3'))
hold off
subplot(3,1,2)
plot(T,PC3_36); hold on; 
ylabel('PC strength, 36-month mov-mean')
plot([1900 2041],[0 0]); hold off;
xlim([1900 2041])
xlabel('time')
subplot(3,1,3)
errorbar(G3,GE3,'Marker','.','MarkerSize',18); hold on;
xlim([0 1692])
ylabel('impulse response')
xlabel(tseries_units)
hold off;
filename = strcat(fig_dir,save_id,'_ERA_20Ci_pred_eof3.fig');
saveas(gcf,filename)


figure

errorbar(T, y4, e4, 'Marker', '.', 'MarkerSize', 24);
hold on
plot(T, y4, 'Color', [0 0.5 0], 'linewidth', 3)
% line(xl,yl, 'color', [0.5 0 1], 'linewidth', 1); plot([1900 2010],[0 0]);
hold on
% plot(ref_data(:,1),ref_data(:,2)-0.5*10^12,'r*')
xlim([1900 2041])
ylabel(strcat(target_id,', ',target_units))
title(strcat('Predicted ',title_id,' response to ERA-20Ci, all eofs'))
hold off
xlabel('time')
filename = '/home/ocean2/samc/MO_CRFs/figures/FWC_anude_ERA_20Ci_pred_all_eofs.fig';
saveas(gcf,filename)


figure
suptitle('Impulse Response Functions')
xl = zeros(360,1);
xt = 1:360;
subplot(3,1,1)
eshade2(xt,f1.*mean_G(:,1),impulse_err(:,1),[1 0.6 0.4])
hold on
plot(xt,f1.*mean_G(:,1),'k','linewidth',2)
plot(xt,xl,'k')
% ylim([-2*10^11 3*10^11])
xlim([0 370])
xlabel(tseries_units)
ylabel(strcat(target_id,', ',target_units))
hold off
legend('PC1')
subplot(3,1,2)
eshade2(xt,f2.*mean_G(:,2),impulse_err(:,2),[0.5 1 0.4])
hold on
plot(xt,f2.*mean_G(:,2),'k','linewidth',2)
plot(xt,xl,'k')
% ylim([-2*10^11 3*10^11])
xlim([0 370])
xlabel(tseries_units)
ylabel(strcat(target_id,', ',target_units))
hold off
legend('PC2')
subplot(3,1,3)
eshade2(xt,f3.*mean_G(:,3),impulse_err(:,3),[0.8 0.8 0.8])
hold on
plot(xt,f3.*mean_G(:,3),'k','linewidth',2)
plot(xt,xl,'k')
% ylim([-2*10^11 3*10^11])
xlim([0 370])
xlabel(tseries_units)
ylabel(strcat(target_id,', ',target_units))
hold off
legend('PC3')

filename = strcat(fig_dir,save_id,'impulse_responses.fig');
saveas(gcf,filename)


