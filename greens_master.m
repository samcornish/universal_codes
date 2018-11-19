% This is a universal script for computing impulse responses, step responses, and reconstructions.
% It uses a statistical technique developed by Kostov et al. (2017), and used by Johnson et al. (2018).
% The basis of the technique is to treat a target timeseries (or response series) as a convolution of a known forcing series and a given number of unknown impulse response functions.
% We solve for the impulse responses using multiple lagged regression.
% Step responses are obtained by integration of the impulse responses through all lags.
% Reconstructions of the target series are obtained by convolving the derived impulse responses with the original forcing.
% Many estimates are made for the impulse responses, to obtain a good understanding of the error and to minimise over-fitting.
% These estimates come from changing the impulse cutoff lag, and scanning through the target series using multiple segments of a shorter length than the whole target series 
% Error is calculated by combining in quadrature the estimates of error derived from the residual between the mean reconstruction and target ("residual") and those from the standard deviation of the estimates ("spread").
% This code was written by Sam B Cornish, 2018. It is based on precursors written by Emma Beer and Yavor K. Kostov.
% Written to be run with Matlab 2016a, may need updating for subsequent editions

% This script allows the user to specify:
% One target timeseries
% A given number of forcing timeseries
% The length of the derived impulses (the number of lags to compute the regression over). The target must be at least twice this length to avoid over-fitting.
% The range of lengths of impulses (to generate more estimates). The minimum length will always be that which is taken forwards.
% The number of different segment schemes used to scan through the target timeseries and generate more estimates.
% The length of each segment scheme.
% The spacing between segments used in the segment scheme.


% Code status: in progress

%% Loading variables and intial processing - user should edit in this section only 
clear
close all

% addpaths here to the location of the variables
% ------------------------------------------------- % vv
addpath /home/ocean2/samc/MO_CRFs/variables/
% addpath to the location of the function eshade2.m 
addpath /home/ocean2/samc/HiGEMArcticCRFs/fw_content_atm_circ/analysis/functions/
% ------------------------------------------------- % ^^

% load forcing timeseries here. This should probably be deseasonalised, if using climate data
% ------------------------------------------------- % vv
load pressure_EOFs_anude.mat
forcing = pc(:,1:end);	% replace variable name as appropriate. 
%IMPORTANT: this should be a vector or a matrix containing the forcing along separate columns.
% ------------------------------------------------- % ^^
% Here you have the option to specify the required number of components
% (separate columns) you want to include in the forcing (replace 'end'
% accordingly).
sdf = std(forcing); % be aware that std will take the standard deviation of each column
szf = size(forcing);
forcing = detrend(forcing,'constant');   % remove the mean from each component. You may wish to use a linear detrend.

% normalise by the standard deviation of each component of the forcing. We
% do this assuming we want to find the response to a one standard deviation
% step in the forcing
% ------------------------------------------------- % vv
% IMPORTANT: comment out these lines if you don't want to normalise the
% forcing by its standard deviation
for i = 1:szf(2)
    forcing(:,i) = forcing(:,i)./sdf(i);
end
% ------------------------------------------------- % ^^

% Load the target series here
% ------------------------------------------------- % vv
load fw_tseries_1_anude.mat fw_ts_ds  % this should probably be deseasonalised, if using climate data
target = fw_ts_ds;  % replace var_name as appropriate
% ------------------------------------------------- % ^^
target = detrend(target,'linear');   % decide if a linear or constant detrend is most appropriate 

% Determining the strings used to save files
% saving variables
% insert the directory name for saving the variables
% ------------------------------------------------- % vv
var_dir = '/home/ocean_personal_data/samc/MO_CRFs/variables/';
% directory for figures:
fig_dir = '/home/ocean_personal_data/samc/MO_CRFs/figures/';
% ------------------------------------------------- % ^^

% filename identifier - this should be a word or several that uniquely
% identifies the saved variables to the input variables. Usually this will
% be the name of the target timeseries and the model or record it is from
% ------------------------------------------------- % vv
save_id = 'FWC_anude';
title_id = 'MO GC2 FWC';    % for the sake of plotting (spaces are fine!)
target_id = 'FWC';          % for y axis label (spaces are fine!)
target_units = 'm^3';
tseries_units = 'months';
% ------------------------------------------------- % ^^

%% Determine cutoff and segment choices
% Choices for the length of memory going back a given number of lagged
% months. We vary the cutoff lag to obtain independent estimates.
tau_cutoff_choices= 360:6:420;     % Change these as appropriate
N = max(tau_cutoff_choices) + 100; % a nominal large number, always larger than the greatest value of tau cutoff. +100 probably unecessary. 
forcing = [zeros(N,szf(2));forcing];  % This should append N zeros to the start of each column of forcing.

number_string = {'one','two','three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten'};  % add more if necessary
comp_no = number_string(1:szf(2));  % component number - these cells will be used later

L = length(target);

warning_limit = 2*max(tau_cutoff_choices);  % This is a warning when the memory you are trying to capture is close to the length of the target timeseries (this will result in over-fitting). 
if L <= warning_limit
    fprintf('length of target is probably insufficient for the chosen cutoff tau. Reconsider, or change the warning limit')
end   
% segment choices - store in a 3d array 
% L - (no_segs-1)*spacing = length_segs
no_segs = 100;  % number of segments, this must be chosen a priori
% it is usually best to choose the spacing of the segments a priori,
% although everything really depends on the relationship between cutoff
% lags used and L.
Nsp = 1000; % max possible value of spacing. 1000 is an arbitrary high value, that would only be used with very long L. Consider changing if L is order 1e7 or above.
l_seg = zeros(Nsp,1);
for poss_sp = 1:Nsp    % possible values of the spacing. 
    l_seg(poss_sp) = L - (no_segs-1)*poss_sp;
    if l_seg(poss_sp) <= warning_limit
        break
    end
end
loc_vals = find(l_seg);    % locate non-zero values
l_seg = l_seg(loc_vals);  % omit zero values
lls = length(l_seg);
l_seg_choices = [L, l_seg(1), l_seg(round(lls/3)), l_seg(round(lls*2/3))];   % sample different possible options of segment lengths
% we also want to include the case where the whole run is sampled at once, hence inclusion of L as first value (spacing = 0).
spacing = (L - l_seg_choices)./(no_segs-1);    
    
l_sp=length(spacing);
spacing_index = 1:l_sp;

all_start_times = zeros(no_segs,l_sp);
all_end_times = zeros(no_segs,l_sp);
segment = zeros(L,no_segs,l_sp); 
for n = 1:l_sp
for nn=1:no_segs
     all_start_times(nn,n)=(nn-1)*spacing(n)+1; % all_start_times gives the starting indices of each segment
end
 all_end_times(:,n)=all_start_times(:,n)+L-all_start_times(end,n);

 for i = 1:no_segs
     segment(1:l_seg_choices(n),i,n) = all_start_times(i,n):all_end_times(i,n);
 end
end
sz_seg = size(segment);

%% Impulse and step response calculation

% Building a five-dimensional array of all estimates for the impulse, G
G_5D = zeros(min(tau_cutoff_choices),length(tau_cutoff_choices),sz_seg(2),l_sp,szf(2)); 
sz_G = size(G_5D);
G_4D = zeros(sz_G(1),sz_G(2)*sz_G(3),sz_G(4),szf(2));    % later on we will reshape the G array, combining estimates from changing the cutoff lag and from different segments within a segment scheme into one dimension
mean_G_all = zeros(sz_G(1),l_sp,szf(2));    % we will then collaps that dimension by taking the mean
std_G_all = zeros(sz_G(1),l_sp,szf(2));
fc = 1; % Forcing component
ssn = 1;    % Segment scheme number

fprintf('starting impulse calculations')
for e = 1:szf(2)    % indexing for the component of the forcing used
for k = spacing_index     % indexing for the segment scheme used
    l = l_seg_choices(k);       % ensuring that we take into account the fact that the segments have different lengths between the segments schemes
    for s = 1:sz_seg(2) % indexing through each segment s in scheme k
        for tau_index = 1:length(tau_cutoff_choices)    % indexing for the lag 
            
            tau_cut = tau_cutoff_choices(tau_index);
            
            clearvars X Y M i G Y   % This might be unnecessary 
            
            X = forcing(all_start_times(s,k)+N-tau_cut+1:all_end_times(s,k)+N,e);    % We use X for simplicity; it represents the forcing series in what follows

            Y = target(segment(1:l_seg_choices(k),s,k));    % We use Y for simplicity; it represents the target series in what follows
            
            M = zeros(tau_cut,l);     %  M is a matrix of the forcing, with no. of rows = number of lags tau used, and no. of columns = total length l of the target
            for i = 1:tau_cut
                M(i,:) = X(tau_cut-(i-1):l+tau_cut-i);  % each row i has l elements
            end

            G = M'\(Y'); % regression of the freshwater onto the forcing, lagged (from M)
            G_5D(:,tau_index,s,k,e) = G(1:min(tau_cutoff_choices));
        end
    end
    fprintf('completed segment scheme number #%d\n', ssn)
        ssn = ssn +1;
end


G_4D(:,:,:,e) = reshape(G_5D(:,:,:,:,e),sz_G(1),sz_G(2)*sz_G(3),sz_G(4));

std_G_all(:,:,e) = squeeze(std(G_4D(:,:,:,e),[],2));
mean_G_all(:,:,e) = squeeze(mean(G_4D(:,:,:,e),2));
fprintf('finished forcing component no.#%d\n', fc)
fc = fc +1;
end

% compute step responses
step_resp = zeros(min(tau_cutoff_choices),l_sp,szf(2));
step_err_spread = zeros(min(tau_cutoff_choices),l_sp,szf(2));
for e = 1:szf(2)
for k = spacing_index
for j = 1:min(tau_cutoff_choices)
            step_resp(j,k,e) = sum(mean_G_all(1:j,k,e));		
            step_err_spread(j,k,e) = sum(std_G_all(1:j,k,e));
end
end
end
fprintf('step responses complete')
%% Reconstruction and error calc - individual segment schemes


tau_cut = min(tau_cutoff_choices);    % choose the minimimum cutoff choice to do reconstructions
step_errs= zeros(tau_cut,l_sp,szf(2));
Recon = zeros(L,l_sp,szf(2));   
Recon_err= zeros(L,l_sp,szf(2));
impulse_err= zeros(tau_cut,l_sp,szf(2));
ssn = 1;
fc = 1;

for e = 1:szf(2)
for k = spacing_index
    
        clearvars X Y M sd G epsilon var G_var S2 recon_ctrl s2...
            recon_ctrl_errs_spread recon_errs_res impulse_err_spread...
            step_err_res X1 i

        X1 = forcing(N+1:end,e); % Now need the whole target series to compute the reconstruction
       
        Y = target;
        
        X = [zeros(tau_cut-1,1);X1]; % The column vector X with tau_cut - 1 zeros appended to start
        
        M = zeros(tau_cut,L);     %  initialising a matrix of the forcing, with tau_cut rows and L columns
        for i = 1:tau_cut
            M(i,:) = X(tau_cut-(i-1):L+tau_cut-i);  
        end
        
        recon_ctrl = mean_G_all(:,k,e)'*M;
        recon_ctrl_errs_spread = std_G_all(:,k,e)'*M;
        
        epsilon = Y-mean_G_all(:,k,e)'*M; %errors
        s2 = (epsilon*epsilon')/(L-tau_cut); % sigma^2 = sum [epsilon] / (n-p) , n-p = degress of freedom = no. of data points minus number of impulse response coefficients = L - tau_cut

        S2 = s2*eye(tau_cut);
        G_var = M*M'\S2; %variance-covariance matrix for the impulse response to a given forcing component
        
        impulse_err_spread = diag(G_var(1:tau_cut,1:tau_cut)).^0.5;
        step_err_res = zeros(tau_cut,1);
        for j = 1:tau_cut
            step_err_res(j) = sum(sum(G_var(1:j, 1:j))).^0.5;
        end
        recon_errs_res = ( diag(  M' * G_var * M ) ) .^0.5 ;
      
        step_errs(:,k,e) = (step_err_res.^2+step_err_spread(:,k,e).^2).^0.5; % The error takes into account the residual from the fits and the spread across the ensemble of fits.

        Recon(:,k,e) = recon_ctrl;
        Recon_err(:,k,e)  = (recon_errs_res'.^2 + recon_ctrl_errs_spread.^2).^0.5;

        impulse_err(:,k,e) = (impulse_err_spread.^2 + std_G_all(:,k,e).^2).^0.5;

        fprintf('seg scheme error calc number complete: no. #%d\n', ssn)
            ssn = ssn +1;
    
end
fprintf('forcing component complete: no. #%d\n', fc)
            fc = fc +1;
end
%% Collating estimates for G and calculating 'grand' step responses. pt.1/2 error calculation.
% error calculation for collated G
G_coll = zeros(sz_G(1),sz_G(2)*sz_G(3)*sz_G(4),szf(2)); % G collated across all estimates, where we change: cutoff lag, segment, segment scheme.
G_grand_mean = zeros(sz_G(1),szf(2));    % the grand mean estimate of G for each forcing component
G_grand_std = zeros(sz_G(1),szf(2));     % the grand std component of error for G for each forcing component
grand_SR = zeros(tau_cut,szf(2));        % grand step response
grand_SE_spread = zeros(tau_cut,szf(2)); % grand step error: spread component 
fc = 1;   

for e = 1:szf(2)
G_coll(:,:,e) = reshape(G_4D(:,:,:,e),sz_G(1),sz_G(2)*sz_G(3)*sz_G(4)); 
G_grand_mean(:,e) = squeeze(mean(G_coll(:,:,e),2));
G_grand_std(:,e) = squeeze(std(G_coll(:,:,e),[],2)); 
figure 
plot(G_grand_mean(:,e))

for j = 1:tau_cut
            grand_SR(j,e) = sum(G_grand_mean(1:j,e));		
            grand_SE_spread(j,e) = sum(G_grand_std(1:j,e));
end
figure 
plot(grand_SR(:,e))
fprintf('error calc collated G pt.1 forcing component no. #%d\n', fc)
            fc = fc +1;
end
%% error calculation for collated G (2/2)
% just do one set of calculations for each forcing component

grand_CR = zeros(L,szf(2));
grand_CR_errs = zeros(L,szf(2));
grand_SE = zeros(tau_cut,szf(2));
grand_impulse_err = zeros(tau_cut,szf(2));
fc = 1;
G_var_save = zeros(tau_cut,tau_cut,szf(2));

for e = 1:szf(2)
    clearvars X X1 Y M sd G epsilon var G_var i S2 resp_ctrl s2 grand_IE_spread grand_CR_err_res grand_SE_res grand_CR_errs_spread
        grand_SE_res = zeros(tau_cut,1);
        X1 = forcing(N+1:end,e); % Now need the whole target series
        Y = target;
        X = [zeros(tau_cut-1,1);X1]; % The column vector X (forcing) with tau_cut - 1 zeros appended to start
        
        M = zeros(tau_cut,L);     %  initialising a matrix of the forcing, with tau_cut rows and L columns
        for i = 1:tau_cut
            M(i,:) = X(tau_cut-(i-1):L+tau_cut-i);  
        end
        
        grand_CR(:,e) = G_grand_mean(:,e)'*M;         % reconstruction
        grand_CR_errs_spread = G_grand_std(:,e)'*M;       % recon errors: spread component
        
        epsilon = Y-G_grand_mean(:,e)'*M; % residual
        s2 = (epsilon*epsilon')/(L-tau_cut); % sigma^2 = sum [epsilon] / (n-p) , n-p = degress of freedom = no. of data points minus number of impulse response coefficients = t-g

        S2 = s2*eye(tau_cut);
        G_var = M*M'\S2; % variance-covariance matrix for the impulse response to a given forcing component
        G_var_save(:,:,e) = G_var;

        grand_IE_spread = diag(G_var(1:tau_cut,1:tau_cut)).^0.5;       % impulse error spread component
        
        for j = 1:tau_cut
            grand_SE_res(j) = sum(sum(G_var(1:j, 1:j))).^0.5;        % step error residual component
        end
        grand_CR_errs_res = ( diag(  M' * G_var * M ) ) .^0.5 ;  % recon errors: residual component
      
        grand_SE(:,e) = (grand_SE_res.^2+grand_SE_spread(:,e).^2).^0.5; % step errors from res and spread. The error takes into account the residual from the fits and the spread across the ensemble of fits.

        grand_CR_errs(:,e) = (grand_CR_errs_res.^2 + grand_CR_errs_spread'.^2).^0.5;    % total recon errors

        grand_impulse_err(:,e) = (grand_IE_spread.^2 + G_grand_std(:,e).^2).^0.5;
        fprintf('error calc collated G pt.2 forcing component complete: #%d\n', fc)
            fc = fc +1;
end

fprintf('plots and saving')
%% plots and saving
%plots

figure
for e = 1:szf(2)
for k = spacing_index
plot(step_resp(:,k,e))
hold on
end
end
legend('eof 1','eof 2','eof 3')
filename = strcat(var_dir,save_id,'_G.mat');
save(filename,'G_5D','std_G_all','mean_G_all','G_grand_mean','G_grand_std','grand_impulse_err','G_var_save');
filename = strcat(var_dir,save_id,'_step_resp.mat');
save(filename,'step_resp', 'step_errs', 'grand_SE', 'grand_SR')
filename = strcat(var_dir,save_id,'_recon.mat');
save(filename,'target','grand_CR', 'grand_CR_errs', 'Recon', 'Recon_err')

leg_info = cell(1,szf(2));
for i=1:szf(2)   % Looping through the forcing components
figure; %individual segment schemes
hold on
for s = fliplr(spacing_index)
errorbar(1:length(step_resp(:,s,i)), step_resp(:,s,i), step_errs(:,s,i));
leg_info{s} = ['spacing = ' num2str(spacing(s))];
end
hold off;
legend(leg_info)
xlim([0 1.1*min(tau_cutoff_choices)])
set(gca,'FontSize',18)
set(gca,'FontName','Arial')
xlabel(tseries_units,'FontName','Arial','FontSize',18); ylabel(strcat(target_id,' change, ',target_units),'FontName','Arial','FontSize',18); grid on; hold off
name = strcat(title_id,' response to forcing component ',comp_no{i}); 
title(name);
filename = strcat(fig_dir,save_id,'_step_resp_segments_forcing_comp_',comp_no{i});
saveas(gcf, filename)
clearvars leg_info
figure % collated segments 
errorbar(1:length(grand_SR(:,i)), grand_SR(:,i), grand_SE(:,i), 'b'); hold on;
plot(grand_SR(:,i),'k'); hold off
set(gca,'FontSize',18)
set(gca,'FontName','Arial')
name = strcat(title_id,' step resp to forcing component: ',comp_no{i});
title(name,'FontName','Arial','FontSize',18); xlabel(tseries_units,'FontName','Arial','FontSize',18); ylabel(strcat(target_id,' change, ',target_units),'FontName','Arial','FontSize',18); grid on; 
filename = strcat(fig_dir,save_id,'step_resp_collated_forcing_comp_',comp_no{i});
saveas(gcf, filename)

figure  % collated recon 
eshade2(1:length(grand_CR(:,i)), grand_CR(:,i), grand_CR_errs(:,i), 'b')
hold on
plot(grand_CR(:,i), 'g', 'LineWidth', 2.6)
plot(target,'r','LineWidth',2.8);hold on
set(gca,'FontSize',18)
set(gca,'FontName','Arial')
name = strcat(title_id,' recon, forcing comp: ',comp_no{i});
title(name,'FontName','Arial','FontSize',18)
xlabel(tseries_units,'FontName','Arial','FontSize',18); ylabel(strcat(target_id,' change, ',target_units),'FontName','Arial','FontSize',18); grid on; hold off

filename = strcat(fig_dir,save_id,'_recon_fc_',comp_no{i});
saveas(gcf, filename)
end

figure % step responses to all forcing components in one plot 

for f = fliplr(1:szf(2))    % we do this in reverse order so that the leading forcing components are at the front of the plot
errorbar(1:length(grand_SR(:,f)), grand_SR(:,f), grand_SE(:,f)); hold on;
leg_info{f} = ['forcing component ',num2str(f)];
plot(grand_SR(:,f),'k','LineWidth',2);
end
legend(leg_info)
xlabel(tseries_units); ylabel(strcat(target_id,', ',target_units)); title(strcat(title_id,' step responses')); xlim([0 1.1*min(tau_cutoff_choices)])

filename = strcat(fig_dir,save_id,'_step_resps');
saveas(gcf,filename)


%% 
% correlations

R_fc_seg = zeros(l_sp,szf(2));
P_fc_seg = zeros(l_sp,szf(2));
figure
for i = 1:szf(2)
for k = spacing_index
[r,p] = corrcoef(detrend(target(1,360:end)),detrend(Recon(360:end,k,i)));
    R_fc_seg(k,i) = r(1,2);
    P_fc_seg(k,i) = p(1,2);
end

plot(R_fc_seg(:,i),'--o')
hold on 
% plot(P(:,i))


xlim([1 7])
xlabel('segment scheme')
name = 'target/recon correlations across forcing components';
title(name,'FontName','Arial','FontSize',18)
ylabel('R')

[r_fc, p_fc] = corrcoef(detrend(target(1,360:end)),detrend(grand_CR(360:end,i)));

end
[r_overall, p_overall] = corrcoef(detrend(target(1,360:end)),detrend(sum(grand_CR(360:end,:),2)));
legend(comp_no)
filename = strcat(fig_dir,save_id,'_recon_correlations');
saveas(gcf, filename)
filename = strcat(var_dir,save_id,'_recon_corrcoefs.mat');
save(filename,'R_fc_seg','P_fc_seg','r_fc','p_fc','r_overall','p_overall')
