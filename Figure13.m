%% Test 2 - predicted Mmax given b-value
clear

%% Load data
data = readtable('Basel_Data.csv'); 
mws = data.Mw; 
t = datenum(data.DateTime); 

% shut-in date
injStopDate = datenum(2006,12,8,11,33,0); % this is the actual shut-in time
index = t <= injStopDate;
nindex = logical(1-index);

%% number of events before and after shut-in

% completeness threshold
Mc = 0.9;

% get events above Mc before and after shut-in
InjEvents = mws(index);
PostInjEvents = mws(nindex);
InjEvents(isnan(InjEvents))=[];
InjEvents = InjEvents(InjEvents >= Mc); 
PostInjEvents(isnan(PostInjEvents)) = []; 
PostInjEvents = PostInjEvents(PostInjEvents>=Mc);

Ninj = numel(InjEvents); 
Npostinj = numel(PostInjEvents); 

%% Estimated b-values
% Use the Tinti-Malargia b-value (same as Aki MLE but adjusted for discrete
% magnitude grid size
binj_tm = bval_TM(InjEvents(InjEvents>Mc)); 

% bootstrap b-value to get uncertainty
% b-value for all events 
bootdat = InjEvents(InjEvents>Mc); 
bootstat = bootstrp(1000, @bval_TM, bootdat(:)); 

% pre- and post- shut-in events separately
bootdat_post = PostInjEvents(PostInjEvents>Mc); 
bootstat2 = bootstrp(1000, @bval_TM, bootdat_post(:)); 

%% Compute theoretical Mmax for various b-values
Mmax_obs = max(mws); 
bvals = [0.6:.05:1.7]; 
[Mmaxpi, Mmaxlpi, Mmaxupi] = deal(zeros(length(bvals), 1)); 
for loop = 1:length(bvals)
    b = bvals(loop); 
    [Mmaxpi(loop), Mmaxlpi(loop), Mmaxupi(loop)] = compute_Mmax_GR (Mc,b, Npostinj, 95);
end

%% Plot the results
figure; 
plot(bvals, Mmaxpi, '--k', 'LineWidth', 1)
hold on
plot(bvals, Mmaxlpi, '--g', 'LineWidth', 1)
plot(bvals, Mmaxupi, '--r', 'LineWidth', 1)
refline(0, Mmax_obs)
histogram(bootstat, 'normalization', 'pdf')
histogram(bootstat2, 'normalization', 'pdf')
xlabel('b-value')
ylabel('Mmax')
legend('Mean', 'Lower 95% bound', 'Upper 95% bound', 'Observed', 'Pre shut-in', 'Post shut-in')
xlim([min(bvals), max(bvals)])
set(gca, 'FontSize', 16)


