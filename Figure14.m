%% Test 3 - Probability of N_M events >= M given N, b
% This script produces the full probability distribution of observing N_M
% events with magnitude >= M given the total number of events, N, assuming
% GRD with a specific b-value. 

% magnitude M
Mtest = 2.5;

b_pre = 1.56;
b_post = 1; 

% pick which b-values to use
b = b_pre; 

%% Load data
data = readtable('Basel_Data.csv'); 
mws = data.Mw; 
t = datenum(data.DateTime); 

% shut-in date
injStopDate = datenum(2006,12,8,11,33,0); % this is the actual shut-in time

% completeness threshold
% get events above Mc before and after shut-in
Mc = 0.9;
Mc_ind = mws >= Mc; 
AllEvents = mws(Mc_ind);
Allt = t(Mc_ind);
tindex = Allt  > injStopDate;
PostInjEvents = AllEvents(tindex);
Npostinj = numel(PostInjEvents); 

%% Numerical simulation
% simulate a bunch of seismicity sequences that follow GRD using the given
% b-value, and count the number of simulations that contain at least N_M
% events >= M, given the number of events that occurred after shut-in. 
Niter = 1e4; 
Ntot = Npostinj; 

% initialize and run simulations
n_M = zeros(Niter, 1); 
nbins = 20; 
tindex = zeros(Niter, 1); 
for loop = 1:Niter
    
    % simulate Ntot total events from GRD using b-value b and Mc
    simmags = zeros(Ntot, 1); 
    for k = 1:Ntot
        simmags(k) = rs2mw(genRsGR(mw2rs(4), mw2rs(Mc), b)); 
    end
    % accumulate events >= M
    n_M(loop) = sum(simmags>=Mtest); 
    
    % demonstrate that the index of Mmax is uniform, consistent with van
    % der Elst et al., 2017
    tindex(loop) = find(simmags==max(simmags)); 
end

%% Analytically compute the pdf, using the binomial distribution
Nms = [0:15];
pMgm = @(m, mc,b) 10.^(-b*(m-mc));
prob2p5 = pMgm(2.5, Mc, b);
probNm = binopdf(Nms, Ntot, prob2p5);

%% Plot the results
figure; 
histogram(n_M, 'normalization', 'pdf')
xlabel('# of Events post-shut-in')
ylabel('Probability')
hold on
plot(Nms, probNm)
% vline(sum(mws>=2.5))
vline(sum(PostInjEvents>=2.5), '--');
axis tight
legend('Simulated','Theoretical', 'Observed')

figure;
histogram(tindex(:)./Ntot, 'Normalization', 'pdf')
xlabel('Relative Sequence Location of Mmax')
title('Location of Mmax in random synthetic GR catalogs')
