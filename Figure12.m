%% Test 1 for time-dependent FMD: Mw vs event number

clear

% Load data
data = readtable('Basel_Data.csv'); 
mws = data.Mw; 

% plot
figure; 
scatter([1:length(mws)], mws, '.k')
xlabel('Event #')
ylabel('M_W')
title("Changing FMD Test 1")
