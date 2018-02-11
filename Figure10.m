%% Script for testing on Paul's numerical simulations

clear
clc

%% Load data, set parameters
load EqMagData_fromSegallLupaper.mat
times = s.tt; 
a_of_t = sqrt(4*s.props.c*s.tt);
dtau = 3e6;
shear_mod = s.props.G; 
T = s.T; 
b = 1.5; 

%% compute and normalize the rate of events
% note that there should be a dV here, but because we will normalize it
% actually doesn't matter
Rt = sum(s.Rr,2); %--> rate of events at each time
ntot= sum(Rt); 

% normalize to the desired total number of events between a range of
% magnitudes
Mwmin = 0; 
Mwmax = 6; 

Ndes = 1e4; 
Rt = round(Rt./(ntot/Ndes)); 
maxRate = 2^(nextpow2(max(Rt)));
NumEventsTotal = zeros(length(times),1); 
for loop = 2:length(times)
    NumEventsTotal(loop) = trapz(times(1:loop), Rt(1:loop)); 
end

%% Compute event times and perturbed region radius with time
% compute event times for each time using the rate of events, also compute
% the size of the perturbed region at each set of event times.
eventtimes = []; avec = [];
for loop = 1:length(times)-1
   t = linspace(times(loop), times(loop+1), maxRate); 
   tvec = sort(randi(length(t), 1,Rt(loop))); 
   t2add = t(tvec); 
   eventtimes = [eventtimes, t2add];
   avec = [avec, interp1([min(t), max(t)], [a_of_t(loop), a_of_t(loop+1)], t2add)];
end

% plot basic stuff
% figure; plot(times, Rt)
% figure; plot(times, a_of_t)

N = length(eventtimes); 
% figure; scatter(eventtimes,avec, '.')

%% Generate from GRD
events_GRD = zeros(N,1); 
rmin = mw2rs(Mwmin, dtau); 
rmax = mw2rs(Mwmax, dtau); 
for k = 1:N
    events_GRD(k) = genRsGR(rmax, rmin, b); 
end
Mws_GRD = rs2mw(events_GRD, dtau); 

%% Generate from Maurer-Segall
% Run the loop
Nx = 600; 
Mws_JP = zeros(N,1);
[Rx, Prs] = deal(zeros(Nx, N)); 
for loop = 1:N
    a = avec(loop); 
    rhomin = mw2rs(Mwmin, dtau)./a; 
    rhomax = mw2rs(Mwmax, dtau)./a; 
    
    % compute P(r_s)
%     [Frs, rx, Prs(:,loop)]  = Compute_Frs_an(rhomax, rhomin,b, Nx);
    [Prs(:,loop), rx, Frs]  = Compute_Prs(b, rhomin, rhomax, Nx);
    
    % generate random samples
    Mws_JP(loop) = genEqsPrs(rx, Frs, 1,a, dtau);     
    Rx(:,loop) = rx.*a; 
end

%% Compute pre- and post-shut-in FMDs

% bins
Nbins = 20; 

% separate between events before and after shut-in
shutinindex = eventtimes < T; 

% GRD pre-shut-in
[hGRp,xGRp] = hist(Mws_GRD(shutinindex),Nbins);

% Maurer-Segall post-shut-in
[hJPp] = hist(Mws_JP(shutinindex),xGRp); 

% theoretical Maurer-Segall pre-shut-in
[Prs_totp,RXp] = compute_timeintegrated_prs(Rx(:,shutinindex),Prs(:,shutinindex),eventtimes(shutinindex)); 
[mwsp,fMwp] = rspdf2mwpdf(RXp, Prs_totp, dtau, 7*pi/16);
ind = find(abs(mwsp-min(xGRp)) == min(abs(mwsp-min(xGRp))));
fMwp = fMwp./(fMwp(ind)/max(hJPp)); 

% GRD post-shut-in
[hGRpo,xGRpo] = hist(Mws_GRD(~shutinindex),Nbins);

% Maurer-Segall post-shut-in
[hJPpo] = hist(Mws_JP(~shutinindex),xGRpo); 

% theoretical Maurer-Segall after shut-in
[Prs_totpo,RXpo] = compute_timeintegrated_prs(Rx(:,~shutinindex),Prs(:,~shutinindex),eventtimes(~shutinindex)); 
[mwspo,fMwpo] = rspdf2mwpdf(RXpo, Prs_totpo, dtau, 7*pi/16);
ind = find(abs(mwspo-min(xGRpo)) == min(abs(mwspo-min(xGRpo))));
fMwpo = fMwpo./(fMwpo(ind)/max(hJPpo)); 

%% Plotting
figure; 
subplot(121)
scatter(eventtimes, Mws_GRD, '.k') % assuming GRD
hold on
scatter(eventtimes(shutinindex), Mws_JP(shutinindex), 'r') % pre-shut-in
scatter(eventtimes(~shutinindex), Mws_JP(~shutinindex), 'b') % post-shut-in
plot(eventtimes, rs2mw(avec), 'k')  % perturbed region radius
plot(s.tt, log10(Rt), 'g', 'LineWidth', 3)  % rate of earthquakes
plot([T, T], [0, 5])
ylim([0,5])
xlim([0, 30])
ylabel('M_W')

subplot(122)
semilogy(xGRp,hJPp, '.r', 'MarkerSize', 20)
hold on
semilogy(mwsp, max(fMwp)*10.^(-b*mwsp), '--k') 
semilogy(mwsp, fMwp, 'r')

semilogy(xGRpo,hJPpo, '.b','MarkerSize', 20)
semilogy(mwspo, fMwpo, 'b')
semilogy(mwspo, max(fMwpo)*10.^(-b*mwspo), '--k') 
refline(0,1)
ylim([.1, max(fMwp)])
xlim([0, 4.5])
title('Pre-shut-in FMD')
