%% Script to generate figures from Maurer & Segall, Fig. 8

% parameters
b = 1; 
Mmax = 5; 
Mc = 0; 
dtau = 3e6;         % stress drop
c = 7*pi/16;        % geometrical constant for stress drop
kf = .1;           % diffusivity
Neq = 1e3;          % number of earthquakes at each time step
Nx = 500; 
Nt = 2000; 

% time in seconds at r = 100 m
tfinal = (100/2)^2 / kf; 
time_vec = linspace(3600, tfinal, Nt); 

% compute a(t) for the range of times needed
at = 2*sqrt(kf*time_vec); 
rmin = mw2rs(Mc, dtau, c)./at; 
rmax = mw2rs(Mmax, dtau, c)./at; 

% Compute P(r_s) for each time
[Rx, Prs, Pin, Pshap_lo, Pshap_hi] = deal(zeros(Nx, sum(tindex))); 
for loop = 1:length(time_vec)
    [prs, rx]  = Compute_Prs(b, rmin(loop), rmax(loop), Nx);
    pgrd = computeGRD(rx, b, rmin(loop)); 
    pin = Compute_Pin(rx).*pgrd; 
    [pshaplo] = Compute_Shapiro_lower(rx);  
    Rx(:,loop) = rx.*at(loop); 
    Prs(:,loop) = prs./at(loop); 
    Pin(:,loop) = pin./at(loop);
    Pshap_lo(:,loop) = pshaplo.*pgrd./at(loop); 
end

% integrate over time
[Frs,RX] = compute_timeintegrated_prs(Rx,Prs,times); 
[Fin] = compute_timeintegrated_prs(Rx,Pin,times); 
[Fshap_lo] = compute_timeintegrated_prs(Rx,Pshap_lo,times); 

% convert the pdfs to magnitude
[mx,fMx] = rspdf2mwpdf(RX, Frs, dtau, c);
[~,fin] = rspdf2mwpdf(RX, Fin, dtau, c);
[~,fshap_lo] = rspdf2mwpdf(RX, Fshap_lo, dtau, c);

% Shapiro upper end-member model
rmin_full = mw2rs(Mc, dtau, c); 
rmax_full = mw2rs(Mmax, dtau, c); 
pshaphi = Compute_Shapiro_upper(rmin_full, rmax_full, Nx); 
Pshap_hi = pshaphi.*computeGRD(mw2rs(mx, dtau, c), b, rmin_full); 
[~,fshap_hi] = rspdf2mwpdf(RX, Pshap_hi, dtau, c);

% GRD
y = 10.^(-b*mx); 
y= y(:)./max(y);


figure; 
semilogy(mx, fin./max(fin))
hold on
semilogy(mx, fshap_lo./max(fshap_lo))
semilogy(mx, y, '--k')
semilogy(mx, fMx./max(fMx))
semilogy(mx, fshap_hi./max(fshap_hi))
ylim([1e-4, 1])
xlim([Mc, max(mx)])
xlabel('Magnitude')
ylabel('Relative Frequency')
title(['Time-integrated FMDs for b = ', num2str(b)])
