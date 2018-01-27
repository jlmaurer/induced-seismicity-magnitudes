%% Script to generate figures from Maurer & Segall, Fig. 8

% parameters
b = 1; 
Mmax = 5; 
Mc = 0; 
dtau = 3e6;         % stress drop
c = 7*pi/16;        % geometrical constant for stress drop
kf = .1;           % diffusivity
Neq = 1e3;          % number of earthquakes at each time step
Nx = 3000; 
Nt = 2000; 

% time in seconds at r = 100 m
tfinal = (100/2)^2 / kf; 
time_vec = linspace(3600, tfinal, Nt); 

% compute a(t) for the range of times needed
at = 2*sqrt(kf*time_vec); 
rmin = mw2rs(Mc, dtau, c)./at; 
rmax = mw2rs(Mmax, dtau, c)./at;

% Compute P(r_s) for each time
[Rx, Rxu, Prs, Pin, Pshap_lo, Pshap_hi] = deal(zeros(Nx, length(time_vec))); 
for loop = 1:length(time_vec)
    [prs, rx]  = Compute_Prs(b, rmin(loop), rmax(loop), Nx);
    pgrd = computeGRD(rx, b, rmin(loop)); 
    pin = Compute_Pin(rx); 
    [pshaplo] = Compute_Shapiro_lower(rx);  
    [pshaphi, rxu] = Compute_Shapiro_upper(rmin(loop), rmax(loop), length(rx)); 
    Rx(:,loop) = rx.*at(loop); 
    Rxu(:,loop) = rxu.*at(loop); 
    Prs(:,loop) = prs./at(loop); 
    Pin(:,loop) = pin.*pgrd./at(loop);
    Pshap_lo(:,loop) = pshaplo.*pgrd./at(loop); 
    Pshap_hi(:,loop) = pshaphi.*computeGRD(rxu, b, rmin(loop))./at(loop); 
end

% integrate over time
[Frs,RX] = compute_timeintegrated_prs(Rx,Prs,time_vec); 
[Fin] = compute_timeintegrated_prs(Rx,Pin,time_vec); 
[Fshap_lo] = compute_timeintegrated_prs(Rx,Pshap_lo,time_vec); 
[Fshap_hi, RXU] = compute_timeintegrated_prs(Rxu,Pshap_hi,time_vec); 

% convert the pdfs to magnitude
[mx,fMx] = rspdf2mwpdf(RX, Frs, dtau, c);
[~,fin] = rspdf2mwpdf(RX, Fin, dtau, c);
[~,fshap_lo] = rspdf2mwpdf(RX, Fshap_lo, dtau, c);
[mxu,fshap_hi] = rspdf2mwpdf(RXU, Fshap_hi, dtau, c);

% remove first ten elements, these do not get contributions from each time
% step. 
fshap_hi(1:10) = []; 
fshap_lo(1:10) = []; 
fin(1:10) = []; 
fMx(1:10) = []; 
mx(1:10) = [];
mxu(1:10) = []; 
mxu = mxu(:); 

% Ordinary GRD = Shapiro intermediate model (time-invariant)
y = 10.^(-b*mx); 
yu = 10.^(-b*mxu); 

% Scale distributions by first element
fin = fin./max(fin);
y = y(:)./max(y);
fshap_lo = fshap_lo./max(fshap_lo); 
fshap_hi = fshap_hi./max(fshap_hi);
fMx = fMx./max(fMx); 

% Plot the results

figure; 
semilogy(mx, fin)
hold on
semilogy(mx, fshap_lo)
semilogy(mx, y, '--k')
semilogy(mx, fMx)
semilogy(mxu, fshap_hi)
ylim([1e-4, 1])
xlim([Mc, max(mxu)])
legend('P_{in}', 'Shapiro-lower', 'GRD', 'This study', 'Shapiro-upper')
xlabel('Magnitude')
ylabel('Relative Frequency')
title(['Time-integrated FMDs for b = ', num2str(b)])

figure; 
plot(mx, fin./y)
hold on
plot(mx, fshap_lo./y)
plot(mx, y./y, '--k')
plot(mx, fMx./y)
plot(mxu, fshap_hi./yu)
% ylim([1e-4, 1])
xlim([Mc, max(mxu)])
xlabel('Magnitude')
legend('P_{in}', 'Shapiro-lower', 'GRD', 'This study', 'Shapiro-upper')
ylabel('Frequency relative to GRD')
title(['Time-integrated FMDs for b = ', num2str(b)])
set(gca, 'yscale', 'log')
