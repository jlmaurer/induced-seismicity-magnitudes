%% Script to generate plots from Maurer & Segall, Fig. 5
clear

%% Pick parameter values
a = 250;           % 250 or 2500 m
maxrho = 1500;      % 500 or 1500 m

% other parameters
Mwmin = 0.5; 
minrho = mw2rs(Mwmin); 
b = 1;
dtau = 3e6; 
Nx = 1000; 

%% Simulate events numerically
M = 5e6;                    % number of simulations
genrho = @(rhomax,rhomin, b) genRsGR(rhomax,rhomin, b); 
[rsin,rsp1,rsp2, rsr, rho] = SimEq3D(M, genrho, maxrho/a, minrho/a, 1, b); 

% combine events
xlocs = linspace(Mwmin, rs2mw(maxrho, dtau), 100); 
rs = [rsin(:); rsp1(:);rsp2(:); rsr(:)].*a; 
Mws = rs2mw(rs, dtau); 
[hm, xm] = hist(Mws, xlocs); 
ya = trapz(xm, hm); 
hm = hm./ya; 

Mw_pin = rs2mw(rsin(:).*a, dtau); 
[hm_pin, xm_pin] = hist(Mw_pin, xlocs); 
ya = trapz(xm_pin, hm_pin); 
hm_pin = hm_pin./ya; 

%% Compute the analytical expressions
[ Prs, rx, Frs] = Compute_Prs(b, minrho/a,maxrho/a, Nx);
[Pin] = Compute_Pin(rx)./a;
rx = rx.*a; 
Prs = Prs./a; 
[mwx,fMw] = rspdf2mwpdf(rx, Prs, dtau);
Pin_full = Pin.*computeGRD(rx, b, minrho);
[~,fPin] = rspdf2mwpdf(rx, Pin_full, dtau);
ya = trapz(mwx, fPin); 
fPin = fPin./ya; 

%% Plot the results
mwa = rs2mw(a, dtau); 

figure; 
semilogy(mwx, fPin./max(fPin), 'g')
hold on
semilogy(mwx, fMw./max(fMw), ':r', 'LineWidth', 2)
semilogy(xm, hm./max(hm), 'or')
semilogy(xm_pin, hm_pin./max(hm_pin), 'dr')
vline(mwa);
xlabel('M_W')
ylabel('Normalized Frequency')
legend('P_{in} - Theoretical', 'P(r_s) - Theoretical', 'P_{in} - Numerical', ...
    'P(r_s) - Numerical', 'M_W(a(t))', 'Location', 'Southwest')
title('Components of P(r_s)')
axis tight
