%% Script to generate figures from Maurer & Segall, Fig. 5

% pick a(t), rho_max
a = 2500;           % 250 or 2500 m
maxrho = 1500;      % 500 or 1500 m

% other parameters
minrho = mw2rs(0.5); 
b = 1;
Nx = 1000; 

% Compute the analytical expression
[ Prs, rx, Frs] = Compute_Prs(b, minrho/a,maxrho/a, Nx);
[Pin] = Compute_Pin(rx, a);

% Plot the individual contribution
figure; 
loglog(rx.*a, Pin./a, 'b')
hold on
loglog(rx.*a, Prs./a, ':', 'LineWidth', 2)
plot([a, a], [min([Pp(:); Pin(:)]./a), max([Pp(:); Pin(:)]./a)])
xlabel('Normalized Fault Size')
ylabel('Relative Frequency')
legend('P_{in}', 'P_p', 'P_r', 'a(t)')
title('Components of P(r_s)')
axis tight

% Plot the full solution and compare to Segall & Lu, 2015
figure; 
loglog(rx.*a, Prs./a)
hold on
loglog(rx.*a, Pin./a, '--')
plot([a, a], [min([Prs(:); Pin(:)]./a), max([Prs(:); Pin(:)]./a)])
xlabel('M_W')
ylabel('Relative Frequency')
legend('Maurer-Segall 2017', 'Segall-Lu 2015', 'a(t)')
title('Comparing Maurer & Segall 2017 to Segall & Lu 2015')
axis tight
