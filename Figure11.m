%% Figure 11 a & c in Maurer & Segall
% Use data from Kraft, T. & Deichmann, N., (2014).
% High-precision relocation and focal mechanism of the injection 
% induced seismicity at the Basel EGS, Geothermics, 52, 59-73.
clear

%% FOR BASEL DATA: 
% use 2D solution, a(t) = r(t) = 100 m/day. 
data = readtable('Basel/Basel_Data.csv'); 
mws = data.Mw; 
t = datenum(data.DateTime); 
tinj = t - t(1); 
% %cumulative number of events
% Ncum = cumsum([1:length(t)]); 
rho = sqrt(data.e.^2 + data.n.^2 + data.z.^2); 

% shut-in date
injStopDate = datenum(2006,12,8,11,33,0); % this is the actual shut-in time
index = t <= injStopDate;
nindex = logical(1-index);

%% Seismicity R-T plot
% estimated diffusivity m^2/day (fit by eye) = 0.2025 m/sec
Df = 1.75e4;

Nd2show = 10; 
figure; 
scatter(tinj(index), rho(index), 10.^mws(index), mws(index), 'filled')
hold on
scatter(tinj(nindex), rho(nindex), 10.^mws(nindex), mws(nindex))
xlabel(['First ',num2str(Nd2show), ' days of events'])
ylabel('Distance (m)')
xlim([0, Nd2show])
colorbar
caxis([0.5, 3])
plot(tinj, 2*sqrt(Df*tinj), 'LineWidth', 1)

%% Compute Mw^max from the estimated volume of the perturbed region
dtau = 3e6;                 % assume 3 MPa stress drop
r_est = 2*sqrt(Df*tinj);    % use radius estimated from seismicity cloud 
Mmax_est = rs2mw(r_est, dtau); 

figure; 
plot(tinj, Mmax_est, 'LineWidth', 1)
hold on
scatter(tinj(nindex), mws(nindex), 'filled')
scatter(tinj(index), mws(index), 'filled')
xlim([0, Nd2show])
xlabel(['First ',num2str(Nd2show), ' days of events'])
ylabel('M_W')
legend('Estimated Mmax', 'Pre shut-inM_W', 'Post shut-in M_W')
