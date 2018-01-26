%% Make figures using the Basel data
% This script generates plot using seismicity data from Kraft and 
% Deichmann, 2014, shown in Maurer and Segall (2018). 

%% Load Basel Data 
data = readtable('Basel_Data.csv'); 
mws = data.Mw; 
t = datenum(data.DateTime); 
tdays = t - t(1); 
nsec = 24*3600; 
rho = sqrt(data.e.^2 + data.n.^2 + data.z.^2); 

% shut-in time
injStopDate = datenum(2006,12,8,11,33,0) - t(1); 

% use Mc = 0.9 (Bachmann et al., 2011)
Mc = 0.9; 
m = mws(mws >=Mc); 
t_mc = tdays(mws>=Mc);

%% R-T plot
% we use the seismicity cloud to estimate an approximate stimulated radius
tt = linspace(0, 10*nsec, 1000); 

% use diffusivity Df = 0.202
Df = 0.202; 
r_seis = 2*sqrt(Df*tt); 

% plot first 10 days
day10ind = tdays <=10; 
t10 = tdays(day10ind); 
rho10 = rho(day10ind); 
m10 = mws(day10ind);
sind = t10 > injStopDate; 

figure; 
scatter(t10(~sind), rho10(~sind), 10.^m10(~sind), m10(~sind), 'filled')
hold on
scatter(t10(sind), rho10(sind), 10.^m10(sind), m10(sind))
plot(tt./nsec, r_seis)
xlabel('First 10 days of events')
ylabel('Distance (m)')
colorbar
caxis([0.5, 3])

%% Mmax with Time
shutinInd = t > injStopDate; 

%% Assuming r_seis

% assume 2 MPa stress drop
dtau = 2e6; 

% compute Mmax assuming r_seis
Mmax_seis = rs2mw(r_seis, dtau); 

% plot seismicity and Mmax prediction
figure; 
scatter(t10(~sind), m10(~sind), 'filled')
hold on
scatter(t10(sind), m10(sind), 'filled')
plot(tt./nsec, Mmax_seis, 'LineWidth', 1)
xlabel('First 10 days of events')
ylabel('M_W')
% xlim([0, 10])

%% Assuming an approximation to injected volume
% Use the McGarr and Barbour (2017) model to estimate Mmax with time

% first estimate the approximate rate of injection from Bachmann et al.
% 2011 and integrate to get the total injected volume with time. Can check
% against the total (overall) injected volume of 11,570 m^3. 
fr_init = 0;
fr_final = 3200*(1/60)*(1/1000); 
fr_accel = (fr_final - fr_init)/(5*nsec); 
[total_vol_m3 ,rate_m3 ] = deal(zeros(size(tt))); 
rate_m3(1) = fr_accel*(tt(2)-tt(1));
for k = 2:length(tt)
    rate_m3(k) = trapz(tt(1:k), fr_accel*ones(1,k));
end
total_vol_m3(1) = rate_m3(1)*(tt(2)-tt(1)); 
for k = 2:length(tt)
    total_vol_m3(k) = trapz(tt(1:k), rate_m3(1:k));
end
obs_total = 11570; 

% compute the Mcgarr, 2017 Mmax estimate
mu = 30e9; 
M0_mcgarr = 2*mu*total_vol_m3; 
Mmax_mcgarr = mo2mw(M0_mcgarr); 

plot(tt./nsec, Mmax_mcgarr, 'LineWidth', 1)

%% Expected assuming GRD and the cumulative number of events

% use the pre-shut-in b-value 
b_preshutin = 1.56;
t_mc_10 = t_mc(t_mc <= 10);

% cumulative number of events
Ncum = cumsum(m10(m10 >= Mc)); 

% compute Mmax using van der Elst et al. (2017)
[Mmax_GR, Mlo, Mhi] = compute_Mmax_GR(Mc, b_preshutin, Ncum, 95); 

plot(t_mc_10, Mmax_GR, 'LineWidth', 1)
plot(t_mc_10, Mlo, ':k', t_mc_10, Mhi, ':k')
legend('Pre shut-in', 'Post shut-in', 'Mmax\_seis (\Delta\tau = 2 MPa)', ...
    'Mmax\_McGarr', 'Mmax\_GRD', 'Mmax\_GRD\_95%\_bounds')

