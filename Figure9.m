%% Script to generate figures from Maurer & Segall, Fig. 8

% specify b-value as 1 or 1.4 to replicate Fig 8a or b, respectively
b = 1; 

% other parameters
Mmax = 5; 
Mc = 0; 
dtau = 3e6;         % stress drop
c = 7*pi/16;        % geometrical constant for stress drop
kf = .01;           % diffusivity
Neq = 1e3;          % number of earthquakes at each time step

% discretize in time, Mw
Nx = 500; 
Nt = 2000; 

% compute a(t) for the range of times needed
% days2plot = [2:2:10, 20, 30, 50];
% use r = 100 m
tfinal = (100/2)^2 / kf; 
nsec = 24*3600;     % seconds in a day
time_vec = linspace(1, tfinal*nsec, Nt); 
at = 2*sqrt(kf*time_vec); 
rmin = mw2rs(Mc, dtau, c)./at; 
rmax = mw2rs(Mmax, dtau, c)./at; 

for oloop = 1:length(days2plot)
    tplot = days2plot(oloop)*nsec;
    tindex = time_vec <= tplot;
    times = time_vec(tindex); 
    
    % Compute P(r_s) for each time
    [Rx, Prs, Pin] = deal(zeros(Nx, sum(tindex))); 
    for loop = 1:sum(tindex)
        [prs, rx]  = Compute_Prs(b, rmin(loop), rmax(loop), Nx);
        pin = Compute_Pin(rx, b, rmin(loop), 1); 
        Rx(:,loop) = rx.*at(loop); 
        Prs(:,loop) = prs./at(loop); 
        Pin(:,loop) = pin./at(loop); 
    end

    % integrate over time
    [Prs_tot,RX] = compute_timeintegrated_prs(Rx,Prs,times); 

    % convert the pdf to magnitude
    [mx,fMx] = rspdf2mwpdf(RX, Prs_tot, dtau, c);

    % normalize and plot
    Prs_norm = Prs./repmat(max(Prs), size(Prs, 1), 1); 
    Rxm = rs2mw(Rx); 
    cs = varycolor(length(days2plot)); 
    fMx_norm = fMx./max(fMx); 
    fMx_norm(1) = []; 
    fMx_norm = [fMx_norm; 0]; 
    semilogy(mx, fMx_norm, 'color', cs(oloop,:))
    hold on
    ylim([1e-4, 1])
    xlim([0.5, ceil(max(mx))])
end
y = 10.^(-b*mx); 
y= y(:)./max(y);
semilogy(mx, y, '--k')
xlabel('Magnitude')
ylabel('Relative Frequency')
title(['Time-integrated FMDs for b = ', num2str(b)])
