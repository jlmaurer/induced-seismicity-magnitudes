%% Script to replicate Fig. 7 in Maurer & Segall, 2018
clear

%% Parameters

% diffusivity
c = 0.01; 

% seismicity parameters
Mwmax = 5; 
Mwmin = 0.5; 
b = 1; 
dtau = 3e6;

%% Compute a(t)
tdays = [1, 4:4:150]; 
nsec = 24*3600; 
ts = tdays*nsec; 
at = sqrt(c*ts); 

%% Compute P(r_s) for each a(t)
Nx = 300; 
Nt = length(ts); 
[mx, Rx, fM1, fM2, fM3] = deal(zeros(Nx, Nt)); 
for loop = 1:Nt
    a = at(loop); 
    rhomin = mw2rs(Mwmin, dtau)./a; 
    rhomax = mw2rs(Mwmax, dtau)./a; 
%     all_rhomin(loop) = rhomin; 
%     all_rhomax(loop) = rhomax; 
    
    % compute P(r_s)
    [Prs, rx, ~, Pin, Pp] = Compute_Prs_num (b,rhomin,rhomax, Nx); 
    Rx(:,loop) = rx.*a; 
    P1 = Pin./a; 
    P2 = (Pin + Pp)./a;
    P3 = Prs./a; 
    [mx(:,loop),fm2] = rspdf2mwpdf(Rx(:,loop), P2, dtau);
    [~,fm3] = rspdf2mwpdf(Rx(:,loop), P3, dtau);
    [~, fm1] = rspdf2mwpdf(Rx(:,loop), P1, dtau);
    if fm1(1) < fm1(2)
        fM1(:,loop) = [fm1(2:end)./fm1(2), 0];     
    else
        fM1(:,loop) = fm1./fm1(1); 
    end
    
    if fm2(1) < fm2(2)
        fM2(:,loop) = [fm2(2:end)./fm2(2), 0];        
    else
        fM2(:,loop) = fm2./fm2(1);
    end
    
    if fm3(1) < fm3(2)        
        fM3(:,loop) = [fm3(2:end)./fm3(2), 0]; 
    else
        fM3(:,loop) = fm3./fm3(1); 
    end
end

%% Plot results
MX = mx(:,end); 
cs = varycolor(length(tdays)); 

figure; 
for k = 1:size(fM1,2)
    semilogy(mx(:,k), fM1(:,k), 'color', cs(k,:))
    hold on
end
axis tight
colorbar;
y = 10.^(-b*MX); 
y= y(:)./max(y);
ylim([1e-4, 1])
semilogy(MX, y, 'k--')
xlabel('M_W')
ylabel('Relative Frequency')
title('P_{in}')

figure; 
for k = 1:size(fM2,2)
    semilogy(MX, fM2(:,k), 'color', cs(k,:))
    hold on
end
axis tight
colorbar;
y = 10.^(-b*mx(:,end)); 
y= y(:)./max(y);
ylim([1e-4, 1])
semilogy(MX, y, 'k--')
xlabel('M_W')
ylabel('Relative Frequency')
title('P_{p}')

figure; 
for k = 1:size(fM3,2)
    semilogy(MX, fM3(:,k), 'color', cs(k,:))
    hold on
end
axis tight
colorbar;
y = 10.^(-b*MX); 
y= y(:)./max(y);
ylim([1e-4, 1])
semilogy(MX, y, 'k--')
xlabel('M_W')
ylabel('Relative Frequency')
title('P(r_s)')
