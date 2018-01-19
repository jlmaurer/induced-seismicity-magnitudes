function [Prs, rx,Frs] = Compute_Prs (b, minrho,maxrho, Nx)
%Compute_Frs_an: This function computes the cdf and pdf for the Maurer & 
% Segall solution to predicting FMDs. Note that all dimensions (maxrho, 
% minrho) should be normalized by a(t). 

    % normalized range of event sizes
    rx = logspace(log10(minrho), log10(min(1,maxrho)), Nx+1);  
    rx(end) = []; 

    % scaled parameters
    rdrm = rx./minrho; 
    Rm = maxrho/minrho;
    alpha = -2*b; 

    % use the appropriate expressions for b=1, b>1, or numerical approach
    % for b < 1
    if b == 1
        F1 = FB(1, rx, maxrho); 
        H1 = HB(1, rx, maxrho); 
        gam1 = Gam(rx); 
        Prs = (rdrm.^alpha).*(gam1 + F1 + H1); 
    elseif b > 1
        fb = FB(b, rx, maxrho);
        hb = HB(b, rx, maxrho);
        Prs = (rdrm.^alpha).*fb + (Rm.^alpha).*hb; 
    else
        Nx = 1000; 
        [ Prs, rx] = Compute_Prs_num (maxrho,minrho, b,Nx);
    end

    
    %% Normalize to integer area
    % this is necessary because a number of sources end up located fully
    % outside the perturbed region
    Prs = Prs./(trapz(rx, Prs)); 
    
    %% Compute the CDF
    Frs = zeros(size(Prs)); 
    for loop = 2:length(Frs)
        Frs(loop) = trapz(rx(1:loop), Prs(1:loop)); 
    end
    
    rx = rx(:); Prs = Prs(:); Frs = Frs(:); 
    
end

function [fb] = FB(b, x, max_rho)
% this function computes Fb
    if b==1 && max_rho == inf
        error('Cannot have b = 1, max_rho = inf, please supply a finite max_rho')
    end
    [~, rz, zdr] = ZC(x);
    alpha = -2*b; 
    
    if b == 1
        R = x./max_rho; 
        fb = (acos(x) - 3*rz).*(1-(R.^2)) + 4*rz.*(1-R);
    elseif b > 1
        t1 = (b/3)*zdr.*(2 + x.^2);
        t2 = (1-b)*acos(x); 
        t3 = (-(alpha + 3)/(1 + alpha))*rz;
        t4 = ((x.^2)./(2*zdr))*(1./((1-b)*(1-2*b))); 
        fb = t1 + t2 + t3 + t4;
    else
        error('Warning: cannot handle b < 1 currently')
    end
end

function [hb] = HB(b, x, max_rho)
% this function computes Hb(x)

    % default case is max_rho = inf, hb = 0
    if nargin < 3
        hb = 0; 
        return
    end
    [zc, zr, zdr] = ZC(x);

    if b==1 
        R = x./max_rho;
        ff = (x.^2)./zdr;
        hb = ff.*(-log(R) - .5*(R.^2) + 2*R - 1.5); 
    elseif b > 1
        t1 = -acos(x); 

        tt2 = (b/(1-b))*(max_rho.^2); 
        
        tt3a = (2-3*(x.^2))./x;
        tt3b = (4*b./(1-2*b)); 
        tt3 = tt3a.*tt3b*max_rho;
        
        tt4 = (6 - 7*(x.^2)); 
        
        hb = t1 + (tt2 + tt3 + tt4)./(2*zdr); 
    else
        error('Warning: cannot handle b < 1 currently')
    end
end

function [gamp] = Gam(x)
% this function computes gamma'_b
    [~, ~, zdr] = ZC(x);
    gamp = zdr.*(2 + x.^2)/3 - acos(x); 
end

function [zc, rz, zdr] = ZC(x)
    if nanmin(x) > 1
        error('rx must be normalized by a(t)')
    end
    zc= sqrt(1 - x.^2); 
    rz = x.*zc; 
    zdr = zc./x;
end

