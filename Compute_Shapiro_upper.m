function [Prs_shapiro_hi, rx] = Compute_Shapiro_upper(rho_min, rho_max, Nx)
% this function returns the probability corresponding to the upper
% end-member model of Shapiro et al. (2011, 2013, Shapiro 2015). It is
% required for this fcn that rho_min and rho_max are normalized by a(t). 

    if nargin < 3
        Nx = 200; 
    end
    
    if numel(Nx)==1
        rx = linspace(rho_min, rho_max, Nx); 
    else
        rx = Nx; 
    end
    
    % third, P_s/P_c = 1/P_c
    Prs_shapiro_hi = 1 + 3*(rx.^2)/2 + (3*pi*rx/4); 
    Prs_shapiro_hi = Prs_shapiro_hi(:); 

end