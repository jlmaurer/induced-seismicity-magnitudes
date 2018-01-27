function [Prs_shapiro_hi] = Compute_Shapiro_upper(rho_min, rho_max, Nx)
% this function returns the probability corresponding to the upper
% end-member model of Shapiro et al. (2011, 2013, Shapiro 2015). 

    if nargin < 3
        Nx = 200; 
    end
    
    rx = linspace(rho_min, rho_max, Nx); 
    % third, P_s/P_c = 1/P_c
    Psu = 1 + 3*(rx.^2)/2 + (3*pi*rx/4); 


end