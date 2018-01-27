function [Psl] = Compute_Shapiro_lower (rx)
% This function returns the geometrical modifying factors from Shapiro et
% al. (2011, 2013, 2015), by default assuming a uniform distribution on
% fault size rho. Can be multiplied by GRD to get the corresponding
% distribution assuming GRD. 

    if max(rx) <=1
        zc = sqrt(1 - (rx.^2)); 
    else
        error('Need to normalize r_s by a(t)')
    end

    % first: P_s/P_c = P_vol
    asr = asin(rx); 
    t1 = 1 + .5*(rx.^2); 
    t2 = zc; 
    t3 = 3*pi*rx/4; 
    t4 = 1.5*rx.*asr; 
    Psl = t1.*t2 - t3 + t4; 

end
