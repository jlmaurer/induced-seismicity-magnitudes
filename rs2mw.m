function mw = rs2mw (rs, deltatau, c)
% Function converts a source radius to an equivalent magnitude given a
% stress drop. 
    if nargin < 3
        c = 7*pi/16; 
    end
    if nargin < 2
        deltatau = 3e6; 
    end

    t1 = deltatau*(pi/c)*(rs.^3); 
    mw = (2/3)*log10(t1) - 6.03; 
end