function rs = mw2rs (mw, deltatau, c)
% Function converts a magnitude to an equivalent source radius given a
% stress drop. 
    if nargin < 3
        c = 7*pi/16; 
    end
    if nargin < 2
        deltatau = 3e6; 
    end

    t1 = (mw +6.03)*(3/2); 
    t2 = 10.^t1; 
    t3 = t2./(deltatau*(pi/c)); 
    rs = t3.^(1/3); 
end