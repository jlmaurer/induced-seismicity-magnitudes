function [mws, fMw] = rspdf2mwpdf(rs, frs, dtau, c)
% this function converts from the function on rs to the function on
% magnitude

    mws = rs2mw(rs, dtau, c);  
    ex = .5*(mws(:) + 6.03); 
    numer = log(10)/2; 
    denom = (pi*dtau/c)^(1/3);
    jac = (numer/denom)*10.^(ex); 
    fMw = frs(:).*jac;
end