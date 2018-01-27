function [Pin] = Compute_Pin(rx, b, rmin, a)
% This function computes Pin, as presented by Maurer & Segall (2018), and
% originally formulated in Segall & Lu (2015). 

    rxda = rx./a; 
    Pgrd = computeGRD(rxda, b, rmin); 
    zc = sqrt(1 - rxda.^2); 
    racr = rxda.*acos(rxda); 
    
    Pin = 0.5*(zc - racr - ((zc.^3)/3)).*Pgrd; 
end
