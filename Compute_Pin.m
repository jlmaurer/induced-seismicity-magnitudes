function [Pin] = Compute_Pin(rx, a)
% This function computes Pin, as presented by Maurer & Segall (2018), and
% originally formulated in Segall & Lu (2015). Does NOT include the GRD
% term, so this term must be multiplied by GRD to get the full
% distribution. 

    if nargin < 2
        rxda = rx; 
    else
        rxda = rx./a; 
    end
    
    zc = sqrt(1 - rxda.^2); 
    racr = rxda.*acos(rxda); 
    
    Pin = 0.5*(zc - racr - ((zc.^3)/3)); 
end
