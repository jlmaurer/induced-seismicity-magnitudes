function [frs, Frs] = computeGRD(rstest, b, rmin)
% computeGRD: This function computes the probability of an r_s-sized source
% based on the GRD with b given.

    ind = rstest < rmin; 
    
    alpha = (2*b);   
    lf = log(alpha) + alpha*log(rmin) - (alpha+1)*log(rstest); 
    frs = exp(lf); 
    frs(ind) = 0; 
    Frs = 1 - ((rstest./rmin).^(-alpha)); 
end