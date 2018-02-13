function [rs] = genRsGR(rmax,rmin, b)
% this function takes inputs that have been normalized by rreal.
    
    u = rand; 
    rs = rmin*(u^(-1/(2*b))); 
    
    while rs > rmax
        u = rand; 
        rs = rmin*(u^(-1/(2*b))); 
    end
end
