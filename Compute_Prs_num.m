function [ Prs, rx, Frs,Pin, Pp, Pr] = Compute_Prs_num (b,minrho,maxrho, Nx)
%Compute_Frs_JP: This function computes the cdf for the Maurer & 
% Segall solution to predicting earthquake magnitudes, uses GRD
% distribution for rho. Note that maxrho and minrho should be normalized by
% a. 
   
    % rhos to integrate over
    rhos = logspace(log10(minrho), log10(maxrho), 2*Nx)';
    frho = @(x) computeGRD(x, b, minrho); 
    rhogrd = frho(rhos); 
    
    % rx - r_s values to evaluate
    rx = logspace(log10(minrho), log10(min(1,maxrho)), Nx+1);  
    rx(end) = []; 
    fgrd= frho(rx)'; 
    zc = sqrt(1 - rx.^2); 

    % Pin
    tin = zc - rx.*acos(rx) - (zc.^3./3);
    Pin = tin.*fgrd/2; 
    Pin = Pin(:); 

    % check which terms are needed
    intcond = logical(max(0, maxrho - 1));
       
    switch intcond
        case true
            ft3b = rx./zc;
            [Pp, Pr] = deal(zeros(length(rx),1)); 
            for k = 1:length(rx)
                
                % P_p
                gamma = rx(k); 
                zgam = zc(k); 
                ind = rhos >= gamma;
                geom_term_Pp = acos(gamma) + (2*rhos(ind)-3*gamma).*zgam; 
                Pp(k) = trapz(rhos(ind), geom_term_Pp(:).*rhogrd(ind)); 
    
                % P_r
                geom_term_Pr = (rhos(ind) - gamma).^2; 
                Pr(k) = .5*ft3b(k)*trapz(rhos(ind), geom_term_Pr.*rhogrd(ind));
            end
            
        case false
            %Pp
            [Pp, Pr] = deal(zeros(length(rx),1)); 
            for k = 1:length(rx)-1
                
                % P_p
                gamma = rx(k); 
                zgam = zc(k); 
                ind = rhos > gamma;
                geom_term_Pp = acos(gamma) + (2*rhos(ind)-3*gamma).*zgam; 
                if sum(ind) > 1
                    Pp(k) = trapz(rhos(ind), geom_term_Pp(:).*rhogrd(ind)); 
                else
                    Pp(k) = 0; 
                end
            end
    end
       
    % put everything together
    Prs = Pin + Pp + Pr; 
    
    %% Normalize to integer area
    Prs = Prs./(trapz(rx, Prs)); 
    
    %% Compute the CDF
    Frs = zeros(size(Prs)); 
    for loop = 2:length(Frs)
        Frs(loop) = trapz(rx(1:loop), Prs(1:loop)); 
    end
    
    rx = rx(:); 
end
