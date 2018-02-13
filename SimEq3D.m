function [rsin,rsp1,rsp2, rsr, rho, distrs,distrho, rc, r, z,Zcount] = SimEq3D (M, genrho, maxrho, minrho,r_sphere, varargin)
% This function runs a simulation to generate r_s-sized events for a given
% r-sized perturbed zone, source distribution function frho, and maximum
% event size rhomax. 
rng('shuffle')

[z,rsin,rsp1,rsp2, r, rc, rho,rsr, distrho, distrs]= deal(zeros(M, 1)); 
Zcount = 0; 

genz = @(r_sphere) 2*r_sphere*rand - r_sphere; 
genrc = @(maxrho, r_sphere) (maxrho + r_sphere)*sqrt(rand); 

% loop over values of rs
for k = 1:M
    z(k) = genz(r_sphere); 
    r(k) = sqrt(r_sphere^2 - z(k)^2); 

    rho(k) = feval(genrho, maxrho,minrho, varargin{:}); 
    rc(k) = genrc(maxrho, r_sphere); 
    
    distrho(k) = sqrt(rc(k).^2 + z(k).^2); 
    
    if rho(k) <= r(k)
        if rc(k) < r(k) - rho(k)
            rsin(k) = rho(k); 
            distrs(k) = distrho(k); 
        elseif rc(k) > r(k)+rho(k)
            Zcount = Zcount +1; 
            distrs(k) = nan; 
        else
            rsp1(k) =(r(k) + rho(k) - rc(k))/2;
            distrs(k) = r(k) - rsp1(k); 
        end
        
    elseif rho(k) > r(k)
        if rc(k) < rho(k) - r(k)
            rsr(k) = r(k); 
            distrs(k) = 0; 
        elseif rc(k) > r(k)+rho(k)
            Zcount = Zcount +1; 
            distrs(k) = nan; 
        else
            rsp2(k) =(r(k) + rho(k) - rc(k))/2;
            distrs(k) = r(k) - rsp2(k); 
        end    
    else
        error('Somethings wrong')
    end
end

% % delete parallel pool
% poolobj = gcp('nocreate'); % If no pool, do not create new one.
% delete(poolobj);

% remove non-events from rc, z to compute 
R = rsin + rsp1 + rsp2 + rsr; 
indR = R < minrho; 
nanR = isnan(R); 
rc(indR) = []; 
rc(nanR) = []; 
z(indR) = []; 
z(nanR) = []; 
distrho(indR) = []; 
distrho(nanR) = []; 
r(indR) = []; 
r(nanR) = []; 
distrs(indR) = []; 
distrs(nanR) = []; 

% remove non-events from rs's
rsin(isnan(rsin)) = []; 
rsin(rsin < minrho) = []; 
rsp1(isnan(rsp1)) = []; 
rsp1(rsp1 < minrho) = []; 
rsp2(isnan(rsp2)) = []; 
rsp2(rsp2 < minrho) = []; 
rsr(rsr < minrho) = []; 

end