function [Prs_total, RX] = compute_timeintegrated_prs (rx,Prs,times)
% this function computes the time-integrated probability of an r_s-sized
% event.

    % rx, Prs - matrices of Nx x Nt
    % times - vector Nt x 1 of times

    minrx = min(rx(:)); 
    maxrx = max(rx(:)); 
    RX = logspace(log10(minrx), log10(maxrx), 5000); 
    p = zeros(length(RX), size(rx,2)); 
    for loop = 1:size(rx,2)
        try
            p(:,loop) = interp1(rx(:,loop), Prs(:,loop), RX); 
        catch
            keyboard; 
        end
        
    end
    p(isnan(p)) = 0; 
    lr = length(RX); 
    Prs_total = zeros(lr,1); 
    for loop2 = 1:lr
        Prs_total(loop2) = trapz(times,p(loop2,:));
%         Prs_total(loop2) = sum(p(loop2,:));
%         if Prs_total(loop2) < 0
%             keyboard; 
%         end
    end
    ta = trapz(RX, Prs_total); 
    Prs_total = Prs_total./ta; 
end