function [b_TM] = bval_TM(mags, dM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function take an earthquake catalog in a specific format and
% calculates a b-value using the method of Tinti and Mulargia (1987), as 
% related by Marzocchi and Sandri (2003). 
%
% Input:
% mag       = Earthquake catalog with magnitudes only 
% dM        = magnitude bin used in the catalog, e.g. 0.1. 
%
% Output:
% bvalue    = b-value estimate from Tinti & Mulargia (1987)
%
% Written by Jeremy Maurer
%
% Last modified on  Wed March 1 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % estimate the magnitude bin if not supplied
    if nargin < 2
        diffs = diff(sort(mags)); 
        diffs(diffs==0) = []; 
        dM = min(diffs); 
    end

    if sum(isnan(mags))~=0
        warning('NANs are present')
    end

    % mean and minimum magnitude
    mean_mag = nanmean(mags);
    min_mag = nanmin(mags);

    p = 1 + (dM/(mean_mag - min_mag)); 
    const = 1/log(10); 

    b_TM = (const/dM)*log(p); 

end
