function [mw] = mo2mw(Mo, unit)
% This function converts moment Mo to moment magnitude. Default unit is N m

if nargin<2, unit='Nm'; end

if strcmp(unit, 'Nm')==1
    mw = (2/3)*log10(Mo*1e7) - 10.7; 
elseif strcmp(unit, 'dc')==1
    mw = (2/3)*log10(Mo) - 10.7; 
else
    error('Unit must be either Nm or dc')
end

end