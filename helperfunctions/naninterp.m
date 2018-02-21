function X = naninterp(X)
% Interpolate over NaNs
% See INTERP1 for more info
% X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)),'pchip'); % original
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)),'pchip','extrap'); % edited by TOB 20171215
return