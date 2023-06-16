%--------------------------------------------------------------------------
% File: thresholding_coef.m
%
% Goal: Given a threshold value tau and a vector, find the thresholding
% coefficient and the tresholded vector
%
% Use: [best_basis_thresholded] = thresholding_coef(tau,best_basis)
%
% Inputs: tau - threshold value
%         best_basis - input vector         
%
% Outputs: best_basis_thresholded - thresholded vector
%
% Recalls: abs.m, sign.m, tmp.m
%
% Authors: IM Bulai, S Saliani
% Date last modified: September, 2022
%
% This file is part of the SGWPT toolbox
% (Spectral Graph Wavelet Packet Transform Toolbox)
% Copyright (C) 2022, Iulia Martina Bulai and Sandra Saliani.
%
% The SGWPT toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The SGWPT toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the SGWPT toolbox.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------
function [best_basis_thresholded] = thresholding_coef(tau,best_basis)
for i=2:numel(best_basis) % the scale coefficients are not thresholded 
    % (in order to be able to compare with Irion's results)
    best_basis_thresholded{1}  = best_basis{1};
    tmp{i} = abs(best_basis{i})-tau;
    best_basis_thresholded{i} = sign(best_basis{i}).*tmp{i}.*(tmp{i}>0);
end
end