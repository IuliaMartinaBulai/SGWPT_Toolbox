%--------------------------------------------------------------------------
% File: SNR.m
%
% Goal: compute the SNR between f (original signal) and g (noisy signal)
%
% Use: [value, sigma] = SNR(g, f)
%
% Inputs:  g - restored or noisy signal
%          f - original reference signal
%
% Outputs: value - signal/Noise ratio
%
% Recalls: norm.m
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
function [value] = SNR(g, f) 
value = 20*log10(norm(f,'fro')/norm(g-f,'fro'));
end