% sgwpt_setpath : Set paths for SGWPT toolbox
%
% user should set SGWPT_ROOT to be directory where sgwpt_toolbox is installed
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

SGWPT_ROOT=[fileparts(mfilename('fullpath')),filesep];

sgwpt_relpathlist={'',...
                  'utils',...
                  'utils\spm12',...
                  'utils\cbiNifti',...
                  'data',...
                  'examples',...
                  'examples\Toronto',...
                  'examples\Brain',...
                 };
fprintf('Welcome to sgwpt_toolbox. SGWPT root directory is %s\n',SGWPT_ROOT);
for k=1:numel(sgwpt_relpathlist)
  sgwpt_tmp_pathname=[SGWPT_ROOT,sgwpt_relpathlist{k}];
  fprintf('adding path %s\n',sgwpt_tmp_pathname);
  addpath(sgwpt_tmp_pathname);

end
clear sgwt_relpathlist sgwt_tmp_pathname
