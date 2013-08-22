% Copyright (c) 2008, Department of Information Engineering, University of Padova.
% All rights reserved.
% 
% This file is part of YourTubes.
% 
% YourTubes is free software: you can redistribute it and/or modify it under the terms
% of the GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% YourTubes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
% PURPOSE.  See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with YourTubes.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
%  @date 28/04/2008
%  @author Filippo Zanella <filippo.zanella@dei.unipd.it>
%
% ===================================================================================

function about(varargin)
%ABOUT About the YourTubes Toolbox.
%   ABOUT displays the version number of the YourTubes Toolbox
%  and the copyright notice in a modal dialog box.
 
tlbxName = 'YourTubes';
tlbxVersion = '08.04.29 beta';
str = sprintf('%s %s\nCopyright 2008 DEI, Unipd.', ...
              tlbxName, tlbxVersion);
msgbox(str,tlbxName,'help','modal');

