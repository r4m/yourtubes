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

function varargout = poly(varargin)
%POLY Select polygonal region of interest.
%   Use POLY to select a polygonal region of interest within an
%   image. POLY returns a binary image that you can use as a mask for
%   masked filtering.
%
%   [x,y,BW,xi,yi] = POLY returns the XData and YData in x and y; the
%   mask image in BW; and the polygon coordinates in xi and yi.
%
%   [x,y,BW,xi,yi,roix,roiy] = POLY returns the XData and YData in x and y; the
%   mask image in BW; the polygon coordinates in xi and yi, the polygon
%   coordinates in xi and yi in pixels.
%
%   If POLY is called with no output arguments, the resulting image is
%   displayed in a new figure.
%
%   Class Support
%   -------------
%   The input image I can be uint8, uint16, int16, single or double.  The
%   output image BW is logical. All other inputs and outputs are double.

[xdata,ydata,num_rows,num_cols,xi,yi] = parse_inputs(varargin{:});

if length(xi)~=length(yi)
    eid = sprintf('Images:%s:xiyiMustBeSameLength',mfilename);
    error(eid,'%s','XI and YI must be the same length.');
end

% Make sure polygon is closed.
if (~isempty(xi))
    if ( xi(1) ~= xi(end) || yi(1) ~= yi(end) )
        xi = [xi;xi(1)];
        yi = [yi;yi(1)];
    end
end
% Transform xi,yi into pixel coordinates.
roix = axes2pix(num_cols, xdata, xi);
roiy = axes2pix(num_rows, ydata, yi);

d = poly2mask(roix, roiy, num_rows, num_cols);

switch nargout
    case 0
        figure
        if (~isequal(xdata, [1 size(d,2)]) || ~isequal(ydata, [1 size(d,1)]))
            imshow(d,'XData',xdata,'YData',ydata);  % makes tick labels visible
        else
            imshow(d)
        end

    case 1
        varargout{1} = d;
    case 2
        varargout{1} = d;
        varargout{2} = xi;

    case 3
        varargout{1} = d;
        varargout{2} = xi;
        varargout{3} = yi;

    case 4
        varargout{1} = xdata;
        varargout{2} = ydata;
        varargout{3} = d;
        varargout{4} = xi;

    case 5
        varargout{1} = xdata;
        varargout{2} = ydata;
        varargout{3} = d;
        varargout{4} = xi;
        varargout{5} = yi;
    case 6
        varargout{1} = xdata;
        varargout{2} = ydata;
        varargout{3} = d;
        varargout{4} = xi;
        varargout{5} = yi;
        varargout{6} = roix;
    case 7
        varargout{1} = xdata(2:end);
        varargout{2} = ydata(2:end);
        varargout{3} = d;
        varargout{4} = xi(2:end);
        varargout{5} = yi(2:end);
        varargout{6} = round(roix(2:end));
        varargout{7} = round(roiy(2:end));        
    otherwise
        eid = sprintf('Images:%s:tooManyOutputArgs',mfilename);
        error(eid,'%s','Too many output arguments');
end

%-------------------------------------------------------------------
function [x,y,nrows,ncols,xi,yi] = parse_inputs(varargin)

switch nargin
    case 0,
        % ROIPOLY
        %  Get information from the current figure
        [x,y,a,hasimage] = getimage;
        if ~hasimage,
            eid = sprintf('Images:%s:needImageInFigure',mfilename);
            error(eid,'%s',...
                'The current figure must contain an image to use ROIPOLY.');
        end
        [xi,yi] = getline(gcf,'closed'); % Get rect info from the user.
        nrows = size(a,1);
        ncols = size(a,2);
    case 1
        [x,y,a,hasimage] = getimage(varargin{1});
        if ~hasimage,
            eid = sprintf('Images:%s:needImageInFigure',mfilename);
            error(eid,'%s',...
                'The current figure must contain an image to use ROIPOLY.');
        end
        [xi,yi] = getline(varargin{1},'closed'); % Get rect info from the user.
        nrows = size(a,1);
        ncols = size(a,2);
    case 2
        a = varargin{1};
        nrows = size(a,1);
        ncols = size(a,2);
        x = [1 ncols];
        y = [1 nrows];
        %imshow(a);
        [xi,yi] = getline(varargin{2},'closed');
    otherwise,
        eid = sprintf('Images:%s:invalidInputArgs',mfilename);
        error(eid,'%s','Invalid input arguments.');
end

xi = cast_to_double(xi);
yi = cast_to_double(yi);
x = cast_to_double(x);
y = cast_to_double(y);
nrows= cast_to_double(nrows);
ncols = cast_to_double(ncols);

%-------------------------------------------------------------------
function a = cast_to_double(a)
  if ~isa(a,'double')
    a = double(a);
  end
  