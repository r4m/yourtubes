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

function varargout = table(varargin)
%IPTTABLE Display data in a two-column table.
%   HP = TABLE(HPARENT,DATA) displays the information in DATA
%   in a two-column table. TABLE creates the table in a uipanel
%   object and returns HP, a handle to this object. HPARENT specifies
%   the figure or uipanel object that will contain the table uipanel. 
%   
%   DATA can be a structure or a cell array.
%
%   If data is a structure, TABLE displays the structure fieldnames
%   in the first column of the table and the corresponding values in 
%   the second column of the table. The columns have the headings
%   "Fieldname" and "Value", respectively. 
%
%   Note: TABLE displays only the first level of the structure. 
%
%   If data is a cell array, the cell array must be N-by-2. TABLE 
%   displays the first element in each row of the cell array in the 
%   the first column of the table and the second element in each row
%   of the cell array in the second column of the table. The columns
%   have the headings "Attribute" and "Value", respectively.
%
%   Positioning
%   -----------
%   TABLE positions the uipanel object it creates in the lower-left
%   corner of HPARENT. TABLE sizes the uipanel to fit the amount of
%   information in the table, up to the maximum size of HPARENT. If 
%   the table doesn't fit, TABLE adds scroll bars to the table.
%
%   Depending on the amount of information to display, TABLE can 
%   appear to take over the entire figure. By default, HPANEL has
%   'Units' set to 'normalized' and 'Position' set
%   to [0 0 1 1]. If you want to see the other children of HPARENT, you
%   must manually set the 'Position' property of HPANEL.
%
%   [HP, HUIT] = TABLE(...) returns a handle to the uitable, HUIT,
%   contained in HP for testing purposes. HUIT is not accessible in the
%   graphics hierarchy.

[hparent,tableData,columnNames] = parseInputs(varargin{:});

%%
hp = uipanel('Parent', hparent,...
    'Units','pixels',...
    'BorderType','none',...
    'Visible','off');
%%
[huit,huic] = uitable('Data',tableData,...
    'ColumnNames',columnNames);

set(huic,'Parent',hp,...
    'Visible','off');
set(huit,'Editable', 0);  % must set separately b/c uitable not HG

% declare fudge factors so it has scope. these work b/c uitable font size
% cannot be changed.
fudge = 15;

setPositionOfPanel;
setPositionOfTable;

%set(hp,'ResizeFcn',@resizeTable);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setPositionOfPanel
        
        pos = getpixelposition(hparent);
        FigWidth = pos(3);
        FigHeight = pos(4);

        spaceWidth = 0.1*FigWidth;
        spaceHight = 0.05*FigHeight;
        
        %sizes = cellfun('prodofsize',columnNames);
        %maxWidth = sum(sizes);
        maxWidth = FigWidth - 2*spaceWidth;
        if isempty(maxWidth)
           maxWidth = 10;  %in event of empty tableData.
            % works b/c uitable's font size cannot be changed
        end
        
        NumRows = get(huit,'NumRows');
        if(NumRows < 10)
            maxHeight = NumRows  + 1; %include header row
        else
            maxHeight = 11; %include header row
        end
        
        if ispc
          fudgeFactor = 16;
        else
          fudgeFactor = 18;
        end
        
        hpPos = get(hp,'Position');
%         set(hp,'Position',[hpPos(1) ...
%             hpPos(2) ...
%             min(maxWidth*fudge,hpPos(3)) ...
%             min(maxHeight*fudgeFactor,hpPos(4))]);
%         end
        set(hp,'Position',[hpPos(1)+spaceWidth ...
            hpPos(2)+spaceHight ...
            maxWidth ...
            maxHeight*fudgeFactor+1]);
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setPositionOfTable

        %oldUnits = get(hp,'Units');
        %set(hp,'Units','pixels');
        hpPos = get(hp,'Position');
        set(huic,'Position',[1 1 hpPos(3) hpPos(4)]);

         columnOneWidth = floor((hpPos(3))/6) - 3;
         huit.setColumnWidth(0,columnOneWidth);
         columnTwoWidth = columnOneWidth;
         huit.setColumnWidth(1,columnTwoWidth);
         columnThreeWidth = columnOneWidth;
         huit.setColumnWidth(2,columnThreeWidth);
         columnFourWidth = columnOneWidth;
         huit.setColumnWidth(3,columnFourWidth);
         columnFiveWidth = columnOneWidth;
         huit.setColumnWidth(4,columnFiveWidth);
         columnSixWidth = columnOneWidth;
         huit.setColumnWidth(5,columnSixWidth);
        %set(hp,'Units',oldUnits);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function resizeTable(obj,evt)
%         setPositionOfTable;
%     end

set(huic,'Visible','on');
set(hp,'Visible','on');
varargout{1} = hp;
varargout{2} = huit;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hparent,tableData,colNames] = parseInputs(varargin)

iptchecknargin(3,3,nargin,mfilename);

hparent = varargin{1};
iptcheckhandle(hparent,{'figure','uipanel','uicontainer'},mfilename, ...
    'HPARENT',1);

tableDataProp = [];
tableDataReal = [];

if iscell(varargin{2}) && size(varargin{2},2) == 3     
    tableDataProp = createTableDataPropFromCellArray(varargin{2});    
% else
%     eid = sprintf('Images:%s:invalidInputArgument',mfilename);
%     msg = 'The second input argument must be a valid structure or ';
%     msg2 = 'N-by-3 cell array.';
%     error(eid,'%s%s',msg,msg2);
end

if iscell(varargin{3}) && size(varargin{3},2) == 3     
    tableDataReal = createTableDataRealFromCellArray(varargin{2});
%  else
%     eid = sprintf('Images:%s:invalidInputArgument',mfilename);
%     msg = 'The third input argument must be a valid structure or ';
%     msg2 = 'N-by-3 cell array.';
%     error(eid,'%s%s',msg,msg2);
end

tableData = [tableDataProp, tableDataReal];
colNames = {'Diameter [%]', 'x [%]', 'y [%]','Diameter [cm]', 'x [cm]', 'y [cm]'};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tableData = createTableDataPropFromCellArray(c)

diameterProp = c(:,1);
abscissaProp = c(:,2);
ordinateProp = c(:,3);
charArray1 = evalc('disp(diameterProp)');
dispOfValues1 = strread(charArray1,'%s','delimiter','\n');
charArray2 = evalc('disp(abscissaProp)');
dispOfValues2 = strread(charArray2,'%s','delimiter','\n');
charArray3 = evalc('disp(ordinateProp)');
dispOfValues3 = strread(charArray3,'%s','delimiter','\n');

numFields = length(diameterProp);
tableData = cell(numFields,6);

% First column of tableData contain diameter. Second column of tableData
% contains the string representation of abscissa. We use the values or
% dispOfValues depending on whether each element of values is a vector of
% characters.

for idx = 1: numFields
    val = diameterProp{idx};
    if ischar(val) && size(val,1) == 1
        tableData{idx,1} = val;
    else
        val = dispOfValues1{idx};
        spaces = isspace(val);  % Remove extra whitespace,e.g, [    8].
        val(spaces)= '';
        tableData{idx,1} = val;
    end
end
for idx = 1: numFields
    val = abscissaProp{idx};
    if ischar(val) && size(val,1) == 1
        tableData{idx,2} = val;
    else
        val = dispOfValues2{idx};
        spaces = isspace(val);  % Remove extra whitespace,e.g, [    8].
        val(spaces)= '';
        tableData{idx,2} = val;
    end
end
for idx = 1: numFields
    val = ordinateProp{idx};
    if ischar(val) && size(val,1) == 1
        tableData{idx,3} = val;
    else
        val = dispOfValues3{idx};
        spaces = isspace(val);  % Remove extra whitespace,e.g, [    8].
        val(spaces)= '';
        tableData{idx,3} = val;
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tableData = createTableDataRealFromCellArray(c)

diameterReal = c(:,4);
abscissaReal = c(:,5);
ordinateReal = c(:,6);
charArray4 = evalc('disp(diameterReal)');
dispOfValues4 = strread(charArray4,'%s','delimiter','\n');
charArray5 = evalc('disp(abscissaReal)');
dispOfValues5 = strread(charArray5,'%s','delimiter','\n');
charArray6 = evalc('disp(ordinateReal)');
dispOfValues6 = strread(charArray6,'%s','delimiter','\n');

numFields = length(diameterReal);
tableData = cell(numFields,6);

% First column of tableData contain diameter. Second column of tableData
% contains the string representation of abscissa. We use the values or
% dispOfValues depending on whether each element of values is a vector of
% characters.

for idx = 1: numFields
    val = diameterReal{idx};
    if ischar(val) && size(val,1) == 1
        tableData{idx,4} = val;
    else
        val = dispOfValues4{idx};
        spaces = isspace(val);  % Remove extra whitespace,e.g, [    8].
        val(spaces)= '';
        tableData{idx,4} = val;
    end
end
for idx = 1: numFields
    val = abscissaReal{idx};
    if ischar(val) && size(val,1) == 1
        tableData{idx,5} = val;
    else
        val = dispOfValues5{idx};
        spaces = isspace(val);  % Remove extra whitespace,e.g, [    8].
        val(spaces)= '';
        tableData{idx,5} = val;
    end
end
for idx = 1: numFields
    val = ordinateReal{idx};
    if ischar(val) && size(val,1) == 1
        tableData{idx,6} = val;
    else
        val = dispOfValues6{idx};
        spaces = isspace(val);  % Remove extra whitespace,e.g, [    8].
        val(spaces)= '';
        tableData{idx,6} = val;
    end
end
end
