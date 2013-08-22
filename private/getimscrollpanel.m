function hpanel = getimscrollpanel(himage)
%GETIMSCROLLPANEL Get image scrollpanel.
%   HPANEL = GETIMSCROLLPANEL(HIMAGE) returns the imscrollpanel associated
%   with HIMAGE. If no imscrollpanel is found, GETIMSCROLLPANEL returns an
%   empty matrix.
  
%   Copyright 2004-2005 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2005/03/31 16:33:21 $

iptcheckhandle(himage,{'image'},mfilename,'HIMAGE',1);
hScrollable = ancestor(himage,'uipanel');

% Not using ancestor here because of a possible HG bug.
% hpanel = ancestor(hScrollable,'uipanel') is returning hScrollable
hpanel = get(hScrollable,'Parent');  

% validate hpanel is a scrollpanel by checking hierarchy
if ~isempty(hpanel)

  kids = get(hpanel,'children');

  if numel(kids)~=4 
    hpanel = [];
    return
  end

  types = get(kids,'Type');
  if ~isequal(types,{'uipanel','uicontrol','uicontrol','uicontrol'}')
    hpanel = [];
    return
  end
  
  styles = get(kids(2:4),'Style');
  if ~isequal(styles,{'frame','slider','slider'}')
    hpanel = [];
    return
  end
  
  grandkid = get(kids(1),'children');
  if numel(grandkid)~=1 || ~strcmp('axes',get(grandkid,'Type'))
    hpanel = [];
    return
  end
  
  greatgrandkids = get(grandkid,'children');
  greatgrandkid_image = findobj(greatgrandkids,'Type','image');
  if numel(greatgrandkid_image)~=1 || greatgrandkid_image~=himage
    hpanel = [];
    return
  end

end
