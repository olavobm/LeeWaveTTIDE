function hcb = addColorbar(haxs, dcb, wcb)
% hcb = ADDCOLORBAR(haxs, dcb, wcb)
%
%   inputs
%       - haxs:
%       - dcb
%       - wcb
%
%   outputs
%       - hcb:
%
%
%
% Olavo Badaro Marque, 25/Sep/2020.


%%

if length(haxs)==1
    
elseif length(haxs)==2
    
else
    
    error('Unaccepted number of axes handles.')
    
end

%% Check arrangement of axes, if multiple are given.
%
% (for now, assume they are vertically stacked)


%%

if length(haxs)==1
	%
    indaxscb = 1;
    indaxscb_other = 1;
    
elseif length(haxs)==2
    
    %
    indaxscb = 2;
    indaxscb_other = 1;
    
end





%% assume they are vertically stacked

%
hcb = colorbar(haxs(indaxscb));

%
hcb.Position = [sum(haxs(indaxscb).Position([1, 3])) + dcb, ...
                haxs(indaxscb).Position(2), wcb, ...
                sum(haxs(indaxscb_other).Position([2, 4])) - haxs(indaxscb).Position(2)];





