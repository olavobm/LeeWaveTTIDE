function JKsection = JK_slicebotdepth(JKsection)
% JKsection = JK_SLICEBOTDEPTH(JKsection)
%
%   inputs
%       - JKsection: structure with the variables saved
%                    in the Slice****km.mat files.
%
%   outputs
%       - JKsection: same structure, with "bottomdepth"
%                    as an additional field.
%
%
%
%
%
% Olavo Badaro Marques, 30/Jul/2018.


%%

%
JKsection.bottomdepth = NaN(1, length(JKsection.x));

%
bla = squeeze(JKsection.U(1, :, :));

%
for i = 1:size(bla, 2)
    
    %
	indbottom_aux = find(bla(:, i)==0, 1, 'first');
   
    %
    JKsection.bottomdepth(i) = JKsection.z(indbottom_aux);
end
