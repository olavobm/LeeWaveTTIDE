function jkmodel = JK_xy2lonlat(jkmodel)
% jkmodel = JK_XY2LONLAT(jkmodel)
%
%
%
% Uses jkmodel.x and jkmodel.y to compute
% longitude and latitude.
%
%
%

[X, Y] = meshgrid(jkmodel.x, jkmodel.y);
XX = X+sqrt(-1)*Y;
XX = XX.*exp(-sqrt(-1)*12*pi/180);

kmpernm = 1.8532;
dtheta = 12;
cenlon = 148.;
cenlat = -44.;

%
jkmodel.lon = real(XX)/kmpernm/60./cos(cenlat*pi/180.)+cenlon;
jkmodel.lat = imag(XX)/kmpernm/60.+cenlat;
