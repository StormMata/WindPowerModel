clearvars -except I AeroDyn A Output Hannah HannahSS SSs SSsn speed uRel vRel Hphi aoa Hchord Htwist Hpower HW sweep
clc

addpath(['/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/' ...
         'Courses/Research/WindPowerModel/5 MW Turbine']);

%% Inputs

load('Turbine.mat')                                                         % Load turbine parameters

load('Airfoils.mat')                                                        % Load airfoil parameters

% load('WindProfiles.mat')                                                    % Load atmospheric conditions

%% Interpolate Blade Parameters

R         = 1.5:1:T.R;                                                      % Evenly spaced radial positions    [m]

Chord     = interp1(AeroDyn.r,AeroDyn.Chord,R,'linear','extrap');           % Interpolated chord length         [m]

Twist     = interp1(AeroDyn.r,AeroDyn.Twist,R,'linear','extrap');           % Interpolated twist angle          [deg]

AeroIndex = interp1(AeroDyn.r,AeroDyn.AeroIndex,R,'nearest','extrap');      % Interpolated airfoil index        [-]

%% Mesh Rotor Area

[ycoords,zcoords] = meshgrid(-T.R:T.R,(T.Hub+T.R):-1:(T.Hub-T.R));          % Generate grid of points in cartesian coordinates

[azimuth,R]       = cart2pol(ycoords,zcoords-T.Hub);                        % Convert to polar coordinates

azimuth           = fliplr(mod(azimuth - pi/2,2*pi));                       % Fix azimuth angle 0 position

%% Generate Freestream Velocity Field

a = 1/3;

U      = interp1(I.Heights,I.SpeedShear,((T.Hub+T.R):-1:(T.Hub-T.R))','linear') .* ones(size(azimuth));

LocYaw = deg2rad(interp1(I.Heights,I.DirShear,((T.Hub+T.R):-1:(T.Hub-T.R))','linear') .* ones(size(azimuth)));

U_axi  = (U*(1-a) .* cos(LocYaw .* sin(azimuth)) .* cos(LocYaw .* cos(azimuth)));

U_tan  = (T.Omega .* R - U*(1-a) .* cos(LocYaw * sin(azimuth)) .* sin(LocYaw .* cos(azimuth)));

W      = U_axi.^2 + U_tan.^2;

Phi    = atan2(U_axi, U_tan);

Chord  = interp1(AeroDyn.r,AeroDyn.Chord,R,'linear','extrap');

Twist     = interp1(AeroDyn.r,AeroDyn.Twist,R,'linear','extrap');

AeroIndex = interp1(AeroDyn.r,AeroDyn.AeroIndex,R,'nearest','extrap');      % Interpolated airfoil index        [-]

alpha  = rad2deg(Phi) - Twist;



% plot(AeroDyn.r,AeroDyn.Chord,'b')
% hold on
% yline(0)
% axis equal
% plot(0:63,Cint,'r')

