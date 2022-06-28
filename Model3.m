clearvars -except Hannah HannahSS speed uRel vRel Hphi aoa Hchord Htwist Hpower HW sweep I
% clc

addpath(['/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/' ...
         'Courses/Research/WindPowerModel/5 MW Turbine']);

%% Inputs

load('Turbine.mat')                                                         % Load turbine parameters

load('Airfoils.mat')                                                        % Load airfoil parameters

load('WindProfiles.mat')                                                    % Load atmospheric conditions

%% Mesh Rotor Area

% spacing = 1/0.5;
% N = 128;
% 
% [ycoords,zcoords] = meshgrid(linspace(-1,1,N*2+1),linspace(1,-1,N*2+1));          % Generate grid of points in cartesian coordinates
% 
% [azimuth,R]       = cart2pol(ycoords,zcoords);                        % Convert to polar coordinates
% 
% azimuth           = fliplr(mod(azimuth - pi/2,2*pi));                       % Fix azimuth angle 0 position

%% Generate Freestream Velocity Field

a         = 1/3;
rho       = 1.225;

for i = 1:size(I.SpeedShear,2)
    i

U_Hub     = interp1(I.Heights,I.SpeedShear(:,i),T.Hub,'linear');
U         = interp1(I.Heights,I.SpeedShear(:,i),((T.Hub+T.R):-spacing:(T.Hub-T.R))','linear') .* ones(size(azimuth));

LocYaw    = deg2rad(interp1(I.Heights,I.DirShear(:,i),((T.Hub+T.R):-spacing:(T.Hub-T.R))','linear') .* ones(size(azimuth)));

U_axi     = (U*(1-a) .* cos(LocYaw .* sin(azimuth)) .* cos(LocYaw .* cos(azimuth)));

U_tan     = (T.Omega .* R - U*(1-a) .* cos(LocYaw * sin(azimuth)) .* sin(LocYaw .* cos(azimuth)));

U_rel     = U_axi.^2 + U_tan.^2;

Phi       = atan2(U_axi, U_tan);

Chord     = interp1(AeroDyn.r,AeroDyn.Chord,R,'linear','extrap');

Twist     = interp1(AeroDyn.r,AeroDyn.Twist,R,'linear','extrap');

AeroIndex = interp1(AeroDyn.r,AeroDyn.AeroIndex,R,'nearest','extrap');      % Interpolated airfoil index        [-]

alpha     = rad2deg(Phi) - Twist;

for m = 1:size(AeroIndex,1)
    for n = 1:size(AeroIndex,2)
        CL(m,n) = interp1(A.(sprintf('A%i',AeroIndex(m,n))).Alpha,A.(sprintf('A%i',AeroIndex(m,n))).CL,alpha(m,n),'linear');
        CD(m,n) = interp1(A.(sprintf('A%i',AeroIndex(m,n))).Alpha,A.(sprintf('A%i',AeroIndex(m,n))).CD,alpha(m,n),'linear');
    end
end

Q = T.B .* 0.5 .* rho .* Chord .* U_rel .* R .* (CL .* sin(Phi) - CD .* cos(Phi));

P = T.Omega .* Q;

OutterFilter = R >= T.R_hub;   % Exclude hub area

InnerFilter  = R <= T.R;   % Only capture values within the outer rotor radius

P = P .* OutterFilter .* InnerFilter;

area = sum(OutterFilter .* InnerFilter, 'all') * spacing^2

ri     = linspace(1.5,T.R,length(R));
thetai = linspace(0,2*pi,length(R));

Power(i) = trapz(thetai,trapz(ri,P,2))/(2*pi);
Power(i) = trapz(1:length(R),trapz(1:length(R),P,2))/(pi*T.R^2-pi*T.R_hub^2);
Power(i) = sum(P,'all')/(2*pi);
trapz(1:length(R),trapz(1:length(R),P,2))/(2*pi)
trapz(trapz(P))/area
% Power(i) = trapz(thetai,trapz(ri,P,2))/(pi*T.R^2-pi*T.R_hub^2);
% Power(i) = sum(P,'all')/(pi*T.R^2-pi*T.R_hub^2);

end

% P(P==0) = NaN;
% 
% % trapz(linspace(-T.R,T.R,length(R)),trapz(linspace(-T.R,T.R,length(R)),P))/(pi*T.R^2-pi*T.R_hub^2)
% 
% trapz(linspace(1,T.R,length(R)),trapz(linspace(T.Hub-T.R,T.Hub+T.R,length(R)),P,1),2)/(pi*T.R^2-pi*T.R_hub^2)
% 
% trapz(linspace(-T.R,T.R,length(R)),trapz(linspace(-T.R,T.R,length(R)),P,1),2)/(pi*T.R^2-pi*T.R_hub^2)
% 
% trapz(linspace(T.Hub+T.R,T.Hub-T.R,length(R))',P,1)/(pi*T.R^2-pi*T.R_hub^2)
% 
% trapz(trapz(P.*azimuth,2).*R,1)/(pi*T.R^2-pi*T.R_hub^2)

% ri     = linspace(1.5,T.R,length(R));
% thetai = linspace(0,2*pi,length(R));
% 
% trapz(thetai,trapz(ri,P,2))/(2*pi)

% 2*pi*trapz(ri,ri.*trapz(thetai,P,1))/(pi*T.R^2-pi*T.R_hub^2)

% P = sum(P,'all')*2*pi/(pi*T.R^2-pi*T.R_hub^2)

% plot(AeroDyn.r,AeroDyn.Chord,'b')
% hold on
% yline(0)
% axis equal
% plot(0:63,Cint,'r')

figure;
surf(ycoords(1,:),zcoords(:,1),P,'EdgeColor','none')
xlim([min(min(ycoords)) max(max(ycoords))])
ylim([min(min(zcoords)) max(max(zcoords))])
xlabel('y (m)')
ylabel('z (m)')