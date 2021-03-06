clearvars -except
clc

addpath(['/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/' ...
         'Courses/Research/WindPowerModel/5 MW Turbine']);

%% Inputs

load('Turbine.mat')                                                         % Load turbine parameters

load('Airfoils.mat')                                                        % Load airfoil parameters

% load('Uniform_Rated.mat')                                                   % Load atmospheric conditions
% load('R2_Sweep.mat')                                                        % Load atmospheric conditions

sp = ('/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Research/WindPowerModel/5 MW Turbine/DShear Sensitivity');
% sp = ('/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Research/WindPowerModel/5 MW Turbine/TSR Sensitivity');
addpath(sp)
fname = ('AlphaBetaSweep_Linear.mat');
load(fname)

%% Calculate

yaw     = 0.0;
azimuth = linspace(0,2 * pi,50);
r       = linspace(T.R_hub,T.R,100);
a       = 1/3 * ones(1, length(r));

[U, U_axi, U_tan, W, Phi, Chord, Twist, Alpha, AeroID, CL, CD, dQ, dP] = deal(zeros(length(azimuth),length(r)));
P = zeros(1, size(I.SpeedShear,2));
[speed, direc] = meshgrid(-0.1:0.1:0.7,0.6:-0.1:-0.1);

for index = 1:size(I.SpeedShear,2)
    index
    for j = 1:length(azimuth)
        for i = 1:length(r)
    
            zloc        = T.Hub + r(i) * cos(azimuth(j));
    
            betar       = deg2rad(interp1(I.Heights,I.DirShear(:,index),zloc,'linear','extrap'));
        
            U(j,i)      = interp1(I.Heights,I.SpeedShear(:,index),zloc,'linear');
    
            W(j,i)      = (U(j,i) * (1 - a(i)) * cos(betar * sin(azimuth(j))) * cos(betar * cos(azimuth(j))));

            vel(j,i)    = (U(j,i)*cos(betar));
    
%             U_tan(j,i)  = (I.Omega(index) * r(i) - U(j,i) * (1 - a(i)) * cos(betar * sin(azimuth(j))) * sin(betar * cos(azimuth(j))));
        
%             W(j,i)      = sqrt(U_axi(j,i)^2 + U_tan(j,i)^2)^3;
        
%             Phi(j,i)    = atan2(U_axi(j,i), U_tan(j,i));
%         
%             Chord(j,i)  = interp1(AeroDyn.r,AeroDyn.Chord,r(i),'linear','extrap');
%     
%             Twist(j,i)  = interp1(AeroDyn.r,AeroDyn.Twist,r(i),'linear','extrap');
%         
%             Alpha(j,i)  = rad2deg(Phi(j,i)) - Twist(j,i);
%     
%             AeroID(j,i) = interp1(AeroDyn.r,AeroDyn.AeroIndex,r(i),'nearest','extrap');      % Interpolated airfoil index        [-]
%     
%             CL(j,i)     = interp1(A.(sprintf('A%i',AeroID(j,i))).Alpha,A.(sprintf('A%i',AeroID(j,i))).CL,Alpha(j,i),'linear');
%             CD(j,i)     = interp1(A.(sprintf('A%i',AeroID(j,i))).Alpha,A.(sprintf('A%i',AeroID(j,i))).CD,Alpha(j,i),'linear');
%     
%             dQ(j,i)     = T.B * 0.5 * T.rho * Chord(j,i) * W(j,i) * r(i) * (CL(j,i) * sin(Phi(j,i)) - CD(j,i) * cos(Phi(j,i)));
%     
%             dP(j,i)     = dQ(j,i) * I.Omega(index);

        end

    end
    
    P(index) = ((trapz(r,trapz(azimuth,vel,1))/(pi*T.R^2 - pi*T.R_hub^2)))^3;

%     P(index) = (trapz(r,trapz(azimuth,W,1))^3/(pi*T.R^2 - pi*T.R_hub^2));

end

P = reshape(P,[size(speed,2),size(speed,1)])';

%% Make Alpha-Beta Plot

% AB = figure;
figure;
    imagesc(speed(1,:),direc(:,1),P)
    axis xy
    axis equal
    xlim([min(speed(1,:))-0.05 max(speed(1,:))+0.05])
    colormap('parula')
    cbar = colorbar;
    xlabel('Shear Exponent')
    ylabel('Direction Shear (deg/m)')
    cbar.Label.String = 'Aerodynamic Power (MW)';
    title(fname, 'Interpreter', 'none');
    savefig(AB,sprintf('%s/%s_%s',sp,erase(fname,'.mat'),'AB'))

[Ra, Th] = meshgrid(r,azimuth);

X = Ra.*sin(Th);
Y = Ra.*cos(Th);

figure;
surf(X,Y+T.Hub,W,'LineStyle','none')
xlabel('y (m)')
ylabel('z (m)')
axis equal
view(gca,0,90)

