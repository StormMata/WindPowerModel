clearvars -except
clc

addpath(['/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/' ...
         'Courses/Research/WindPowerModel/5 MW Turbine']);

%% Inputs

load('Turbine.mat')                                                         % Load turbine parameters

load('Airfoils.mat')                                                        % Load airfoil parameters

% load('Uniform_Rated.mat')                                                   % Load atmospheric conditions
% load('R2_Sweep.mat')                                                        % Load atmospheric conditions

% sp = ('/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Research/WindPowerModel/5 MW Turbine/DShear Sensitivity');
sp = ('/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Research/WindPowerModel/5 MW Turbine/TSR Sensitivity');
addpath(sp)
fname = ('TSRsweep_U9.mat');
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
    
            U_axi(j,i)  = (U(j,i) * (1 - a(i)) * cos(betar * sin(azimuth(j))) * cos(betar * cos(azimuth(j))));
    
            U_tan(j,i)  = (I.Omega(index) * r(i) - U(j,i) * (1 - a(i)) * cos(betar * sin(azimuth(j))) * sin(betar * cos(azimuth(j))));
        
            W(j,i)      = U_axi(j,i)^2 + U_tan(j,i)^2;
        
            Phi(j,i)    = atan2(U_axi(j,i), U_tan(j,i));
        
            Chord(j,i)  = interp1(AeroDyn.r,AeroDyn.Chord,r(i),'linear','extrap');
    
            Twist(j,i)  = interp1(AeroDyn.r,AeroDyn.Twist,r(i),'linear','extrap');
        
            Alpha(j,i)  = rad2deg(Phi(j,i)) - Twist(j,i);
    
            AeroID(j,i) = interp1(AeroDyn.r,AeroDyn.AeroIndex,r(i),'nearest','extrap');      % Interpolated airfoil index        [-]
    
            CL(j,i)     = interp1(A.(sprintf('A%i',AeroID(j,i))).Alpha,A.(sprintf('A%i',AeroID(j,i))).CL,Alpha(j,i),'linear');
            CD(j,i)     = interp1(A.(sprintf('A%i',AeroID(j,i))).Alpha,A.(sprintf('A%i',AeroID(j,i))).CD,Alpha(j,i),'linear');
    
            dQ(j,i)     = T.B * 0.5 * T.rho * Chord(j,i) * W(j,i) * r(i) * (CL(j,i) * sin(Phi(j,i)) - CD(j,i) * cos(Phi(j,i)));
    
            dP(j,i)     = dQ(j,i) * I.Omega(index);

        end

    end
    
    P(index) = trapz(r,trapz(azimuth,dP,1)) / (1e6 * 2 * pi);

end

P = reshape(P,[size(speed,2),size(speed,1)])';

P_LH = P;

P_uni = P(7,2);

low  = min(P,[],'all');                                       % Find minimum
high = max(P,[],'all');                                       % Find maximum

P_LH(P < P_uni & ~isnan(P)) = low;       % Convert all values >1 to maximum
P_LH(P >= P_uni)            = high;                               % Convert all values <1 to minimum

% if abs(high-1) > abs(1-low)                                             % Calibrate colormap values
%     P_LH(end,end) = 1-abs(high-1);
% else
%     P_LH(end,end) = 1+abs(low-1);
% end

%% Make Alpha-Beta Plot

AB = figure;
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

AB_LH = figure;
    imagesc(speed(1,:),direc(:,1),P_LH)
    axis xy
    axis equal
    xlim([min(speed(1,:))-0.05 max(speed(1,:))+0.05])
    colormap([0.66 0.66 0.66; .33 .33 .33])
    cbar = colorbar;
    xlabel('Shear Exponent')
    ylabel('Direction Shear (deg/m)')
    cbar.Label.String = 'Aerodynamic Power (MW)';
    title(fname, 'Interpreter', 'none');
    savefig(AB_LH,sprintf('%s/%s_%s',sp,erase(fname,'.mat'),'AB_LH'))

%% Make Surface Plot

AB_surf = figure;
    % surf(speed,direc,P/P(7,2),'LineStyle','none')
    surf(speed,direc,P,'LineStyle','none')
    colormap('parula')
    xlabel('Shear Exponent')
    ylabel('Direction Shear (deg/m)')
    zlabel('Aerodynamic Power (MW)')
    title(fname, 'Interpreter', 'none');
    savefig(AB_surf,sprintf('%s/%s_%s',sp,erase(fname,'.mat'),'AB_surf'))

% [Ra, Th] = meshgrid(r,azimuth);
% 
% X = Ra.*sin(Th);
% Y = Ra.*cos(Th);
% 
% figure;
% surf(X,Y+T.Hub,P,'LineStyle','none')
% xlabel('y (m)')
% ylabel('z (m)')
% axis equal
% view(gca,0,90)

% %% Make Direction Shear Plot
% 
% Dplot = figure;
%     plot(I.DirShear(:,[1 10 19 27 37 41 51 61 71]),I.Heights)
%     hold on
%     yline(T.Hub+T.R)
%     yline(T.Hub-T.R)
%     title(fname, 'Interpreter', 'none');
%     savefig(Dplot,sprintf('%s/%s_%s',sp,erase(fname,'.mat'),'DShear'))
%     legend('-0.1 deg/m','','','','','','','0.6 deg/m','Location','northwest')
%     xlabel('Wind Direction (deg)')
%     ylabel('z (m)')
% 
% Dplot = figure;
%     plot(I.SpeedShear(:,1:10:71),I.Heights)
%     hold on
%     yline(T.Hub+T.R)
%     yline(T.Hub-T.R)
%     title(fname, 'Interpreter', 'none');
%     savefig(Dplot,sprintf('%s/%s_%s',sp,erase(fname,'.mat'),'DShear'))
%     legend('-0.1 deg/m','','','','','','','0.6 deg/m','Location','northwest')
%     xlabel('Wind Direction (deg)')
%     ylabel('z (m)')
