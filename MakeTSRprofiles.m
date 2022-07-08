%% TSR Sensitivity

sp = ('/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Research/WindPowerModel/5 MW Turbine/TSR Sensitivity');

addpath(sp);

load('Turbine.mat')

% U = 5 m/s

    U = 5;
    
    TSR = 10;

    [speed, direc] = meshgrid(-0.1:0.1:0.7,0.6:-0.1:-0.1);
    index          = 1;

    I.Heights = (200:-1:1)';
    
    for i = 1:size(direc,1)
        for j = 1:size(speed,2)
            I.SpeedShear(:,index) = U * (I.Heights./T.Hub).^speed(i,j);
            I.DirShear(:,index)   = direc(i,j)*I.Heights - direc(i,j)*T.Hub;
            I.Omega(index)        = U*TSR/T.R;
            index = index + 1;
        end
    end

    save(sprintf('%s/TSRsweep_U%1.0i.mat',sp,U),'I')

% U = 6 m/s

    U = 6;
    
    TSR = 10 - 5 * (0.3544/1.3922);

    [speed, direc] = meshgrid(-0.1:0.1:0.7,0.6:-0.1:-0.1);
    index          = 1;

    I.Heights = (200:-1:1)';
    
    for i = 1:size(direc,1)
        for j = 1:size(speed,2)
            I.SpeedShear(:,index) = U * (I.Heights./T.Hub).^speed(i,j);
            I.DirShear(:,index)   = direc(i,j)*I.Heights - direc(i,j)*T.Hub;
            I.Omega(index)        = U*TSR/T.R;
            index = index + 1;
        end
    end

    save(sprintf('%s/TSRsweep_U%1.0i.mat',sp,U),'I')

% U = 7 m/s

    U = 7;
    
    TSR = 10 - 5 * (0.5667/1.3922);

    [speed, direc] = meshgrid(-0.1:0.1:0.7,0.6:-0.1:-0.1);
    index          = 1;

    I.Heights = (200:-1:1)';
    
    for i = 1:size(direc,1)
        for j = 1:size(speed,2)
            I.SpeedShear(:,index) = U * (I.Heights./T.Hub).^speed(i,j);
            I.DirShear(:,index)   = direc(i,j)*I.Heights - direc(i,j)*T.Hub;
            I.Omega(index)        = U*TSR/T.R;
            index = index + 1;
        end
    end

    save(sprintf('%s/TSRsweep_U%1.0i.mat',sp,U),'I')

% U = 8 - 10 m/s

    U = 9;
    
    TSR = 10 - 5 * (0.6847/1.3922);

    [speed, direc] = meshgrid(-0.1:0.1:0.7,0.6:-0.1:-0.1);
    index          = 1;

    I.Heights = (200:-1:1)';
    
    for i = 1:size(direc,1)
        for j = 1:size(speed,2)
            I.SpeedShear(:,index) = U * (I.Heights./T.Hub).^speed(i,j);
            I.DirShear(:,index)   = direc(i,j)*I.Heights - direc(i,j)*T.Hub;
            I.Omega(index)        = U*TSR/T.R;
            index = index + 1;
        end
    end

    save(sprintf('%s/TSRsweep_U%1.0i.mat',sp,U),'I')

% U = 11 m/s

    U = 11;
    
    TSR = 10 - 5 * (0.8027/1.3922);

    [speed, direc] = meshgrid(-0.1:0.1:0.7,0.6:-0.1:-0.1);
    index          = 1;

    I.Heights = (200:-1:1)';
    
    for i = 1:size(direc,1)
        for j = 1:size(speed,2)
            I.SpeedShear(:,index) = U * (I.Heights./T.Hub).^speed(i,j);
            I.DirShear(:,index)   = direc(i,j)*I.Heights - direc(i,j)*T.Hub;
            I.Omega(index)        = U*TSR/T.R;
            index = index + 1;
        end
    end

    save(sprintf('%s/TSRsweep_U%1.0i.mat',sp,U),'I')

% U = 11.4 m/s

    U = 11.4;
    
    TSR = 10 - 5 * (0.8499/1.3922);

    [speed, direc] = meshgrid(-0.1:0.1:0.7,0.6:-0.1:-0.1);
    index          = 1;

    I.Heights = (200:-1:1)';
    
    for i = 1:size(direc,1)
        for j = 1:size(speed,2)
            I.SpeedShear(:,index) = U * (I.Heights./T.Hub).^speed(i,j);
            I.DirShear(:,index)   = direc(i,j)*I.Heights - direc(i,j)*T.Hub;
            I.Omega(index)        = U*TSR/T.R;
            index = index + 1;
        end
    end

    save(sprintf('%s/TSRsweep_U114.mat',sp),'I')

% U = 12 m/s

    U = 12;
    
    TSR = 10 - 5 * (0.9443/1.3922);

    [speed, direc] = meshgrid(-0.1:0.1:0.7,0.6:-0.1:-0.1);
    index          = 1;

    I.Heights = (200:-1:1)';
    
    for i = 1:size(direc,1)
        for j = 1:size(speed,2)
            I.SpeedShear(:,index) = U * (I.Heights./T.Hub).^speed(i,j);
            I.DirShear(:,index)   = direc(i,j)*I.Heights - direc(i,j)*T.Hub;
            I.Omega(index)        = U*TSR/T.R;
            index = index + 1;
        end
    end

    save(sprintf('%s/TSRsweep_U%1.0i.mat',sp,U),'I')

% U = 13 m/s

    U = 13;
    
    TSR = 10 - 5 * (1.0859/1.3922);

    [speed, direc] = meshgrid(-0.1:0.1:0.7,0.6:-0.1:-0.1);
    index          = 1;

    I.Heights = (200:-1:1)';
    
    for i = 1:size(direc,1)
        for j = 1:size(speed,2)
            I.SpeedShear(:,index) = U * (I.Heights./T.Hub).^speed(i,j);
            I.DirShear(:,index)   = direc(i,j)*I.Heights - direc(i,j)*T.Hub;
            I.Omega(index)        = U*TSR/T.R;
            index = index + 1;
        end
    end

    save(sprintf('%s/TSRsweep_U%1.0i.mat',sp,U),'I')

% U = 14 m/s

    U = 14;
    
    TSR = 10 - 5 * (1.2039/1.3922);

    [speed, direc] = meshgrid(-0.1:0.1:0.7,0.6:-0.1:-0.1);
    index          = 1;

    I.Heights = (200:-1:1)';
    
    for i = 1:size(direc,1)
        for j = 1:size(speed,2)
            I.SpeedShear(:,index) = U * (I.Heights./T.Hub).^speed(i,j);
            I.DirShear(:,index)   = direc(i,j)*I.Heights - direc(i,j)*T.Hub;
            I.Omega(index)        = U*TSR/T.R;
            index = index + 1;
        end
    end

    save(sprintf('%s/TSRsweep_U%1.0i.mat',sp,U),'I')