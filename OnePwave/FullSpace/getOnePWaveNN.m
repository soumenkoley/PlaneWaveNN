clear; close all;
fAll = 2:0.2:7.0;% frequency vector in Hz
nR = 300; % number of points along radial direction
R = 8000; % maximum radius of integration in meters
nTheta = 200; % numbver of points along theta from [0,pi]
nPhi = 500; % number of points in phi from [0, 2pi]
rCav = 20; % radius of cavity in meters
zCav = 0; % location of cavity in z, (x,y) assumed to be zero
% note that since it is full space, the cavity is at the center
% so zCav = 0

% wave direction in degrees
thetaW = 90; % this can be different but will lead the same
% because of spherical spymmetry
phiW = 0; % again can be anything, does not metter, unit in degrees
vP = 4000; % P-wave velocity in m/s
rho = 2800; % same as the lower limit paper in kg/m^3
G = 6.67430*10^-11;
r1 = rCav; r2 = R; % these values of r1 and r2 will be used for the
% theoretical formula

bounds = [0,rCav,R];
% note that 10+290 = 300 = nR, please respect this
pointsPerRegion = [10, 290];  % finer sampling near cavity

%pointsPerRegion = [30,30,240];

rVec = [];
wR = [];

for i = 1:2
    [r_sub, w_sub] = lgwt(pointsPerRegion(i), bounds(i), bounds(i+1));
    rVec = [rVec; r_sub];
    wR = [wR; w_sub];
end

%[rVec, wR] = lgwt(nR, 0, R);

% -- Gauss–Legendre in [0, pi]
[theta, wTheta] = lgwt(nTheta, 0, pi);

% -- Trapezoidal in [0, 2pi)
phi = linspace(0, 2*pi, nPhi+1); phi(end) = [];
wPhi = 2*pi / nPhi;

% -- Build flat theta–phi grid
[TH, PH] = meshgrid(theta, phi);      % size Nphi × Ntheta
thetaFlat = TH(:);                   % (Nphi*Ntheta) × 1
phiFlat = PH(:);

[thetaPhi, weights] = createSurfPoints(max(rVec), nTheta, nPhi);
sThetaCPhi = sin(thetaPhi(:,1)).*cos(thetaPhi(:,2));
sThetaSPhi = sin(thetaPhi(:,1)).*sin(thetaPhi(:,2));
cTheta = cos(thetaPhi(:,1));
dsUnitVec = [sThetaCPhi,sThetaSPhi,cTheta];

% -- Precompute angular components
sinTheta = sin(thetaFlat);
cosTheta = cos(thetaFlat);
cosPhi   = cos(phiFlat);
sinPhi   = sin(phiFlat);

% -- Repeat Gauss–Legendre weights without interpolation
wtTheta = repmat(wTheta(:)', nPhi, 1);
wtTheta = wtTheta(:);  % Flattened to match thetaFlat

% -- Angular weights (excluding r², added later)
angularWeights = sinTheta .* wtTheta * wPhi;  % size (Nphi*Ntheta) × 1

xyzCav = [];
IVolTotAll = [];
for fNo = 1:1:length(fAll)
    IVolTot = [0,0,0];
    f = fAll(fNo);
    for i = 1:nR
        r = rVec(i);
        wr = wR(i);
        
        % Cartesian coordinates of current shell
        x = r * sinTheta .* cosPhi;
        y = r * sinTheta .* sinPhi;
        z = r * cosTheta;
        
        % Mask out cavity
        insideCavity = (x.^2 + y.^2 + (z - zCav).^2) < rCav^2;
%         xyzCav = [xyzCav;[x(insideCavity,1),y(insideCavity,1),z(insideCavity,1)]];
        % Exclude points inside cavity
        x(insideCavity) = [];
        y(insideCavity) = [];
        z(insideCavity) = [];
        w = angularWeights(~insideCavity);
        
        % Evaluate function
        IVol = getVolNN([x, y, z],zCav,vP,f,thetaW,phiW);  % returns vector same size as x
        
        % Accumulate contribution from this r shell
        IVolTot(1,1) = IVolTot(1,1) + sum(IVol(:,1) .* r^2 .* w) * wr;
        IVolTot(1,2) = IVolTot(1,2) + sum(IVol(:,2) .* r^2 .* w) * wr;
        IVolTot(1,3) = IVolTot(1,3) + sum(IVol(:,3) .* r^2 .* w) * wr;
        
        if(i==nR)
            % evaluate the surface integral
            [ampOut] = getPWaveAmp(x,y,z,vP,f,thetaW,phiW);
            [ISurf] = getSurfNN(ampOut,x,y,z,zCav,dsUnitVec);
    
            % Perform the integral
            ISurfTot(1,1) = sum(ISurf(:,1).*weights(:,1));
            ISurfTot(1,2) = sum(ISurf(:,2).*weights(:,1));
            ISurfTot(1,3) = sum(ISurf(:,3).*weights(:,1));
        end
        %disp('Stop');
    end
    IVolTotAll(fNo,:) = [IVolTot(1,1),IVolTot(1,2),IVolTot(1,3)];
    ISurfTotAll(fNo,:) = [ISurfTot(1,1),ISurfTot(1,2),ISurfTot(1,3)];
    ITotAll(fNo,:) = [IVolTot(1,1)-ISurfTot(1,1),IVolTot(1,2)-ISurfTot(1,2),...
        IVolTot(1,3)-ISurfTot(1,3)];
    k = 2*pi*f/vP;
    kr1 = k*r1; kr2 = k*r2;
    
    volNNTheo(fNo,1) = (cos(kr2)/(kr2)^2-sin(kr2)/(kr2)^3)-(cos(kr1)/(kr1)^2-sin(kr1)/(kr1)^3);
    
%     figure(1);
%     hold on;
%     plot(f,G*rho*abs(IVolTot(1,1)),'b*');
%     hold off;
end

figure(1);
subplot(1,3,1)
plot(fAll,G*rho*abs(IVolTotAll(:,1)),'b','LineWidth',2)
hold on;
plot(fAll,G*rho*abs(ISurfTotAll(:,1)),'r','LineWidth',2);
plot(fAll,G*rho*abs(ITotAll(:,1)),'m','LineWidth',2);
plot(fAll,8*pi/3*G*rho*ones(length(fAll),1),'k','LineWidth',2);
plot(fAll,G*rho*4*pi/3*ones(length(fAll),1),'g','LineWidth',2);
plot(fAll,8*pi*G*rho*volNNTheo,'ko','MarkerSize',6,'MarkerFaceColor','m');
legend({'Volume contribution','Surface contribution','Total','8\pi/3G\rho',...
    '4\pi/3G\rho','Theoretical'})
xlabel('Frequency (Hz)');
ylabel('NN ASD (m/s^2/\surd(Hz))');
title('X component');

subplot(1,3,2)
plot(fAll,G*rho*abs(IVolTotAll(:,2)),'b','LineWidth',2)
hold on;
plot(fAll,G*rho*abs(ISurfTotAll(:,2)),'r','LineWidth',2);
plot(fAll,G*rho*abs(ITotAll(:,2)),'m','LineWidth',2);
plot(fAll,8*pi/3*G*rho*ones(length(fAll),1),'k','LineWidth',2);
plot(fAll,G*rho*4*pi/3*ones(length(fAll),1),'g','LineWidth',2);
plot(fAll,8*pi*G*rho*volNNTheo,'ko','MarkerSize',6,'MarkerFaceColor','m');
legend({'Volume contribution','Surface contribution','Total','8\pi/3G\rho',...
    '4\pi/3G\rho','Theoretical'})
xlabel('Frequency (Hz)');
ylabel('NN ASD (m/s^2/\surd(Hz))');
title('Y component');

subplot(1,3,3)
plot(fAll,G*rho*abs(IVolTotAll(:,3)),'b','LineWidth',2)
hold on;
plot(fAll,G*rho*abs(ISurfTotAll(:,3)),'r','LineWidth',2);
plot(fAll,G*rho*abs(ITotAll(:,3)),'m','LineWidth',2);
plot(fAll,8*pi/3*G*rho*ones(length(fAll),1),'k','LineWidth',2);
plot(fAll,G*rho*4*pi/3*ones(length(fAll),1),'g','LineWidth',2);
plot(fAll,8*pi*G*rho*volNNTheo,'ko','MarkerSize',6,'MarkerFaceColor','m');
legend({'Volume contribution','Surface contribution','Total','8\pi/3G\rho',...
    '4\pi/3G\rho','Theoretical'});
xlabel('Frequency (Hz)');
ylabel('NN ASD (m/s^2/\surd(Hz))');
title('Z component');
