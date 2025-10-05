function [ampOut] = getSWaveAmp(x,y,z,vS,f,thetaW,phiW,wC)
% for a given y,z and x grids obtained using getXGrid we use this to
% get the Plane P-wave amplitude
% x, y, z are column vectors
% vS is the W-wave velocity in m/s, f is frequency in Hz
% thetaW, phiW are the directions of propagation of the wave
% wC = 'SV' or 'SH'
lambda = vS/f;
k = 2*pi/lambda;
kVec = [k*sind(thetaW)*cosd(phiW);k*sind(thetaW)*sind(phiW);k*cosd(thetaW)];
rVec = [x,y,z];
i = sqrt(-1);
kr = rVec*kVec;
ampOut = exp(-i*(kr));
if(wC=='SH')
    kUnitVec = [-sind(phiW);cosd(phiW);0];
end
if(wC=='SV')
    kUnitVec = [cosd(thetaW)*cosd(phiW);cosd(thetaW)*sind(phiW);-sind(thetaW)];
end
%kUnitVec = kVec./(sqrt(kVec(1,1)^2+kVec(2,1)^2+kVec(3,1)^2));
ampOut = [ampOut*kUnitVec(1,1),ampOut*kUnitVec(2,1),ampOut*kUnitVec(3,1)];

end