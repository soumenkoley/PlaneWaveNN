function [ampOut] = getPWaveAmp(x,y,z,vP,f,thetaW,phiW)
% for a given y,z and x grids
% get the Plane P-wave amplitude
% x, y, z are column vectors, z are scalar values
% vP is the P-wave velocity in m/s, f is frequency in Hz
% thetaW and phiW are the dierction of the plane wave
lambda = vP/f;
k = 2*pi/lambda;
kVec = [k*sind(thetaW)*cosd(phiW);k*sind(thetaW)*sind(phiW);k*cosd(thetaW)];
rVec = [x,y,z];
i = sqrt(-1);
kr = rVec*kVec;
ampOut = exp(-i*(kr));
kUnitVec = kVec./(sqrt(kVec(1,1)^2+kVec(2,1)^2+kVec(3,1)^2));
ampOut = [ampOut*kUnitVec(1,1),ampOut*kUnitVec(2,1),ampOut*kUnitVec(3,1)];

end