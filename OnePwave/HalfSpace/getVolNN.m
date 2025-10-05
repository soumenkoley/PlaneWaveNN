function [IFull] = getVolNN(rVec,zCav,vP,f,thetaW,phiW)
%this function computes the volume contribution of NN from all nodes
% defined by rVec
% vP is the P-wave velocity in m/s, f is frequency in Hz
% thetaW and phiW are the dierction of the plane wave
% first get the amplitudes at all node points on the plane
% note that getPWaveAmp is only used prior to getSurfNN
% you could modify getSurfNN like this to altogether remove getPWaveAmp
% in this code the amp generation happens within this function

lambda = vP/f;
k = 2*pi/lambda;
kVec = [k*sind(thetaW)*cosd(phiW);k*sind(thetaW)*sind(phiW);k*cosd(thetaW)];
%rVec = [x',y*ones(length(x),1),z*ones(length(x),1)];
i = sqrt(-1);
kr = rVec*kVec;
ampOut = exp(-i*(kr));
kUnitVec = kVec./(sqrt(kVec(1,1)^2+kVec(2,1)^2+kVec(3,1)^2));
ampOut = [ampOut*kUnitVec(1,1),ampOut*kUnitVec(2,1),ampOut*kUnitVec(3,1)];

% update the rVec for the current calculation
rVec(:,3) = rVec(:,3)-zCav; % (r-r0) vector
rDist = rVec.^2;
rDist = sum(rDist,2);
rDist = sqrt(rDist);
rCapX = rVec(:,1)./rDist;
rCapY = rVec(:,2)./rDist;
rCapZ = rVec(:,3)./rDist;
rDist3 = rDist.^3;

I1X = ampOut(:,1)./(rDist3);
I1Y = ampOut(:,2)./(rDist3);
I1Z = ampOut(:,3)./(rDist3);

% now get I2
I2 = (rCapX).*ampOut(:,1)+(rCapY).*ampOut(:,2)+(rCapZ).*ampOut(:,3);
I2 = 3*[I2.*(rCapX),I2.*(rCapY),I2.*(rCapZ)];
I2(:,1) = I2(:,1)./(rDist3);
I2(:,2) = I2(:,2)./(rDist3);
I2(:,3) = I2(:,3)./(rDist3);

IFull = [I1X-I2(:,1),I1Y-I2(:,2),I1Z-I2(:,3)];

end

