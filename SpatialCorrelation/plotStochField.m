% plot ASD of stochastic field per realization
clear; close all;

nR = 100; % numver of points along radial vector
R = 4000; % units in meters, radius of the hemisphere
nTheta = 90; % number of points in theta for [pi/2,pi]
Nwaves = 100; % number of plane waves in one relaization
f = 3; % frequency in Hz
vP = 4000; % P-wave veloxcity in m/s
vS = vP*0.75; % S-wave velocity m/s
xC = 0; yC = 0; zC = -250; % location of cavity

[rVec, wR] = lgwt(nR, 0, R);

[theta, wTheta] = lgwt(nTheta, pi/2, pi);

phi = [0;pi];

% generate x,y,z
x = []; z = [];

for i = 1:1:length(phi)
    for j = 1:1:length(rVec)
        x = [x;rVec(j)*sin(theta)*cos(phi(i))];
        z = [z;rVec(j)*cos(theta)];
    end
end
y = zeros(length(x),1);
% note the above is just a xz slice at y = 0;

distToCav = sqrt(sum((x-xC).^2 + (y-yC).^2 + (z-zC).^2,2));
[minVal,minInd] = min(distToCav);

% then at all (x,y,z)
[uPAll, uSAll] = synthStochFieldInde([x,y,z], Nwaves, f, vP, vS);
uAll = uPAll+uSAll;
uAllASDX(:,1) = abs(uAll(:,1));
uAllASDY(:,1) = abs(uAll(:,2));
uAllASDZ(:,1) = abs(uAll(:,3));

% use the index to get the ASD of a point closest to cavity
uCavASDX(1,1) = abs(uAll(minInd,1));
uCavASDY(1,1) = abs(uAll(minInd,2));
uCavASDZ(1,1) = abs(uAll(minInd,3));

figure(1);
scatter(x, z, 40, uAllASDZ./uCavASDZ(1,1), 'filled');
xlabel('X-coordinate (m)');
ylabel('Z-coordinate (m)');
title(['ASD of Z-component of stochastic field at f = ',num2str(f),' Hz']);
colorbar;
caxis([min(uAllASDZ./uCavASDZ(1,1)),max(uAllASDZ./uCavASDZ(1,1))]);
% note that at some points the ASD is greater than 1, but thr idea is
% you have ASD =1 at the closest point to 250 meters depth
% the pattern will vary every time you run, that is why its stochastic
% end of code