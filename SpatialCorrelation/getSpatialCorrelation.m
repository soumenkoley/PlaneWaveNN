% this script performs the stiochastic field generation in half space
% it thyen ncomputes the correlation between a point on the top surface
% and a point at the center of the cavity
% the reason is to demonstrate that the spatial corelation between
% two points in the medium which is exciting by a set a plane waves
% converges to a Sin(kr)/krn type of result where is k = wavenumber
% and r is the distance between the points
% Lets say there was only one plane wave, in that case the coherence
% distance is infinite implying correlation between any two points is
% always 1, but that is not what happens in practice
% Spatial corre;lation from reaql field measuremnets show a Sin(kr)/kr
% behavior

clear; close all;

% no of realizations
nRea = 100; % increasing number of relaizations will increase the smoothness mof simulated result 
nWaves = 100; % number of plane waves per realization

% the two points in this case are (0,0,0) and (0,0,250)
x = [0;0]; y = [0;0]; z = [0;250];

vP = 4000; vS = vP*0.75; % P wave and S-wave velocity in m/s
xyz = [x(:), y(:), z(:)];  % [Npoints x 3]

% loop over these frequency values, so 1-20 Hz
fAll = 1:0.2:20; % in Hz

for fNo = 1:1:length(fAll)
    [freqs, gammaReal(fNo,1)] = computeSpatialCoherence(xyz, 1, 2, nRea, nWaves, fAll(fNo), vP, vS);
end

r = norm(xyz(1,:) - xyz(2,:));
kr = 2*pi*fAll*r./((vP+vS)/2);

figure(1);
hold on;
plot(fAll, gammaReal, 'b', 'DisplayName', 'Simulated coherence');
plot(fAll, sin(kr)./kr, 'r--', 'DisplayName', 'sin(kr)/kr');
legend;
xlabel('Frequency (Hz)'); ylabel('Real[Coherence]');
hold off;