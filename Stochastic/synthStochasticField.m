function [uP, uS] = synthStochasticField(xyz, Nwaves, f, vP, vS,theta,phi,phase)
% xyz: [Npoints x 3] array of positions (x, y, z)
% Nwaves: number of plane waves (per type)
% f: frequency in Hz
% vP, vS: P and S wave speeds
% Returns:
%   uP, uS: [Npoints x 3] complex displacement field vectors from P and S waves

% Constants
omega = 2*pi*f;
kP = omega / vP;
kS = omega / vS;
Npts = size(xyz,1);


% Direction unit vectors (Nwaves x 3)
khat = [sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];

% Random phases
%phase = 2*pi*rand(Nwaves,1);

%% -------- P-Waves -------- %%
% Displacement is parallel to propagation direction
uP = zeros(Npts, 3);

for i = 1:Nwaves
    kvec = kP * khat(i,:);              % [1x3]
    arg  = xyz * kvec.' + phase(i);     % [Npts x 1]
    amp  = exp(-1i * arg);              % [Npts x 1]
    uP   = uP + amp .* khat(i,:);       % broadcast to 3 cols
end

uP = uP / sqrt(Nwaves);  % Normalize for energy consistency

%% -------- S-Waves -------- %%
% Displacement is perpendicular to propagation direction
% We use SH waves only for simplicity (you can mix SV too)

uS = zeros(Npts, 3);
for i = 1:Nwaves
    % Pick a perpendicular unit vector (simple method)
    n = khat(i,:);  % propagation direction
    % Get a perpendicular unit vector (robust cross-product trick)
    if abs(n(3)) < 0.99
        perp = cross(n, [0,0,1]);
    else
        perp = cross(n, [0,1,0]);
    end
    perp = perp / norm(perp);  % normalize
    
    kvec = kS * n;
    arg  = xyz * kvec.' + phase(i);     % [Npts x 1]
    amp  = exp(-1i * arg);              % [Npts x 1]
    uS   = uS + amp .* perp;            % add componentwise
end

uS = uS / sqrt(Nwaves);  % Normalize

end