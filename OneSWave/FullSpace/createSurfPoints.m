function [thetaPhi, weights] = createSurfPoints(R, Ntheta, Nphi)
    % Generate flattened list of [theta, phi] points and weights for
    % integrating over the surface of a hemisphere of radius R
    % theta ? [pi/2, pi], phi ? [0, 2pi)

    % Gauss�Legendre nodes and weights in ? ? [?/2, ?]
    [theta, wtheta] = lgwt(Ntheta, pi/2, pi);  % both are Ntheta x 1

    % Uniform quadrature for ? ? [0, 2?)
    phi = linspace(0, 2*pi, Nphi + 1); phi(end) = [];
    wphi = (2*pi) / Nphi;  % scalar

    % Create all (theta, phi) combinations
    [THETA, PHI] = meshgrid(theta, phi);         % Nphi � Ntheta
    thetaPhi = [THETA(:), PHI(:)];               % (Nphi*Ntheta) � 2

    % Evaluate weights: R� sin(?) w? w? at each grid point
    [WTHETA, ~] = meshgrid(wtheta, phi);         % Nphi � Ntheta
    [SIN_THETA, ~] = meshgrid(sin(theta), phi);  % Nphi � Ntheta

    weights = R^2 * wphi * (SIN_THETA .* WTHETA);  % Nphi � Ntheta
    weights = weights(:);                         % (Nphi*Ntheta) � 1
end