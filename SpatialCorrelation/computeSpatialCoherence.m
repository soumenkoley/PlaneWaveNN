function [freqs, gammaReal] = computeSpatialCoherence(xyz, idx1, idx2, Nreal, nWaves, f, vP, vS)
% xyz: [Npts x 3]
% idx1, idx2: indices of two observation points
% Nreal: number of realizations
% f: frequency
% vP: wave speed (e.g., for P waves)
    
    omega = 2*pi*f;
    k = omega / vP;
    gammaSum = 0;
    norm1 = 0;
    norm2 = 0;

    for i = 1:Nreal
        % One realization: P waves only, or use total field
        [uP, uS] = synthStochFieldInde(xyz, nWaves, f, vP, vS);  % same as before

        u1 = uP(idx1,:)+uS(idx1,:);  % [1x3]
        u2 = uP(idx2,:)+uS(idx2,:);  % [1x3]

        gammaSum = gammaSum + u1 *(u2');  % dot product
        norm1 = norm1 + norm(u1)^2;
        norm2 = norm2 + norm(u2)^2;
    end

    % Normalize
    gamma = gammaSum / sqrt(norm1 * norm2);
    gammaReal = real(gamma);
    freqs = f;  % just a placeholder if you want to loop over frequencies
end