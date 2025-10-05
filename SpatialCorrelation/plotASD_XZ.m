function plotASD_XZ(uTotal, xyz, f, y0)
% uTotal: [Npts x 3] complex field (P+S)
% xyz: [Npts x 3]
% f: frequency (Hz)
% y0: fixed y value (scalar)

    % Step 1: Select XZ plane at y = y0
    tol = 1e-6;
    sel = abs(xyz(:,2) - y0) < tol;
    xyz_sel = xyz(sel, :);         % [Nsel x 3]
    u_sel   = uTotal(sel, :);      % [Nsel x 3]

    % Step 2: Compute ASD at each point
    % ASD = sqrt(2) * |amplitude| (for one-sided PSD)
    asdX = sqrt(2) * abs(u_sel(:,1));
    asdZ = sqrt(2) * abs(u_sel(:,3));

    % Step 3: Reshape into XZ grid
    xVals = unique(xyz_sel(:,1));
    zVals = unique(xyz_sel(:,3));

    [X, Z] = meshgrid(xVals, zVals);
    Nx = length(xVals); Nz = length(zVals);

    % Convert vector to matrix form
    asdXgrid = reshape(asdX, [Nx, Nz])';
    asdZgrid = reshape(asdZ, [Nx, Nz])';

    % Step 4: Plot
    figure;
    subplot(1,2,1);
    imagesc(xVals, zVals, 10*log10(asdXgrid)); set(gca, 'YDir', 'normal');
    title(['ASD X-comp at f = ', num2str(f), ' Hz']);
    xlabel('x (m)'); ylabel('z (m)'); colorbar;
    
    subplot(1,2,2);
    imagesc(xVals, zVals, 10*log10(asdZgrid)); set(gca, 'YDir', 'normal');
    title(['ASD Z-comp at f = ', num2str(f), ' Hz']);
    xlabel('x (m)'); ylabel('z (m)'); colorbar;
end