def getSWaveAmp(x, y, z, vS, f, thetaW, phiW, wC):
   
    x, y, z = np.asarray(x), np.asarray(y), np.asarray(z)
    rVec = np.column_stack((x, y, z))

    lambda_ = vS / f
    k = 2 * np.pi / lambda_
    kx = k * np.sin(np.deg2rad(thetaW)) * np.cos(np.deg2rad(phiW))
    ky = k * np.sin(np.deg2rad(thetaW)) * np.sin(np.deg2rad(phiW))
    kz = k * np.cos(np.deg2rad(thetaW))
    kVec = np.array([kx, ky, kz])

    kr = rVec @ kVec
    amp_phase = np.exp(-1j * kr)

    if wC.upper() == 'SH':
        kUnitVec = np.array([-np.sin(np.deg2rad(phiW)),
                              np.cos(np.deg2rad(phiW)),
                              0.0])
    elif wC.upper() == 'SV':
        kUnitVec = np.array([ np.cos(np.deg2rad(thetaW)) * np.cos(np.deg2rad(phiW)),
                              np.cos(np.deg2rad(thetaW)) * np.sin(np.deg2rad(phiW)),
                             -np.sin(np.deg2rad(thetaW))])
    else:
        raise ValueError("wC must be 'SH' or 'SV'")

    ampOut = amp_phase[:, None] * kUnitVec[None, :]

    return ampOut
