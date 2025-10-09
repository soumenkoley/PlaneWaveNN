def getPWaveAmp(x, y, z, vP, f, thetaW, phiW):
    
    lambda_ = vP / f
    k = 2 * np.pi / lambda_

    kVec = np.array([
        k * np.sin(np.deg2rad(thetaW)) * np.cos(np.deg2rad(phiW)),
        k * np.sin(np.deg2rad(thetaW)) * np.sin(np.deg2rad(phiW)),
        k * np.cos(np.deg2rad(thetaW))
    ])

    rVec = np.column_stack((x, y, z))

    kr = rVec @ kVec
    amp = np.exp(-1j * kr)
#normalize
    kUnitVec = kVec / np.sqrt(np.sum(kVec**2))

    
    ampOut = np.column_stack((
        amp * kUnitVec[0],
        amp * kUnitVec[1],
        amp * kUnitVec[2]
    ))

    return ampOut
