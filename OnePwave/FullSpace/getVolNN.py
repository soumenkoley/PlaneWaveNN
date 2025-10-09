def getVolNN(rVec, zCav, vP, f, thetaW, phiW):
    
    #Compute the volume contribution of Newtonian Noise (NN).

    lambda_ = vP / f
    k = 2 * np.pi / lambda_

    kVec = np.array([
        k * np.sin(np.deg2rad(thetaW)) * np.cos(np.deg2rad(phiW)),
        k * np.sin(np.deg2rad(thetaW)) * np.sin(np.deg2rad(phiW)),
        k * np.cos(np.deg2rad(thetaW))
    ])

    kr = rVec @ kVec
    ampOut = np.exp(-1j * kr)

    kUnitVec = kVec / np.linalg.norm(kVec)
    ampOut = np.column_stack((ampOut * kUnitVec[0],
                              ampOut * kUnitVec[1],
                              ampOut * kUnitVec[2]))

    rVec[:, 2] -= zCav

    rDist = np.sqrt(np.sum(rVec**2, axis=1))
    rCapX = rVec[:, 0] / rDist
    rCapY = rVec[:, 1] / rDist
    rCapZ = rVec[:, 2] / rDist
    rDist3 = rDist**3

    I1X = ampOut[:, 0] / rDist3
    I1Y = ampOut[:, 1] / rDist3
    I1Z = ampOut[:, 2] / rDist3

    I2_scalar = rCapX * ampOut[:, 0] + rCapY * ampOut[:, 1] + rCapZ * ampOut[:, 2]
    I2 = 3 * np.column_stack((I2_scalar*rCapX,
                              I2_scalar*rCapY,
                              I2_scalar*rCapZ))
    I2[:, 0] /= rDist3
    I2[:, 1] /= rDist3
    I2[:, 2] /= rDist3

    IFull = np.column_stack((I1X - I2[:, 0],
                             I1Y - I2[:, 1],
                             I1Z - I2[:, 2]))
    return IFull
