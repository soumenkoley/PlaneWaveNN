def getSurfNN(ampOut, x, y, z, zCav, dsUnitVec):
    """
    Compute Newtonian Noise contribution from the outer surface.
    
    rVec = np.column_stack((x, y, z - zCav))

    rDist = np.sqrt(np.sum(rVec**2, axis=1))
    rDist3 = rDist**3

    Ids = np.sum(ampOut * dsUnitVec, axis=1)
    Ids = Ids / rDist3

    I1X = Ids * rVec[:, 0]
    I1Y = Ids * rVec[:, 1]
    I1Z = Ids * rVec[:, 2]

    ISurf = np.column_stack((I1X, I1Y, I1Z))

    return ISurf
