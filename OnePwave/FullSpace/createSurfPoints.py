def createSurfPoints(R, Ntheta, Nphi):
   
    theta, wtheta = lgwt(Ntheta, pi/2, pi)
    
    phi = np.linspace(0, 2*pi, Nphi, endpoint=False)
    wphi = 2*pi / Nphi
    
    THETA, PHI = meshgrid(theta, phi)  # shape (Nphi, Ntheta)
    thetaPhi = column_stack((THETA.ravel(), PHI.ravel()))
    
    WTHETA, _ = meshgrid(wtheta, phi)
    SIN_THETA, _ = meshgrid(sin(theta), phi)
    
    weights = R**2 * wphi * (SIN_THETA * WTHETA)
    weights = weights.ravel()
    
    return thetaPhi, weights
