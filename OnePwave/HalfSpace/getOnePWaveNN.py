fAll = np.arange(2, 7.01, 0.2)
nTheta = 100
nPhi = 500
rCav = 20
zCav = -250
vP = 4000
rho = 2800
G = 6.67430e-11
r1, r2 = rCav, 8000


bounds = [0, -zCav - 2*rCav, -zCav + 2*rCav, r2]
pointsPerRegion = [30, 30, 240]

rVec = np.array([])
wR = np.array([])
for i in range(3):
    r_sub, w_sub = lgwt(pointsPerRegion[i], bounds[i], bounds[i+1])
    rVec = np.concatenate((rVec, r_sub))
    wR = np.concatenate((wR, w_sub))

nR = len(rVec)
nFreq = len(fAll)

theta, wTheta = lgwt(nTheta, np.pi/2, np.pi)
phi = np.linspace(0, 2*np.pi, nPhi, endpoint=False)
wPhi = 2*np.pi / nPhi

TH, PH = np.meshgrid(theta, phi)
thetaFlat = TH.ravel()
phiFlat = PH.ravel()
nPoints = len(thetaFlat)

sinTheta = np.sin(thetaFlat)
cosTheta = np.cos(thetaFlat)
cosPhi = np.cos(phiFlat)
sinPhi = np.sin(phiFlat)
wtTheta = np.tile(wTheta, nPhi)
angularWeights = sinTheta * wtTheta * wPhi

thetaPhi, weights = createSurfPoints(max(rVec), nTheta, nPhi)
sThetaCPhi = np.sin(thetaPhi[:,0]) * np.cos(thetaPhi[:,1])
sThetaSPhi = np.sin(thetaPhi[:,0]) * np.sin(thetaPhi[:,1])
cTheta = np.cos(thetaPhi[:,0])
dsUnitVec = np.column_stack((sThetaCPhi, sThetaSPhi, cTheta))

x_all = rVec[:, None] * sinTheta * cosPhi
y_all = rVec[:, None] * sinTheta * sinPhi
z_all = rVec[:, None] * cosTheta

cavity_mask = (x_all**2 + y_all**2 + (z_all - zCav)**2) < rCav**2

angularWeights = angularWeights[None, :]  # shape 1 x nPoints
weights_surface = weights[:, None]        # shape nPoints x 1

IVolTotAll = np.zeros((nFreq, 3), dtype=complex)
ISurfTotAll = np.zeros((nFreq, 3), dtype=complex)
ITotAll = np.zeros((nFreq, 3), dtype=complex)
volNNTheo = np.zeros(nFreq)

#Now loop over all vectorized frequencies

for f_idx, f in enumerate(fAll):
    lambda_f = vP / f
    k = 2 * np.pi / lambda_f
    kVec = k * np.array([np.sin(np.deg2rad(90)) * np.cos(np.deg2rad(0)),
                         np.sin(np.deg2rad(90)) * np.sin(np.deg2rad(0)),
                         np.cos(np.deg2rad(90))])
    kUnitVec = kVec / np.linalg.norm(kVec)

    IVolTot = np.zeros(3, dtype=complex)

    # and vectorized over all radii
    for i in range(nR):
        mask = ~cavity_mask[i,:]
        x = x_all[i, mask]
        y = y_all[i, mask]
        z = z_all[i, mask]
        r_nodes = np.column_stack((x, y, z))
        w = angularWeights[0, mask]
        r = rVec[i]
        wr = wR[i]

        # Volume vectorized
        IVol = getVolNN(r_nodes, zCav, vP, f, 90, 0)
        IVolTot += np.sum(IVol * (r**2 * w)[:, None], axis=0) * wr

        # Surface integral 
        if i == nR-1:
            ampOut = getPWaveAmp(x, y, z, vP, f, 90, 0)
            ISurf = getSurfNN(ampOut, x, y, z, zCav, dsUnitVec)
            ISurfTot = np.sum(ISurf * weights_surface, axis=0)

    IVolTotAll[f_idx] = IVolTot
    ISurfTotAll[f_idx] = ISurfTot
    ITotAll[f_idx] = IVolTot - ISurfTot

    # Theoretical formula
    kr1, kr2 = k*r1, k*r2
    volNNTheo[f_idx] = (np.cos(kr2)/kr2**2 - np.sin(kr2)/kr2**3) - \
                       (np.cos(kr1)/kr1**2 - np.sin(kr1)/kr1**3)

# --- Plot results ---
components = ['X', 'Y', 'Z']
for idx in range(3):
    plt.subplot(1,3,idx+1)
    plt.plot(fAll, G*rho*np.abs(IVolTotAll[:,idx]), 'b', label='Volume contribution')
    plt.plot(fAll, G*rho*np.abs(ISurfTotAll[:,idx]), 'r', label='Surface contribution')
    plt.plot(fAll, G*rho*np.abs(ITotAll[:,idx]), 'm', label='Total')
    plt.plot(fAll, 8*np.pi/3*G*rho*np.ones_like(fAll), 'k', label='8π/3Gρ')
    plt.plot(fAll, 4*np.pi/3*G*rho*np.ones_like(fAll), 'g', label='4π/3Gρ')
    plt.plot(fAll, 8*np.pi*G*rho*volNNTheo, 'ko', markerfacecolor='m', label='Theoretical')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('NN ASD (m/s^2/√Hz)')
    plt.title(f'{components[idx]} component')
    plt.legend()
    plt.grid(True)

plt.tight_layout()
plt.show()
