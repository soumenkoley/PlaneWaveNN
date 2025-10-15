fAll = np.arange(2.0, 7.01, 0.2)
nR = 300
R = 6000
nTheta = 200
nPhi = 500
rCav = 20
zCav = 0
thetaW, phiW = 90, 0
wComp = 'SV'
vS = 3000
rho = 2800
G = 6.67430e-11
r1, r2 = rCav, R

bounds = [0, rCav, R]
pointsPerRegion = [10, 290]

rVec = np.concatenate([lgwt(pointsPerRegion[i], bounds[i], bounds[i+1])[0] for i in range(2)])
wR = np.concatenate([lgwt(pointsPerRegion[i], bounds[i], bounds[i+1])[1] for i in range(2)])

theta, wTheta = lgwt(nTheta, 0, np.pi)
phi = np.linspace(0, 2*np.pi, nPhi, endpoint=False)
wPhi = 2*np.pi / nPhi

TH, PH = np.meshgrid(theta, phi)
thetaFlat, phiFlat = TH.ravel(), PH.ravel()

sinTheta, cosTheta = np.sin(thetaFlat), np.cos(thetaFlat)
cosPhi, sinPhi = np.cos(phiFlat), np.sin(phiFlat)
unitDirs = np.column_stack([sinTheta*cosPhi, sinTheta*sinPhi, cosTheta])

wtTheta = np.tile(wTheta, nPhi)
angularWeights = sinTheta * wtTheta * wPhi

x_all = rVec[:, None] * unitDirs[:, 0][None, :]
y_all = rVec[:, None] * unitDirs[:, 1][None, :]
z_all = rVec[:, None] * unitDirs[:, 2][None, :]

cavity_mask = (x_all**2 + y_all**2 + (z_all - zCav)**2) < rCav**2

thetaPhi, weights = createSurfPoints(max(rVec), nTheta, nPhi)
sThetaCPhi = np.sin(thetaPhi[:,0]) * np.cos(thetaPhi[:,1])
sThetaSPhi = np.sin(thetaPhi[:,0]) * np.sin(thetaPhi[:,1])
cTheta = np.cos(thetaPhi[:,0])
dsUnitVec = np.column_stack([sThetaCPhi, sThetaSPhi, cTheta])

IVolTotAll = []
ISurfTotAll = []
ITotAll = []
volNNTheo = []

for f in fAll:
    IVolTot = np.zeros(3, dtype=complex)
    
    for i, r in enumerate(rVec):
        # Mask cavity
        mask = ~cavity_mask[i]
        x, y, z = x_all[i, mask], y_all[i, mask], z_all[i, mask]
        w = angularWeights[mask]
        
        # Volume contribution
        r_nodes = np.column_stack((x, y, z))
        IVol = getVolNNSWave(r_nodes, zCav, vS, f, thetaW, phiW, wComp)
        IVolTot += np.sum(IVol * (r**2 * w)[:, None], axis=0) * wR[i]
    
    # Surface integral on outer shell
    outer_mask = ~cavity_mask[-1]
    x, y, z = x_all[-1, outer_mask], y_all[-1, outer_mask], z_all[-1, outer_mask]
    dsUnitVecMask = dsUnitVec[outer_mask]
    weightsMask = weights[outer_mask]
    ampOut = getSWaveAmp(x, y, z, vS, f, thetaW, phiW, wComp)
    ISurfTot = np.sum(getSurfNN(ampOut, x, y, z, zCav, dsUnitVecMask) * weightsMask[:, None], axis=0)
    
    IVolTotAll.append(IVolTot)
    ISurfTotAll.append(ISurfTot)
    ITotAll.append(IVolTot - ISurfTot)

    # Theoretical volume NN term (for comparison)
    k = 2 * np.pi * f / vS
    kr1, kr2 = k * r1, k * r2
    volNNTheo.append((np.cos(kr2)/kr2**2 - np.sin(kr2)/kr2**3)
                     - (np.cos(kr1)/kr1**2 - np.sin(kr1)/kr1**3))

# --- Convert lists to arrays ---
IVolTotAll = np.array(IVolTotAll)
ISurfTotAll = np.array(ISurfTotAll)
ITotAll = np.array(ITotAll)
volNNTheo = np.array(volNNTheo)

# --- Plotting ---
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
labels = ['X component', 'Y component', 'Z component']

for i, ax in enumerate(axes):
    ax.plot(fAll, G * rho * np.abs(IVolTotAll[:, i]), 'b', label='Volume', linewidth=2)
    ax.plot(fAll, G * rho * np.abs(ISurfTotAll[:, i]), 'r', label='Surface', linewidth=2)
    ax.plot(fAll, G * rho * np.abs(ITotAll[:, i]), 'm', label='Total', linewidth=2)
    ax.plot(fAll, 8 * np.pi / 3 * G * rho * np.ones_like(fAll), 'k', label='8π/3 Gρ', linewidth=2)
    ax.plot(fAll, 4 * np.pi / 3 * G * rho * np.ones_like(fAll), 'g', label='4π/3 Gρ', linewidth=2)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('NN ASD (m/s²/√Hz)')
    ax.set_title(labels[i])
    ax.legend()

plt.tight_layout()
plt.show()
