fAll = np.arange(2, 7.01, 0.2)  # Hz
R = 8000                        # outer radius (m)
rCav, zCav = 20, -250           # cavity position
thetaW, phiW = 90, 0
wComp = 'SV'                     # polarization: 'SV' or 'SH'
vS, rho = 3000, 2800
G = 6.67430e-11
r1, r2 = rCav, R


f_max = fAll.max()
nR, nTheta, nPhi, grid_info = compgrid(R, vS, f_max)
print(f"\n--- Grid resolution (Half-space S-wave) ---")
print(f"nR={nR}, nTheta={nTheta}, nPhi={nPhi}")
print(f"λ = {grid_info['lambda']:.2f} m, target step ≈ {grid_info['target_step_m']:.2f} m")


bounds = [0, -zCav - 2*rCav, -zCav + 2*rCav, R]
pointsPerRegion = [int(nR*0.1), int(nR*0.1), int(nR*0.8)]  # roughly 10%, 10%, 80%

r_segments, w_segments = zip(*[lgwt(n, a, b) for n, a, b in zip(pointsPerRegion, bounds[:-1], bounds[1:])])
rVec = np.concatenate(r_segments)
wR = np.concatenate(w_segments)


theta, wTheta = lgwt(nTheta, np.pi/2, np.pi)
phi = np.linspace(0, 2*np.pi, nPhi, endpoint=False)
wPhi = 2*np.pi / nPhi

TH, PH = np.meshgrid(theta, phi)
thetaFlat, phiFlat = TH.ravel(), PH.ravel()

sinTheta, cosTheta = np.sin(thetaFlat), np.cos(thetaFlat)
cosPhi, sinPhi = np.cos(phiFlat), np.sin(phiFlat)
unitDirs = np.column_stack([sinTheta*cosPhi, sinTheta*sinPhi, cosTheta])

wtTheta = np.tile(wTheta, nPhi)
angularWeights = sinTheta * wtTheta * wPhi

x_all = rVec[:, None] * unitDirs[:,0]
y_all = rVec[:, None] * unitDirs[:,1]
z_all = rVec[:, None] * unitDirs[:,2]

cavity_mask = (x_all**2 + y_all**2 + (z_all - zCav)**2) < rCav**2

# Surface points
thetaPhi, weights = createSurfPoints(max(rVec), nTheta, nPhi)
sThetaCPhi = np.sin(thetaPhi[:,0]) * np.cos(thetaPhi[:,1])
sThetaSPhi = np.sin(thetaPhi[:,0]) * np.sin(thetaPhi[:,1])
cTheta = np.cos(thetaPhi[:,0])
dsUnitVec = np.column_stack([sThetaCPhi, sThetaSPhi, cTheta])

IVolTotAll, ISurfTotAll, ITotAll, volNNTheo = [], [], [], []

for f in fAll:
    k = 2*np.pi*f / vS

    # Volume integration
    mask = ~cavity_mask
    x, y, z = x_all[mask], y_all[mask], z_all[mask]

    r_index, ang_index = np.divmod(np.flatnonzero(mask), unitDirs.shape[0])
    r_vals = rVec[r_index]
    wr_vals = wR[r_index]
    w_ang = angularWeights[ang_index]

    r_nodes = np.column_stack([x, y, z])
    IVol = getVolNNSWave(r_nodes, zCav, vS, f, thetaW, phiW, wComp)
    full_weights = (r_vals**2 * w_ang * wr_vals)[:, None]
    IVolTot = np.sum(IVol * full_weights, axis=0)

    # Surface integration
    outer_mask = ~cavity_mask[-1]
    xS, yS, zS = x_all[-1, outer_mask], y_all[-1, outer_mask], z_all[-1, outer_mask]
    dsUnitVecMask = dsUnitVec[outer_mask]
    weightsMask = weights[outer_mask]

    ampOut = getSWaveAmp(xS, yS, zS, vS, f, thetaW, phiW, wComp)
    ISurf = getSurfNN(ampOut, xS, yS, zS, zCav, dsUnitVecMask)
    ISurfTot = np.sum(ISurf * weightsMask[:,None], axis=0)

    # Store results
    IVolTotAll.append(IVolTot)
    ISurfTotAll.append(ISurfTot)
    ITotAll.append(IVolTot - ISurfTot)

    # Theoretical
    kr1, kr2 = k*r1, k*r2
    volNNTheo.append(
        (np.cos(kr2)/kr2**2 - np.sin(kr2)/kr2**3) -
        (np.cos(kr1)/kr1**2 - np.sin(kr1)/kr1**3)
    )

# Convert to arrays
IVolTotAll = np.array(IVolTotAll)
ISurfTotAll = np.array(ISurfTotAll)
ITotAll = np.array(ITotAll)
volNNTheo = np.array(volNNTheo)

#plot
components = ['X', 'Y', 'Z']
fig, axes = plt.subplots(1,3,figsize=(16,5))

for idx, ax in enumerate(axes):
    ax.plot(fAll, G*rho*np.abs(IVolTotAll[:,idx]), 'b', label='Volume contribution')
    ax.plot(fAll, G*rho*np.abs(ISurfTotAll[:,idx]), 'r', label='Surface contribution')
    ax.plot(fAll, G*rho*np.abs(ITotAll[:,idx]), 'm', label='Total')
    ax.plot(fAll, 8*np.pi/3*G*rho*np.ones_like(fAll), 'k--', label='8π/3 Gρ')
    ax.plot(fAll, 4*np.pi/3*G*rho*np.ones_like(fAll), 'g--', label='4π/3 Gρ')
    ax.plot(fAll, 8*np.pi*G*rho*volNNTheo, 'mo', markerfacecolor='none', label='Theoretical')

    ax.set_title(f'{components[idx]} Component (Half-space S-wave)')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('NN ASD (m/s²/√Hz)')
    ax.legend()
    ax.grid(True)

plt.tight_layout()
plt.show()
