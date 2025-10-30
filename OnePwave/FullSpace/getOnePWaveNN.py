#This script takes 1 minutes to run for largest grid.

fAll = np.arange(2, 7.01, 0.2)
R = 8000                      
rCav, zCav = 0.1, 0
thetaW, phiW = 90, 0
vP, rho = 4000, 2800
G = 6.67430e-11
r1, r2 = rCav, R

f_max = fAll.max()
nR, nTheta, nPhi, grid_info = compgrid(R, vP, f_max)
print(f"\n--- Grid resolution ---")
print(f"nR={nR}, nTheta={nTheta}, nPhi={nPhi}")
print(f"λ = {grid_info['lambda']:.2f} m, target step ≈ {grid_info['target_step_m']:.2f} m")

bounds = [0, rCav, R]
pointsPerRegion = [int(nR * 0.05), int(nR * 0.95)]  # 5% near cavity

#radial
rVec = np.concatenate([lgwt(pointsPerRegion[i], bounds[i], bounds[i+1])[0] for i in range(2)])
wR = np.concatenate([lgwt(pointsPerRegion[i], bounds[i], bounds[i+1])[1] for i in range(2)])

#angular
theta, wTheta = lgwt(nTheta, 0, np.pi)
phi = np.linspace(0, 2*np.pi, nPhi, endpoint=False)
wPhi = 2*np.pi / nPhi

TH, PH = np.meshgrid(theta, phi)
thetaFlat, phiFlat = TH.ravel(), PH.ravel()

sinTheta, cosTheta = np.sin(thetaFlat), np.cos(thetaFlat)
cosPhi, sinPhi = np.cos(phiFlat), np.sin(phiFlat)
unitDirs = np.column_stack((sinTheta * cosPhi, sinTheta * sinPhi, cosTheta))

wtTheta = np.tile(wTheta, nPhi)
angularWeights = sinTheta * wtTheta * wPhi

x_all = rVec[:, None] * unitDirs[:, 0]
y_all = rVec[:, None] * unitDirs[:, 1]
z_all = rVec[:, None] * unitDirs[:, 2]

cavity_mask = (x_all**2 + y_all**2 + (z_all - zCav)**2) < rCav**2

#surface points
thetaPhi, weights = createSurfPoints(max(rVec), nTheta, nPhi)
sThetaCPhi = np.sin(thetaPhi[:, 0]) * np.cos(thetaPhi[:, 1])
sThetaSPhi = np.sin(thetaPhi[:, 0]) * np.sin(thetaPhi[:, 1])
cTheta = np.cos(thetaPhi[:, 0])
dsUnitVec = np.column_stack((sThetaCPhi, sThetaSPhi, cTheta))

IVolTotAll, ISurfTotAll, ITotAll, volNNTheo = [], [], [], []

#loop over freq
for f in fAll:
    k = 2 * np.pi * f / vP  # wavenumber

    # Volume integration 
    mask = ~cavity_mask
    x, y, z = x_all[mask], y_all[mask], z_all[mask]

    # weights
    r_index, ang_index = np.divmod(np.flatnonzero(mask), unitDirs.shape[0])
    r_vals = rVec[r_index]
    wr_vals = wR[r_index]
    w_ang = angularWeights[ang_index]

    r_nodes = np.column_stack((x, y, z))
    IVol = getVolNN(r_nodes, zCav, vP, f, thetaW, phiW)

    full_weights = (r_vals**2 * w_ang * wr_vals)[:, None]
    IVolTot = np.sum(IVol * full_weights, axis=0)

    # Surface integration
    outer_mask = ~cavity_mask[-1]
    xS, yS, zS = x_all[-1, outer_mask], y_all[-1, outer_mask], z_all[-1, outer_mask]
    dsUnitVecMask = dsUnitVec[outer_mask]
    weightsMask = weights[outer_mask]

    ampOut = getPWaveAmp(xS, yS, zS, vP, f, thetaW, phiW)
    ISurf = getSurfNN(ampOut, xS, yS, zS, zCav, dsUnitVecMask)
    ISurfTot = np.sum(ISurf * weightsMask[:, None], axis=0)

    # Store results
    IVolTotAll.append(IVolTot)
    ISurfTotAll.append(ISurfTot)
    ITotAll.append(IVolTot - ISurfTot)

    # Theoretical
    kr1, kr2 = k * r1, k * r2
    volNNTheo.append(
        (np.cos(kr2)/kr2**2 - np.sin(kr2)/kr2**3)
        - (np.cos(kr1)/kr1**2 - np.sin(kr1)/kr1**3)
    )

# Convert to arrays
IVolTotAll = np.array(IVolTotAll)
ISurfTotAll = np.array(ISurfTotAll)
ITotAll = np.array(ITotAll)
volNNTheo = np.array(volNNTheo)

#plot
components = ['X', 'Y', 'Z']
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

for idx, ax in enumerate(axes):
    ax.plot(fAll, G*rho*np.abs(IVolTotAll[:, idx]), 'b', label='Volume contribution')
    ax.plot(fAll, G*rho*np.abs(ISurfTotAll[:, idx]), 'r', label='Surface contribution')
    ax.plot(fAll, G*rho*np.abs(ITotAll[:, idx]), 'm', label='Total')
    ax.plot(fAll, 8*np.pi/3*G*rho*np.ones_like(fAll), 'k--', label='8π/3 Gρ')
    ax.plot(fAll, 4*np.pi/3*G*rho*np.ones_like(fAll), 'g--', label='4π/3 Gρ')
    ax.plot(fAll, 8*np.pi*G*rho*volNNTheo, 'mo', markerfacecolor='none', label='Theoretical')

    ax.set_title(f'{components[idx]} Component')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('NN ASD (m/s²/√Hz)')
    ax.legend()
    ax.grid(True)

plt.tight_layout()
plt.show()
