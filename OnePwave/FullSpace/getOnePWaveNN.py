#This script takes several minutes to run, optimization is definitely still possible

fAll = np.arange(2, 7.01, 0.2) 
nR = 300
R = 8000
nTheta = 200
nPhi = 500
rCav = 20
zCav = 0

thetaW = 90
phiW = 0
vP = 4000
rho = 2800
G = 6.67430e-11
r1, r2 = rCav, R

bounds = [0, rCav, R]
pointsPerRegion = [10, 290]  

rVec = np.array([])
wR = np.array([])

for i in range(2):
    r_sub, w_sub = lgwt(pointsPerRegion[i], bounds[i], bounds[i+1])
    rVec = np.concatenate((rVec, r_sub))
    wR = np.concatenate((wR, w_sub))

theta, wTheta = lgwt(nTheta, 0, np.pi)
phi = np.linspace(0, 2*np.pi, nPhi, endpoint=False)
wPhi = 2*np.pi / nPhi

TH, PH = np.meshgrid(theta, phi)
thetaFlat = TH.ravel()
phiFlat = PH.ravel()

thetaPhi, weights = createSurfPoints(max(rVec), nTheta, nPhi)
sThetaCPhi = np.sin(thetaPhi[:,0]) * np.cos(thetaPhi[:,1])
sThetaSPhi = np.sin(thetaPhi[:,0]) * np.sin(thetaPhi[:,1])
cTheta = np.cos(thetaPhi[:,0])
dsUnitVec = np.column_stack((sThetaCPhi, sThetaSPhi, cTheta))

sinTheta = np.sin(thetaFlat)
cosTheta = np.cos(thetaFlat)
cosPhi = np.cos(phiFlat)
sinPhi = np.sin(phiFlat)

wtTheta = np.tile(wTheta, nPhi)
angularWeights = sinTheta * wtTheta * wPhi

IVolTotAll = []
ISurfTotAll = []
ITotAll = []
volNNTheo = []

for f in fAll:
    IVolTot = np.zeros(3, dtype=complex)
    for i, r in enumerate(rVec):
        wr = wR[i]

        x = r * sinTheta * cosPhi
        y = r * sinTheta * sinPhi
        z = r * cosTheta

        insideCavity = (x**2 + y**2 + (z - zCav)**2) < rCav**2
        x = x[~insideCavity]
        y = y[~insideCavity]
        z = z[~insideCavity]
        w = angularWeights[~insideCavity]
#volume
        r_nodes = np.column_stack((x, y, z))
        IVol = getVolNN(r_nodes, zCav, vP, f, thetaW, phiW)

        IVolTot += np.sum(IVol * (r**2 * w)[:, None], axis=0)

# Surface 
        if i == len(rVec)-1:
            ampOut = getPWaveAmp(x, y, z, vP, f, thetaW, phiW)
            ISurf = getSurfNN(ampOut, x, y, z, zCav, dsUnitVec)
            ISurfTot = np.sum(ISurf * weights[:, None], axis=0)

    IVolTotAll.append(IVolTot)
    ISurfTotAll.append(ISurfTot)
    ITotAll.append(IVolTot - ISurfTot)

    # Theory
    k = 2*np.pi*f / vP
    kr1, kr2 = k*r1, k*r2
    volNNTheo.append((np.cos(kr2)/kr2**2 - np.sin(kr2)/kr2**3) -
                     (np.cos(kr1)/kr1**2 - np.sin(kr1)/kr1**3))

IVolTotAll = np.array(IVolTotAll)
ISurfTotAll = np.array(ISurfTotAll)
ITotAll = np.array(ITotAll)
volNNTheo = np.array(volNNTheo)

# Plot
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
