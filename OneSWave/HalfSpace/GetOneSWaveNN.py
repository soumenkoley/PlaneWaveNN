fAll = np.arange(2.0, 7.01, 0.2)     
nR = 300                             
R = 6000.0                           
nTheta = 100                         
nPhi = 500                          
rCav = 20.0                       
zCav = -250.0                     

thetaW = 90.0                       
phiW = 0.0                      
wComp = 'SV'                     

vS = 3000.0                       
rho = 2800.0                  
G = 6.67430e-11                  
r1 = rCav
r2 = R

bounds = [0, -zCav - 2*rCav, -zCav + 2*rCav, R]   
pointsPerRegion = [30, 30, 240]                 

rVec = []
wR = []
for i in range(3):
    r_sub, w_sub = lgwt(pointsPerRegion[i], bounds[i], bounds[i+1])
    rVec.append(r_sub)
    wR.append(w_sub)
rVec = np.concatenate(rVec)
wR = np.concatenate(wR)

theta, wTheta = lgwt(nTheta, np.pi/2, np.pi)
phi = np.linspace(0, 2*np.pi, nPhi+1)[:-1]
wPhi = 2*np.pi / nPhi

TH, PH = np.meshgrid(theta, phi)
thetaFlat = TH.ravel()
phiFlat = PH.ravel()

sinTheta = np.sin(thetaFlat)
cosTheta = np.cos(thetaFlat)
sinPhi = np.sin(phiFlat)
cosPhi = np.cos(phiFlat)

wtTheta = np.tile(wTheta, nPhi)
angularWeights = sinTheta * wtTheta * wPhi

thetaPhi, weights = createSurfPoints(np.max(rVec), nTheta, nPhi)
sThetaCPhi = np.sin(thetaPhi[:, 0]) * np.cos(thetaPhi[:, 1])
sThetaSPhi = np.sin(thetaPhi[:, 0]) * np.sin(thetaPhi[:, 1])
cTheta = np.cos(thetaPhi[:, 0])
dsUnitVec = np.column_stack((sThetaCPhi, sThetaSPhi, cTheta))

IVolTotAll = []
ISurfTotAll = []
ITotAll = []
volNNTheo = []


for f in fAll:
    IVolTot = np.zeros(3, dtype=complex)
    ISurfTot = np.zeros(3, dtype=complex)

    for i, (r, wr) in enumerate(zip(rVec, wR)):
        
        x = r * sinTheta * cosPhi
        y = r * sinTheta * sinPhi
        z = r * cosTheta

        
        insideCavity = (x**2 + y**2 + (z - zCav)**2) < rCav**2
        mask = ~insideCavity

        x_f = x[mask]
        y_f = y[mask]
        z_f = z[mask]
        w = angularWeights[mask]

        r_nodes = np.column_stack((x_f, y_f, z_f))
        IVol = getVolNNSWave(r_nodes, zCav, vS, f, thetaW, phiW, wComp)

        IVolTot += np.sum(IVol * (r**2 * w)[:, None], axis=0) * wr

        #  Surface integral 
        if i == len(rVec) - 1:
            dsUnitVecMask = dsUnitVec[mask, :]
            weightsMask = weights[mask]
            ampOut = getSWaveAmp(x_f, y_f, z_f, vS, f, thetaW, phiW, wComp)
            ISurf = getSurfNN(ampOut, x_f, y_f, z_f, zCav, dsUnitVecMask)

            ISurfTot = np.sum(ISurf * weightsMask[:, None], axis=0)

    IVolTotAll.append(IVolTot)
    ISurfTotAll.append(ISurfTot)
    ITotAll.append(IVolTot - ISurfTot)

    # Theoretical comparison 
    k = 2 * np.pi * f / vS
    kr1, kr2 = k * r1, k * r2
    volNNTheo.append((np.cos(kr2)/kr2**2 - np.sin(kr2)/kr2**3)
                     - (np.cos(kr1)/kr1**2 - np.sin(kr1)/kr1**3))

#plot
IVolTotAll = np.array(IVolTotAll)
ISurfTotAll = np.array(ISurfTotAll)
ITotAll = np.array(ITotAll)
volNNTheo = np.array(volNNTheo)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
labels = ['X component', 'Y component', 'Z component']

for i, ax in enumerate(axes):
    ax.plot(fAll, G * rho * np.abs(IVolTotAll[:, i]), 'b', label='Volume contribution', linewidth=2)
    ax.plot(fAll, G * rho * np.abs(ISurfTotAll[:, i]), 'r', label='Surface contribution', linewidth=2)
    ax.plot(fAll, G * rho * np.abs(ITotAll[:, i]), 'm', label='Total', linewidth=2)
    ax.plot(fAll, 8*np.pi/3 * G * rho * np.ones_like(fAll), 'k', label='8π/3 Gρ', linewidth=2)
    ax.plot(fAll, 4*np.pi/3 * G * rho * np.ones_like(fAll), 'g', label='4π/3 Gρ', linewidth=2)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('NN ASD (m/s²/√Hz)')
    ax.set_title(labels[i])
    ax.legend()

plt.tight_layout()
plt.show()
