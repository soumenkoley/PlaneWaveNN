# Function to auto-compute grid resolution from R, vP, and frequency

def compgrid(R, vP, f_or_fmax,
                           samples_per_wavelength=6,
                           hemisphere=False,
                           min_nR=40, max_nR=300,
                           min_nTheta=20, max_nTheta=80,
                           min_nPhi=40, max_nPhi=200):
    lam = vP / float(f_or_fmax)                # smallest wavelength
    s = float(samples_per_wavelength)
    delta = lam / s                            # target physical step

    theta_span = 0.5 * np.pi if hemisphere else np.pi
    phi_span = 2.0 * np.pi

    nR_f = s * R / lam
    nTheta_f = s * theta_span * R / lam
    nPhi_f = s * phi_span * R / lam

    nR = int(min(max(min_nR, math.ceil(nR_f)), max_nR))
    nTheta = int(min(max(min_nTheta, math.ceil(nTheta_f)), max_nTheta))
    nPhi = int(min(max(min_nPhi, math.ceil(nPhi_f)), max_nPhi))

    info = {
        "lambda": lam,
        "target_step_m": delta,
        "nR_float": nR_f,
        "nTheta_float": nTheta_f,
        "nPhi_float": nPhi_f,
        "samples_per_wavelength": s
    }
    return nR, nTheta, nPhi, info
