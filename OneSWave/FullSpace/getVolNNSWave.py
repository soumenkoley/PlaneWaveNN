def getVolNNSwave(r_vec, z_cav, v_s, f, theta_w, phi_w, wC):
    
    theta_r = np.deg2rad(theta_w)
    phi_r = np.deg2rad(phi_w)

    k = 2 * np.pi * f / v_s
    k_vec = np.array([
        np.sin(theta_r) * np.cos(phi_r),
        np.sin(theta_r) * np.sin(phi_r),
        np.cos(theta_r)
    ]) * k

    kr = r_vec @ k_vec
    amp = np.exp(-1j * kr)

    if wC.upper() == 'SH':
        k_unit = np.array([-np.sin(phi_r), np.cos(phi_r), 0.0])
    elif wC.upper() == 'SV':
        k_unit = np.array([
            np.cos(theta_r) * np.cos(phi_r),
            np.cos(theta_r) * np.sin(phi_r),
            -np.sin(theta_r)
        ])
    else:
        raise ValueError("wC must be 'SH' or 'SV'")

    amp_vec = amp[:, None] * k_unit

    r = r_vec.copy()
    r[:, 2] -= z_cav
    r_dist2 = np.sum(r**2, axis=1)
    r_dist = np.sqrt(r_dist2)
    inv_r3 = 1.0 / (r_dist2 * r_dist)

    r_cap = r / r_dist[:, None]

    I1 = amp_vec * inv_r3[:, None]
    I2_dot = np.einsum('ij,ij->i', r_cap, amp_vec)
    I2 = 3 * I2_dot[:, None] * r_cap * inv_r3[:, None]

    return I1 - I2
