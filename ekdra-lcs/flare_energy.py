def Response_TESS():

    """
    Returns the TESS response curve as a function of the wavelength
    lambd. It's a spline interpolation 
    of the TESS response function.
    """

    responsefile = ""
    rs = np.genfromtxt(responsefile, delimiter=",")
    f = interp1d(rs[:, 0], rs[:, 1])
    return f


def B(lambd, T):

    """
    The Planck distribution as a function of wavelength
    """
    l = lambd * NANOMETER
    return (2.0 * C_PLANCK * C_C ** 2.0 / l ** 5.0) / (
        np.exp(C_PLANCK / (l * T * K_B)) - 1.0
    )


def Lp_STAR(R_star, T_star, R_TESS):

    """
    Incident luminosity of the star in the wavelength band
    """

    integrand = lambda lam: B(lam, T_star) * R_TESS(lam)
    integr = integrate.quad(integrand, LAMBDA_MIN, LAMBDA_MAX)

    return np.pi * R_STAR ** 2 * integr


def Lp_FLARE_DENSITY(T_flare, R_TESS):

    """
    Incident luminosity density of the flare in the TESS band, divided by its
    area
    """

    integrand = lambda lam: B(lam, T_star) * R_TESS(lam)
    integr = integrate.quad(integrand, LAMBDA_MIN, LAMBDA_MAX)

    return integr


def A_flare(rel_amp, LpS, LpF):

    """
    Returns the estimated effective area of the flare, given the
    relative amplitude from the normalized lightcurve;
    the luminosity of the star, and the luminosity of the flare over
    its area
    """

    return rel_amp * LpS / LpF


def L_flare(T_flare, A_f):

    """
    The luminosity of the flare at a given time:
    takes in the area of the flare and its temperature
    """

    return SIGMA_SB * T_flare ** 4.0 * A_f


def E_flare(amps, times, T_flare, LpS, LpF):

    """
    The energy of the flare. Takes in the amplitude,
    and pre-calculated values for the luminosity of the 
    star and the luminosity of the flare per unit area
    """

    dt = times[1] - times[0]
    A_f = A_flare(amps[0])
    integral = dt / 3.0 * L_flare(T_flare, A_f)

    for i in range(1, len(amps) - 1):
        A_f = A_flare(amps[i])
        if i % 2 == 0:
            integral += dt / 3.0 * L_flare(T_flare, A_f) * 2.0
        else:
            integral += dt / 3.0 * L_flare(T_flare, A_f) * 4.0

    A_f = A_flare(amps[-1])
    integral = dt / 3.0 * L_flare(T_flare, A_f)
    return integral


def interflare_time(flare_points, time):

    # Returns a series of interarrival times between flares

    if_times = []
    in_flare = False
    ctr = 0

    for i in range(len(time)):
        if flare_points[i] == 1:
            in_flare = True
        elif flare_points[i] == 0 and in_flare:
            in_flare = False
            ctr += 1
        else:
            ctr += 1
