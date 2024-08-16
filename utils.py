from uncertainties.umath import cos, sin

def galactic_velocity(dist, b, l, A, B, C, K, U, V, W):
    """
    Calculate the galactic velocity components based on Wang et al. 2021 (Eqs. 7-9).

    Parameters:
    -----------
    dist : float
        Distance to the star in kiloparsecs (kpc).
    b : float
        Galactic latitude in radians.
    l : float
        Galactic longitude in radians.
    A, B, C, K : float
        Oort constants in km/s/kpc.
    U, V, W : float
        Solar velocity components in the U (radial), V (tangential), and W (vertical) directions in km/s.

    Returns:
    --------
    V_r0 : float
        Radial velocity component in km/s.
    V_l0 : float
        Tangential velocity component along the Galactic longitude (l) in km/s.
    V_b0 : float
        Tangential velocity component along the Galactic latitude (b) in km/s.
    """
    
    # Precompute trigonometric functions
    cb, sb = cos(b), sin(b)
    cl, sl = cos(l), sin(l)
    c2l, s2l = cos(2 * l), sin(2 * l)
    
    # Calculate velocity components
    V_r0 = -U * cl * cb - V * sl * cb - W * sb + dist * (K + A * s2l + C * c2l) * cb * cb
    V_l0 = U * sl - V * cl + dist * (A * c2l + B - C * s2l) * cb
    V_b0 = U * cl * sb + V * sl * sb - W * cb - dist * (K + A * s2l + C * c2l) * sb * cb
    
    return V_r0, V_l0, V_b0


def velocity_to_propermotion(dist, V_l0, V_b0):
    """
    Convert tangential velocity components to proper motions.

    Parameters:
    -----------
    dist : float
        Distance to the star in kiloparsecs (kpc).
    V_l0 : float
        Tangential velocity component along the Galactic longitude (l) in km/s.
    V_b0 : float
        Tangential velocity component along the Galactic latitude (b) in km/s.

    Returns:
    --------
    mu_l0 : float
        Proper motion in the l direction in milliarcseconds per year (mas/yr).
    mu_b0 : float
        Proper motion in the b direction in milliarcseconds per year (mas/yr).
    """
    
    # Conversion factor from velocity to proper motion
    mu_l0 = 0.211 / dist * V_l0
    mu_b0 = 0.211 / dist * V_b0 
    
    return mu_l0, mu_b0


def oort_constants(model):
    """
    Retrieve the Oort constants and solar velocity components for different models.

    Parameters:
    -----------
    model : str
        Model identifier. Accepted values are:
        'C07' : Comerón & Pasquali 2007
        'B19' : Bobylev & Bajkova 2019 (OSCs)
        'W21' : Wang et al. 2021

    Returns:
    --------
    A, B, C, K : float
        Oort constants in km/s/kpc.
    U, V, W : float
        Solar velocity components in km/s.

    Raises:
    -------
    ValueError:
        If the provided model identifier is not recognized.
    """

    if model == 'C07':
        A, B, C, K = 12.5, -12.5, 0.0, 0.0
        U, V, W = 7.0, 14.0, 7.0
    elif model == 'B19':
        A, B, C, K = 16.40, -12.31, 0.0, 0.0
        U, V, W = 8.53, 11.22, 7.83
    elif model == 'W21':
        A, B, C, K = 16.31, -11.99, -3.10, -1.25
        U, V, W = 11.69, 10.16, 7.67
    else:
        raise ValueError(f"ERROR! Model '{model}' is not recognized. Choose from 'C07', 'B19', 'W21'.")

    return A, B, C, K, U, V, W
