import argparse
import numpy as np
import ast
from uncertainties import ufloat, umath
from uncertainties.umath import *

import astropy.coordinates as coord
import astropy.units as u

from utils import *


def main(source_name, const='W21'):
    filename = f"{source_name}.txt"
    try:
        with open(filename) as f:
            data = ast.literal_eval(f.read())
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return
    except (SyntaxError, ValueError):
        print(f"Error: File '{filename}' contains invalid data.")
        return

    # Gaia DR3 data
    source = coord.SkyCoord(
        ra=data['ra'] * u.degree,
        dec=data['dec'] * u.degree,
        distance=(data['parallax'] * u.mas).to(u.kpc, u.parallax()),
        pm_ra_cosdec=data['mu_ra_cosdec'] * u.mas / u.yr,
        pm_dec=data['mu_dec'] * u.mas / u.yr,
        radial_velocity=data['V_r'] * u.km / u.s,
        frame='icrs'
    )

    # Assign and convert data to uncertainties objects (adimensional)
    parallax = ufloat(data['parallax'], data['parallax_err'])
    ra_err = radians(data['ra_err'])
    dec_err = radians(data['dec_err'])
    pm_ra_err = data['mu_ra_cosdec_err']
    pm_dec_err = data['mu_dec_err']
    v_r_err = data['V_r_err']

    ra = ufloat(source.ra.radian, ra_err)
    dec = ufloat(source.dec.radian, dec_err)
    pm_ra = ufloat(source.pm_ra_cosdec.to('mas/yr').value, pm_ra_err)
    pm_dec = ufloat(source.pm_dec.to('mas/yr').value, pm_dec_err)
    v_r = ufloat(source.galactic.radial_velocity.value, v_r_err)

    # Define the transformation matrix between celestial to Galactic coordinates as in: 
    # https://gea.esac.esa.int/archive/documentation/GDR1/Data_processing/chap_cu3ast/sec_cu3ast_intro.html
   
    # Matrix A_G' from Eq. 3.11:
    Agt = np.array([[-0.0548755604162154,-0.8734370902348850,-0.4838350155487132],\
                  [0.4941094278755837,-0.4448296299600112,0.7469822444972189],\
                  [-0.8676661490190047,-0.1980763734312015,0.4559837761750669]])

    dist = 1./parallax

    # Define the position vectors in ICRS and convert to Galactic:
    cra, sra = cos(ra), sin(ra)
    cd, sd = cos(dec), sin(dec)
    r_icrs = np.array([cra * cd, sra * cd, sd])   # Eq. 3.7
    r_gal = np.dot(Agt, r_icrs)                   # Eq. 3.9

    # Galactic longitude (l) and latitude (b), from Eq. 3.13
    l_gal = umath.atan2(r_gal[1], r_gal[0])   
    b_gal = umath.atan2(r_gal[2], sqrt(r_gal[0]**2 + r_gal[1]**2))
    cl, sl = cos(l_gal), sin(l_gal)
    cb, sb = cos(b_gal), sin(b_gal)

    # Auxiliary column matrices from Eqs. 3.14 and 3.15
    p_icrs = np.array([-sra, cra, 0])
    q_icrs = np.array([-cra * sd, -sra * sd, cd])
    p_gal = np.array([-sl, cl, 0])
    q_gal = np.array([-cl * sb, -sl * sb, cb])

    # Transform the proper motions to Galactic coordinates
    mu_icrs = p_icrs * pm_ra + q_icrs * pm_dec                     # Eq. 3.16
    mu_gal = np.dot(Agt, mu_icrs)                                  # Eq. 3.18
    mul, mub = np.dot(p_gal.T, mu_gal), np.dot(q_gal.T, mu_gal)    # Eq. 3.20 

    print(f'D = {dist:.1u} kpc')
    print(f'(l, b) = ({l_gal*180/np.pi:.1u} deg, {b_gal*180/np.pi:.1u} deg)')
    print(f'\nObserved proper motion: \u03BC_l = {mul:.1u} mas/yr, \u03BC_b = {mub:.1u} mas/yr')

    A, B, C, K, U, V, W = oort_constants(const)

    # Calculate Galactic rotation velocity components
    v_r0, v_l0, v_b0 = galactic_velocity(dist, b_gal, l_gal, A, B, C, K, U, V, W)

    # Convert velocities to proper motions
    mu_l0, mu_b0 = velocity_to_propermotion(dist, v_l0, v_b0)

    # Proper motions w.r.t. the surrounding medium:
    mu_l_corr = mul - mu_l0
    mu_b_corr = mub - mu_b0

    print(f'Corrected proper motion: \u03BC_l = {mu_l_corr:.1u} mas/yr, \u03BC_b = {mu_b_corr:.1u} mas/yr')

    #Calculate mu_ra, mu_dec, and V_tan. We use the transpose of the transformation matrix. 
    Ag = Agt.T

    # Transform the proper motions to Galactic coordinates
    mu_gal_corr = p_gal * mu_l_corr + q_gal * mu_b_corr                                       # Eq. 3.17
    mu_icrs_corr = np.dot(Ag, mu_gal_corr)                                                    # Eq. 3.19
    mu_ra_corr, mu_dec_corr = np.dot(p_icrs.T, mu_icrs_corr), np.dot(q_icrs.T, mu_icrs_corr)  # Eq. 3.21 

    print(f'\nObserved proper motion: \u03BC_ra = {pm_ra:.1u} mas/yr, \u03BC_dec = {pm_dec:.1u} mas/yr')
    print(f'Corrected proper motion: \u03BC_ra = {mu_ra_corr:.1u} mas/yr, \u03BC_dec = {mu_dec_corr:.1u} mas/yr')

    # Calculate proper motion angle before and after correction for Galactic rotation 
    angle = degrees( atan(pm_dec/pm_ra) )
    angle_corr = degrees( atan(mu_dec_corr/mu_ra_corr) )

    print(f'\nObserved angle = {angle:.1u}°')
    print(f'Corrected angle = {angle_corr:.1u}°')

    # Calculate the corrected tangential velocity
    v_t = 4.74 * dist * sqrt(mu_l_corr**2 + mu_b_corr**2)
    print(f'\nV_t = {v_t:.2u} km/s')

    # Calculate corrected radial velocity
    v_r_corr = v_r - v_r0
    print(f'V_r = {v_r_corr:.2u} km/s')

    print(f'V = {sqrt(v_r_corr**2 + v_t**2):.2u} km/s')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transform the proper motion of the source to its value with respect to its surrounding medium.")
    
    parser.add_argument("source", type=str, help="The source name, e.g., 'BD+43_3654'")
    parser.add_argument(
        "const", 
        type=str, 
        nargs='?', 
        default='W21', 
        help="The Oort constants to use (available: 'C07', 'B19', 'W21'; default is 'W21')"
    )  
    args = parser.parse_args()
    
    # Call the main function with the parsed arguments
    main(args.source, args.const)