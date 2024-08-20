import argparse
import numpy as np
import matplotlib.pyplot as plt
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
    Ag = np.array([[-0.0548755604162154,-0.8734370902348850,-0.4838350155487132],\
                  [0.4941094278755837,-0.4448296299600112,0.7469822444972189],\
                  [-0.8676661490190047,-0.1980763734312015,0.4559837761750669]])

    dist = 1./parallax

    # Intermediate vector in the Galactic frame
    cra, sra = cos(ra), sin(ra)
    cd, sd = cos(dec), sin(dec)
    vector_gal = np.dot(Ag, np.array([cra * cd, sra * cd, sd]))

    # Galactic longitude (l) and latitude (b)
    l_gal = umath.atan2(vector_gal[1], vector_gal[0])
    b_gal = umath.atan2(vector_gal[2], sqrt(vector_gal[0]**2 + vector_gal[1]**2))

    # Calculate Galactic proper motions using the transformation matrix
    mul = -sin(l_gal) * ((-Ag[0, 0] * sra + Ag[0, 1] * cra) * pm_ra - 
                          (Ag[0, 0] * cra * cd + Ag[0, 1] * sra * sd - Ag[0, 2] * cd) * pm_dec) + \
          cos(l_gal) * ((-Ag[1, 0] * sra + Ag[1, 1] * cra) * pm_ra - 
                        (Ag[1, 0] * cra * cd + Ag[1, 1] * sra * sd - Ag[1, 2] * cd) * pm_dec)

    mub = -cos(l_gal) * sin(b_gal) * ((-Ag[0, 0] * sra + Ag[0, 1] * cra) * pm_ra - 
                                      (Ag[0, 0] * cra * cd + Ag[0, 1] * sra * sd - Ag[0, 2] * cd) * pm_dec) - \
          sin(l_gal) * sin(b_gal) * ((-Ag[1, 0] * sra + Ag[1, 1] * cra) * pm_ra - 
                                    (Ag[1, 0] * cra * cd + Ag[1, 1] * sra * sd - Ag[1, 2] * cd) * pm_dec) + \
          cos(b_gal) * ((-Ag[2, 0] * sra + Ag[2, 1] * cra) * pm_ra - 
                        (Ag[2, 0] * cra * cd + Ag[2, 1] * sra * sd - Ag[2, 2] * cd) * pm_dec)

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
    Agt = Ag.T

    # Intermediate variables for trigonometric functions
    cos_l = cos(l_gal)
    sin_l = sin(l_gal)
    cos_b = cos(b_gal)
    sin_b = sin(b_gal)

    # Calculate mu_ra_corr
    mu_ra_corr = -sra * ( (-Agt[0, 0] * sin_l + Agt[0, 1] * cos_l) * mu_l_corr\
        - (Agt[0, 0] * cos_l * cos_b + Agt[0, 1] * sin_l * sin_b - Agt[0, 2] * cos_b) * mu_b_corr )\
        + cra * ( (-Agt[1, 0] * sin_l + Agt[1, 1] * cos_l ) * mu_l_corr\
        - (Agt[1, 0] * cos_l * cos_b + Agt[1, 1] * sin_l * sin_b - Agt[1, 2] * cos_b) * mu_b_corr )

    # Calculate mu_dec_corr
    mu_dec_corr = -cra * sd * ( (-Agt[0, 0] * sin_l + Agt[0, 1] * cos_l) * mu_l_corr
        - (Agt[0, 0] * cos_l * cos_b + Agt[0, 1] * sin_l * sin_b - Agt[0, 2] * cos_b) * mu_b_corr )\
        - sra * sd * ( (-Agt[1, 0] * sin_l + Agt[1, 1] * cos_l) * mu_l_corr
        - (Agt[1, 0] * cos_l * cos_b + Agt[1, 1] * sin_l * sin_b - Agt[1, 2] * cos_b) * mu_b_corr )\
        + cd * ( (-Agt[2, 0] * sin_l + Agt[2, 1] * cos_l) * mu_l_corr
        - (Agt[2, 0] * cos_l * cos_b + Agt[2, 1] * sin_l * sin_b - Agt[2, 2] * cos_b) * mu_b_corr )

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