# Intrinsic Proper Motion Calculator

This Python package calculates the intrinsic proper motion of a star within the Galaxy, relative to its surrounding medium. The calculation subtracts the Galactic rotation component from the observed proper motion.

## Usage

1. Create a file named `<source_name>.txt` containing the required star data.
2. Run the script:

   ```bash
   python transform.py <source_name> [<constants>]
   ```

Note that there are different options for the Oort constants which can affect the results.

## Dependencies

This package requires the following Python libraries:

    numpy
    uncertainties
    astropy

Ensure these dependencies are installed before running the script.


This tool was developed for investigating stellar bow shocks. If you find this tool beneficial for your research, please consider citing our work:

## Citation:

Martínez, J.~R., del Palacio, S., & Bosch-Ramon, V. (2023). Probing the non-thermal physics of stellar bow shocks using radio observations. *Astronomy & Astrophysics*, 680, A99. [DOI: 10.1051/0004-6361/202347720](https://doi.org/10.1051/0004-6361/202347720).

**BibTeX Entry:**

```bibtex
@ARTICLE{2023A&A...680A..99M,
       author = {Martínez, J.~R. and del Palacio, S. and Bosch-Ramon, V.},
        title = "{Probing the non-thermal physics of stellar bow shocks using radio observations}",
      journal = {Astronomy & Astrophysics},
     keywords = {radiation mechanisms: non-thermal, radiation mechanisms: thermal, acceleration of particles, shock waves, radio continuum: general, Astrophysics - High Energy Astrophysical Phenomena, Astrophysics - Solar and Stellar Astrophysics},
         year = 2023,
        month = dec,
       volume = {680},
          eid = {A99},
        pages = {A99},
          doi = {10.1051/0004-6361/202347720},
archivePrefix = {arXiv},
       eprint = {2310.18669},
 primaryClass = {astro-ph.HE},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023A&A...680A..99M},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```