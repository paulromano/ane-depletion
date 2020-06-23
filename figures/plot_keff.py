import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
from openmc.deplete import ResultsList
import serpentTools
from uncertainties import unumpy as unp

# customizations
rcParams['axes.labelsize'] = 'large'
rcParams['axes.axisbelow'] = True
rcParams['font.family'] = 'serif'
rcParams['pdf.use14corefonts'] = True
rcParams['savefig.bbox'] = 'tight'
rcParams['font.size'] = 12.0
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = [
    r'\usepackage{amsmath}',
]

# Read Serpent results --------------------------------------------------------

def get_serpent_keff(filename):
    res = serpentTools.read(filename)
    return unp.uarray(*res.resdata['absKeff'].T)

serpent_pwr = get_serpent_keff('pwr/serpent_input_res.m')
serpent_sfr = get_serpent_keff('sfr/serpent_met1000_2d_res.m')
serpent_sfr_nop = get_serpent_keff('sfr/no_ptables/serpent_met1000_2d_res.m')

# Read OpenMC results ---------------------------------------------------------

results_full_thermal = ResultsList.from_hdf5("pwr/depletion_results_full_thermal.h5")
results_full_average = ResultsList.from_hdf5("pwr/depletion_results_full_average.h5")
results_reduced = ResultsList.from_hdf5("pwr/depletion_results_reduced_average.h5")
results_casl_thermal = ResultsList.from_hdf5("pwr/depletion_results_casl_thermal.h5")
results_casl_cutoff = ResultsList.from_hdf5("pwr/depletion_results_casl_cutoff.h5")
results_casl_average = ResultsList.from_hdf5("pwr/depletion_results_casl_average.h5")

results_sfr_full = ResultsList.from_hdf5('sfr/depletion_results_full_average.h5')
results_sfr_casl = ResultsList.from_hdf5('sfr/depletion_results_casl_average.h5')
results_sfr_full_nop = ResultsList.from_hdf5('sfr/no_ptables/depletion_results_full_average.h5')
results_sfr_casl_nop = ResultsList.from_hdf5('sfr/no_ptables/depletion_results_casl_average.h5')

def difference(results, serpent_keff):
    # Obtain K_eff as a function of time
    time, k = results.get_eigenvalue()
    openmc_keff = unp.uarray(*k.T)
    diff = 1e5*(openmc_keff - serpent_keff)
    return unp.nominal_values(diff)

# Plot PWR results -------------------------------------------------------------

burnup = [
    0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
    12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5,
    40.0, 42.5, 45.0, 47.5, 50.0
]

fig, ax = plt.subplots()
ax.plot(burnup, difference(results_full_average, serpent_pwr), 'bo', label='Full, Interpolated FPY')
ax.plot(burnup, difference(results_full_thermal, serpent_pwr), 'bx', label='Full, Thermal FPY')
ax.plot(burnup, difference(results_reduced, serpent_pwr), 'ro', label='Reduced, Interpolated FPY')
ax.plot(burnup, difference(results_casl_average, serpent_pwr), 'go', label='CASL, Interpolated FPY')
ax.plot(burnup, difference(results_casl_thermal, serpent_pwr), 'gx', label='CASL, Thermal FPY')
ax.plot(color='k', linestyle='--')
ax.set_xlabel("Burnup [MWd/kg]")
ax.set_ylabel(r"$k_\text{OpenMC} - k_\text{Serpent}$ [pcm]")
ax.grid(True)
ax.legend()
plt.savefig('pwr_keff.pdf', bbox_inches='tight')
plt.close()

# Plot SFR results -------------------------------------------------------------

days = np.linspace(0.0, 360.0, 19)

fig, ax = plt.subplots()
ax.plot(days, difference(results_sfr_full, serpent_sfr), 'bo', label='Full, ptables')
ax.plot(days, difference(results_sfr_casl, serpent_sfr), 'go', label='CASL, ptables')
ax.plot(days, difference(results_sfr_full_nop, serpent_sfr_nop), 'bs', markerfacecolor='none', label='Full, no ptables')
ax.plot(days, difference(results_sfr_casl_nop, serpent_sfr_nop), 'gs', markerfacecolor='none', label='CASL, no ptables')
ax.plot(color='k', linestyle='--')
ax.set_xlabel("Time [days]")
ax.set_ylabel(r"$k_\text{OpenMC} - k_\text{Serpent}$ [pcm]")
ax.grid(True)
ax.legend()
plt.savefig('sfr_keff.pdf', bbox_inches='tight')
plt.close()
