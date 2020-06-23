import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as mtick
import numpy as np
from openmc.data import ATOMIC_SYMBOL, zam
from openmc.deplete import ResultsList
import serpentTools

# customizations
rcParams['axes.autolimit_mode'] = 'round_numbers'
rcParams['axes.labelsize'] = 'large'
rcParams['axes.xmargin'] = 0
rcParams['axes.ymargin'] = 0
rcParams['axes.axisbelow'] = True
rcParams['font.family'] = 'serif'
rcParams['pdf.use14corefonts'] = True
rcParams['savefig.bbox'] = 'tight'
rcParams['font.size'] = 12.0
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = [
    r'\usepackage{amsmath}',
]

print('Loading PWR data...')
results_full_thermal = ResultsList.from_hdf5("pwr/depletion_results_full_thermal.h5")
results_full_cutoff = ResultsList.from_hdf5("pwr/depletion_results_full_cutoff.h5")
results_full_average = ResultsList.from_hdf5("pwr/depletion_results_full_average.h5")
results_reduced = ResultsList.from_hdf5("pwr/depletion_results_reduced_average.h5")
results_casl_thermal = ResultsList.from_hdf5("pwr/depletion_results_casl_thermal.h5")
results_casl_cutoff = ResultsList.from_hdf5("pwr/depletion_results_casl_cutoff.h5")
results_casl_average = ResultsList.from_hdf5("pwr/depletion_results_casl_average.h5")
fuel = serpentTools.read('pwr/serpent_input_dep.m')['fuel']

print('Loading SFR data...')
results_sfr_full = ResultsList.from_hdf5('sfr/depletion_results_full_average.h5')
results_sfr_casl = ResultsList.from_hdf5('sfr/depletion_results_casl_average.h5')
results_sfr_full_nop = ResultsList.from_hdf5('sfr/no_ptables/depletion_results_full_average.h5')
results_sfr_casl_nop = ResultsList.from_hdf5('sfr/no_ptables/depletion_results_casl_average.h5')
fuel_sfr = serpentTools.read('sfr/serpent_met1000_2d_dep.m')['fuel']
fuel_sfr_nop = serpentTools.read('sfr/no_ptables/serpent_met1000_2d_dep.m')['fuel']


day = 24*60*60

colors = rcParams['axes.prop_cycle'].by_key()['color']

actinides = [
    'U234', 'U235', 'U236', 'U238',
    'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242',
    'Am241', 'Am242', 'Am242_m1', 'Am243', 'Am244',
    'Cm242', 'Cm243', 'Cm244', 'Cm245', 'Cm246',
]

fission_products = [
    'Kr85', 'Sr90', 'Y90', 'Mo95', 'Tc99', 'Ru101', 'Ru106', 'Rh103', 'Pd105', 'Pd107', 'Ag109',
    'Te132', 'I129', 'Xe131', 'Xe135', 'Cs133', 'Cs134', 'Cs135', 'Cs137',
    'La139', 'Ce140', 'Ce142', 'Ce144', 'Pm147', 'Nd142', 'Nd143',
    'Nd144', 'Nd145', 'Nd146', 'Nd148', 'Nd150', 'Sm147', 'Sm148', 'Sm149',
#    'Sm150', 'Sm151', 'Sm152', 'Sm154', 'Eu151', 'Eu153', 'Eu154', 'Eu155',
#    'Gd154', 'Gd155', 'Gd156', 'Gd158', 'Gd160'
    'Zr91', 'Zr93', 'Zr95', 'Zr96',
    'Mo96', 'Mo97', 'Mo98', 'Mo99', 'Mo100'
]

def sup_label(name):
    z, a, m = zam(name)
    meta = "\\text{m}" if m > 0 else ""
    sym = ATOMIC_SYMBOL[z]
    return f"$^{{{a}{meta}}}${sym}"

burnup_pwr = [
    0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
    12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5,
    40.0, 42.5, 45.0, 50.0
]
burnup_sfr = np.linspace(0.0, 360.0, 19)

def isotope_bar_plot(results_openmc, results_serpent, nuclides, filename, decimals=2, **kwargs):
    print(f'{filename}...')
    conc_openmc = {}
    conc_serpent = {}

    for nuc in nuclides:
        # Obtain concentration as a function of time
        time, conc_openmc[nuc] = results_openmc.get_atoms('1', nuc, 'atom/b-cm')

        # Read Serpent results
        z, a, m = zam(nuc)
        conc_serpent[nuc] = results_serpent.getValues('days', 'adens', names=f'{ATOMIC_SYMBOL[z]}{a}{"m" if m else ""}')[0]


    fig, ax = plt.subplots(**kwargs)
    if 'pwr' in filename:
        burn_indices = (7, 14, 21, 28)
        burnup = burnup_pwr
        burnup_units = 'MWd/kg'
    else:
        burn_indices = (4, 9, 13, 18)
        burnup = burnup_sfr
        burnup_units = 'days'
    ind = np.arange(len(nuclides))
    height = 1/(len(burn_indices)+1)
    for i, index in enumerate(burn_indices):
        data_openmc = np.array([conc_openmc[nuc][index] for nuc in nuclides])
        data_serpent = np.array([conc_serpent[nuc][index] for nuc in nuclides])
        diff = np.zeros_like(data_openmc)
        nonzero = data_serpent > 0.0
        diff[nonzero] = (data_openmc[nonzero] - data_serpent[nonzero])/data_serpent[nonzero]
        ax.barh(ind + i*height, diff, height, label=f'{burnup[index]} {burnup_units}')

    locs = ind + height*(len(burn_indices)/2 - 0.5)
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(1.0, decimals))
    ax.set_xlabel('(OpenMC - Serpent)/Serpent')
    ax.set_yticks(locs)
    ax.set_yticklabels([sup_label(x) for x in nuclides])
    ax.set_ylim(locs[0] - 0.5, locs[-1] + 0.5)
    ax.invert_yaxis()
    ax.legend()
    ax.grid()
    plt.savefig(filename)

ext = 'pdf'
isotope_bar_plot(results_full_average, fuel, actinides, f'pwr/pwr_actinides_full.{ext}')
isotope_bar_plot(results_reduced, fuel, actinides, f'pwr/pwr_actinides_reduced.{ext}')
isotope_bar_plot(results_casl_average, fuel, actinides, f'pwr/pwr_actinides_casl.{ext}')

isotope_bar_plot(results_sfr_full, fuel_sfr, actinides, f'sfr/sfr_actinides_full.{ext}')
isotope_bar_plot(results_sfr_casl, fuel_sfr, actinides, f'sfr/sfr_actinides_casl.{ext}')
isotope_bar_plot(results_sfr_full_nop, fuel_sfr_nop, actinides, f'sfr/sfr_actinides_full_nop.{ext}', 3)
isotope_bar_plot(results_sfr_casl_nop, fuel_sfr_nop, actinides, f'sfr/sfr_actinides_casl_nop.{ext}', 4)

kwargs = {'figsize': (6.4, 10.0)}
isotope_bar_plot(results_full_thermal, fuel, fission_products, f'pwr/pwr_fp_full_thermal.{ext}', 0, **kwargs)
#isotope_bar_plot(results_full_cutoff, fuel, fission_products, f'pwr/pwr_fp_full_cutoff.{ext}', 1, **kwargs)
isotope_bar_plot(results_full_average, fuel, fission_products, f'pwr/pwr_fp_full_average.{ext}', 1, **kwargs)
isotope_bar_plot(results_reduced, fuel, fission_products, f'pwr/pwr_fp_reduced_average.{ext}', 1, **kwargs)
isotope_bar_plot(results_casl_thermal, fuel, fission_products, f'pwr/pwr_fp_casl_thermal.{ext}', 0, **kwargs)
#isotope_bar_plot(results_casl_cutoff, fuel, fission_products, f'pwr/pwr_fp_casl_cutoff.{ext}', 0, **kwargs)
isotope_bar_plot(results_casl_average, fuel, fission_products, f'pwr/pwr_fp_casl_average.{ext}', 0, **kwargs)

isotope_bar_plot(results_sfr_full, fuel_sfr, fission_products, f'sfr/sfr_fp_full_average.{ext}', 0, **kwargs)
isotope_bar_plot(results_sfr_casl, fuel_sfr, fission_products, f'sfr/sfr_fp_casl_average.{ext}', 1, **kwargs)
isotope_bar_plot(results_sfr_full_nop, fuel_sfr, fission_products, f'sfr/sfr_fp_full_average_nop.{ext}', 0, **kwargs)
isotope_bar_plot(results_sfr_casl_nop, fuel_sfr, fission_products, f'sfr/sfr_fp_casl_average_nop.{ext}', 1, **kwargs)
