import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle, Patch
import openmc.deplete
import openmc.data

# customizations
mpl.rcParams['axes.autolimit_mode'] = 'round_numbers'
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['axes.xmargin'] = 0
mpl.rcParams['axes.ymargin'] = 0
mpl.rcParams['axes.axisbelow'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['font.size'] = 12.0

chain_full = openmc.deplete.Chain.from_xml('/home/romano/openmc-data/depletion/chain_endfb71_pwr.xml')
chain_reduced = openmc.deplete.Chain.from_xml('/home/romano/openmc-data/depletion/chain_reduced2_pwr.xml')
chain_casl = openmc.deplete.Chain.from_xml('/home/romano/openmc-data/depletion/chain_casl_pwr.xml')

def get_patch(z, a, **kwargs):
    n = a - z
    return Rectangle((n-0.5, z-0.5), 1, 1, **kwargs)

fig, ax = plt.subplots()

style_full = {'alpha': 0.3}
style_reduced = {}
style_casl = {'color': 'orange'}

for nuc in chain_full.nuclides:
    z, a, m = openmc.data.zam(nuc.name)
    if m > 0: continue
    if nuc.name not in chain_reduced:
        ax.add_patch(get_patch(z, a, **style_full))

for nuc in chain_reduced.nuclides:
    z, a, m = openmc.data.zam(nuc.name)
    if m > 0: continue
    if nuc.name not in chain_casl:
        ax.add_patch(get_patch(z, a, **style_reduced))

for nuc in chain_casl.nuclides:
    z, a, m = openmc.data.zam(nuc.name)
    if m > 0: continue
    ax.add_patch(get_patch(z, a, **style_casl))

ax.set_xlim(0, 175)
ax.set_ylim(0, 120)
ax.set_xlabel('N')
ax.set_ylabel('Z')
ax.grid()
legend_elements = [
    Patch(label='Full', **style_full),
    Patch(label='Reduced', **style_reduced),
    Patch(label='CASL', **style_casl)
]
ax.legend(handles=legend_elements, framealpha=1.0)

plt.savefig('chains.pdf')
