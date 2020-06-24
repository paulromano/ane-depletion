import openmc.data
from tabulate import tabulate

nfy = openmc.data.FissionProductYields('/opt/data/endf/endf-b-vii.1/nfy/nfy-092_U_235.endf')
energies = nfy.energies
tc109 = [yi['Tc109'].n for yi in nfy.independent]
mo109 = [yi['Mo109'].n for yi in nfy.independent]
sn129 = [yi['Sn129'].n for yi in nfy.independent]

print(tabulate(list(zip(energies, tc109, mo109, sn129)),
               headers=['Energy [eV]', 'Tc109', 'Mo109', 'Sn129']))
