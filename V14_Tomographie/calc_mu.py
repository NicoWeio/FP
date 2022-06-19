import generate_table
import numpy as np
import pint
ureg = pint.UnitRegistry()

# https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html
# Daten bei 661.7 keV = 0.6617 MeV
# NOTE: Unsere Energie kann unter „Additional energies in MeV“ eingegeben werden.
# Wir nehmen an, dass es sich um INKOHÄHRENTE Streuung handelt.

STOFFE = [
    {
        'name': 'Al',
        'ρ':  ureg('2.70 kg/dm³'),
        'σ':  ureg('7.435E-02 cm²/g'),
        'σ_Photo':  ureg('6.565E-05 cm²/g'),
        'σ_Compton':  ureg('7.428E-02 cm²/g'),
    },
    {
        'name': 'Pb',
        'ρ': ureg('11.34 kg/dm³'),
        'σ': ureg('1.035E-01 cm²/g'),
        'σ_Photo': ureg('4.337E-02 cm²/g'),
        'σ_Compton': ureg('6.015E-02 cm²/g'),
    },
    {
        'name': 'Fe',
        'ρ': ureg('7.87 kg/dm³'),
        'σ': ureg('7.248E-02 cm²/g'),
        'σ_Photo': ureg('8.723E-04 cm²/g'),
        'σ_Compton': ureg('7.161E-02 cm²/g'),
    },
    {
        'name': 'Messing',
        # Cu 0.63
        # Zn 0.37
        'ρ': ureg('8.50 kg/dm³'),
        'σ': ureg('7.162E-02 cm²/g'),
        'σ_Photo': ureg('1.340E-03 cm²/g'),
        'σ_Compton': ureg('7.028E-02 cm²/g'),
    },
    {
        'name': 'Delrin',
        # POM-H
        # → CH₂O
        # abw. Quelle für die Dichte: https://www.kundert.ch/kunststoffdb4.aspx?id=21&kurzbezeichnung=POM%20H%20natur
        'ρ': ureg('1.43 g/cm³'),
        'σ': ureg('8.221E-02 cm²/g'),
        'σ_Photo': ureg('6.784E-06 cm²/g'),
        'σ_Compton': ureg('8.221E-02 cm²/g'),
    }
]


for stoff in STOFFE:
    print(f"{stoff['name']}:")
    μ = stoff['ρ'] * stoff['σ']
    print(f"{μ.to('1/cm'):.4f}")
    stoff['μ'] = μ


# tabledata = [
#     s[attr] for s in STOFFE for attr in ['name', 'ρ', 'σ', 'σ_Photo', 'σ_Compton', 'μ']
# ]

name_list = [s['name'] for s in STOFFE]
ρ_list = [s['ρ'] for s in STOFFE]
σ_list = [s['σ'] for s in STOFFE]
σ_Photo_list = [s['σ_Photo'] for s in STOFFE]
σ_Compton_list = [s['σ_Compton'] for s in STOFFE]
μ_list = [s['μ'] for s in STOFFE]

generate_table.generate_table_pint(
    'build/tab/mu.tex',
    ('Stoff', None, name_list),
    ('ρ', ureg.g / ureg.cm**3, ρ_list, 2),
    ('σ', ureg.mm**2 / ureg.g, σ_list, 3),
    (r'σ_\text{Photo}', ureg.mm**2 / ureg.g, σ_Photo_list, 3),
    (r'σ_\text{Compton}', ureg.mm**2 / ureg.g, σ_Compton_list, 3),
    ('μ', ureg.cm**(-1), μ_list, 3),
)
