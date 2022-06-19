import numpy as np
import pint
ureg = pint.UnitRegistry()

# https://physics.nist.gov/cgi-bin/Xcom/xcom3_1
# Daten bei 6.000E-01 ↓

STOFFE = [
    {
        'name': 'Al',
        'ρ':  ureg('2.7 kg/dm^3'),
        'σ':  ureg('7.762E-2 cm^2/g'),
    },
    {
        'name': 'Pb',
        'ρ': ureg('11.34 kg/dm³'),
        'σ': ureg('1.167E-01 cm²/g'),
    },
    {
        'name': 'Fe',
        'ρ': ureg('7.87 kg/dm³'),
        'σ': ureg('7.583E-02 cm²/g'),
    },
    {
        'name': 'Messing',
        # Cu 0.63
        # Zn 0.37
        'ρ': ureg('8.5 kg/dm³'),
        'σ': ureg('7.503E-02 cm²/g'),
    },
    {
        'name': 'Delrin',
        # POM-H
        # abw. Quelle für die Dichte: https://www.kundert.ch/kunststoffdb4.aspx?id=21&kurzbezeichnung=POM%20H%20natur
        'ρ': ureg('1.43 g/cm³'),
        'σ': ureg('8.220E-02 cm²/g'),
    }
]


for stoff in STOFFE:
    print(f"{stoff['name']}:")
    μ = stoff['ρ'] * stoff['σ']
    print(f"{μ.to('1/cm'):.4f}")
