import activity
import background
import efficiency
import energy_calibration
import gamma_spectrum
import nuclide_identification
from common import console

console.rule("0. Untergrund")
background.main()


console.rule("1a. Energiekalibration: Eu-152")
energy_calibration_results = energy_calibration.main()
channel_to_E = energy_calibration_results['channel_to_E']


console.rule("1b. Vollenergienachweiswahrscheinlichkeit / Effizienz: Eu-152")
efficiency_results = efficiency.main(**energy_calibration_results)

console.rule("2. Spektrum von Cs-137")
gamma_spectrum.main(channel_to_E)


console.rule("3. Aktivitätsbestimmung: Ba-133")
activity.main(channel_to_E=channel_to_E, **efficiency_results)


console.rule("4. Nuklididentifikation und Aktivitätsbestimmung: Uran & Zerfallsprodukte")
nuclide_identification.main(channel_to_E, **efficiency_results)
