import requests

data = {
    "character": 'tab',
    "Name": '',
    "Method": '2',
    "Formula": 'CH2O',
    "NumAdd": '1',
    "Energies": '0.6617',
    "Output": '',
    "WindowXmin": '0.001',
    "WindowXmax": '100000',
    "photoelectric": 'on',
    "incoherent": 'on',
    "without": 'on',
}

response = requests.post('https://physics.nist.gov/cgi-bin/Xcom/data.pl', data=data)

print(response.text)
