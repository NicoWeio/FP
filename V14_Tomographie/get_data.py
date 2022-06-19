import requests

common_data = {
    "character": 'tab',
    "Name": '',
    "NumAdd": '1',
    "Energies": '0.6617',
    "Output": '',
    "WindowXmin": '0.001',
    "WindowXmax": '100000',
    "photoelectric": 'on',
    "incoherent": 'on',
    "without": 'on',
}

# data = {
#     **common_data,
#     "Method": '2',
#     "Formula": 'CH2O',
# }

data = {
    **common_data,
    "Method": '1',
    "OutOpt": 'PIC',
    'ZNum': 13, # Al (eigentlich 13!?)
}

response = requests.post('https://physics.nist.gov/cgi-bin/Xcom/data.pl', data=data)
lines = [l.split('\t') for l in response.text.strip().split('\n')]

try:
    assert len(lines) == 4
except AssertionError:
    print(response.text)
    raise

names = [l1.strip() + ' ' + l2.strip() for l1, l2 in zip(lines[0], lines[1]) if l1.strip() != '']
values = lines[3]

data = dict(zip(names, values))

print(data)
