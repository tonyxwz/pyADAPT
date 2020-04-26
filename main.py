from pyADAPT.examples.toy import ToyModel
import numpy as np
import pandas as pd
from pprint import pprint

t = ToyModel()
data_names = ["s2", "s3"]

for n in data_names:
    print(n in t.states['name'])

for n in t.states['name']:
    print(n in data_names)

x = np.arange(0, 10)
mask = np.concatenate([[True] * 5, [False] * 5])

print(mask)

print(x[mask])
j = pd.DataFrame(data={'test': ['de', 1, 'eee']})
print(j.loc[0, 'test'] in ['de', 1, 'eee'])

k = ['s1', 's2', 's3', 's4']
for n in t.states.name:
    print(k.index(n))
k.insert(8, k.pop(0))

pprint(k)
print(list(t.states['name']) + t.flux_order)
print(t.parameters['k1'])