from . import getsyst
import pandas as pd
import os

PYANA_PATH = os.getenv('PYANA')
PYANA_PATH = '/exp/sbnd/app/users/brindenc/analyze_sbnd/pyana'
systs = pd.read_csv(f'{PYANA_PATH}/syst/systematics.csv')
names = [n for n in systs["name"].values if "Geant4" in n]

def g4syst(f, nuind):
    if len(names) == 0:
        raise ValueError("No Geant4 systematics found")
    else:
        print('*'*60)
        print('Storing G4 systematics for: ')
        [print(n) for n in names]
        print('*'*60)
    return getsyst.getsyst(f, names, nuind)