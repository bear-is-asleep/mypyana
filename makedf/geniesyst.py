from . import getsyst
import pandas as pd
import os
PYANA_PATH = os.getenv('PYANA')
PYANA_PATH = '/exp/sbnd/app/users/brindenc/analyze_sbnd/pyana'
systs = pd.read_csv(f'{PYANA_PATH}/syst/systematics.csv')
names = [n for n in systs["name"].values if "GENIE" in n]
names = [n for n in names if "multisim" in n] #keep multim only for now

def geniesyst(f, nuind):
    if len(names) == 0:
        raise Exception("No GENIE systematics found")
    else:
        print('*'*60)
        print('Storing systematics for: ')
        [print(n) for n in names]
        print('*'*60)
    return getsyst.getsyst(f, names, nuind)

