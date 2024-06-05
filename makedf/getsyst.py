import uproot
import numpy as np
import pandas as pd

def getsyst(f, systematics, nuind):
    if "globalTree" not in f:
        return pd.DataFrame(index=nuind.index)

    nuidx = pd.MultiIndex.from_arrays([nuind.index.get_level_values(0), nuind])

    wgt_names = [n for n in f["globalTree"]['global/wgts/wgts.name'].arrays(library="np")['wgts.name'][0]]
    wgt_types = f["globalTree"]['global/wgts/wgts.type'].arrays(library="np")['wgts.type'][0]
    wgt_nuniv = f["globalTree"]['global/wgts/wgts.nuniv'].arrays(library="pd")['wgts.nuniv'][0]
    #Write systematics df
    syst_dict = {i: [None, None, None] for i in range(len(wgt_names))}
    for i, n in enumerate(wgt_names):
        #print(f'{i}: {n} (type {wgt_types[i]}, {wgt_nuniv[i]} universes)')
        syst_dict[i][0] = n
        syst_dict[i][1] = wgt_types[i]
        syst_dict[i][2] = wgt_nuniv[i]
    syst_df = pd.DataFrame.from_dict(syst_dict, orient='index', columns=['name', 'type', 'nuniv'])
    #syst_df.to_csv('systematics.csv', index=False)

    isyst = pd.Series(np.repeat(list(range(len(wgt_nuniv))), wgt_nuniv), name="isyst")
    isyst.index.name = "iwgt"
    nuniv = wgt_nuniv.sum()

    wgts = f["recTree"]['rec.mc.nu.wgt.univ'].arrays(library="pd")
    wgts["inu"] = wgts.index.get_level_values(1) // nuniv
    wgts["iwgt"] = wgts.index.get_level_values(1) % nuniv
    wgts = wgts.reset_index().set_index(["entry", "inu", "iwgt"]).drop(columns="subentry")
    wgts.columns = ["wgt"]
    wgts = wgts.join(isyst)

    systs = []
    for s in systematics:
        isyst = wgt_names.index(s)
        this_systs = []

        # Get weight type
        # +/- 1,2,3 sigma
        if wgt_types[isyst] == 3 and wgt_nuniv[isyst] == 1: # morph unisim
            s_morph = wgts[wgts.isyst == isyst].wgt.groupby(level=[0,1]).first()
            s_morph.name = (s, "morph")

            this_systs.append(s_morph)
        elif wgt_types[isyst] == 3 and wgt_nuniv[isyst] > 1: # +/- sigma unisim
            nsigma = wgts[wgts.isyst == isyst].wgt.groupby(level=[0,1]).size().values[0] // 2
            for isigma in range(nsigma):
                s_ps = wgts[wgts.isyst == isyst].wgt.groupby(level=[0,1]).nth(2*isigma)
                s_ps.name = (s, "ps%i" % (isigma+1))
                s_ms = wgts[wgts.isyst == isyst].wgt.groupby(level=[0,1]).nth(2*isigma+1)
                s_ms.name = (s, "ms%i" % (isigma+1))
 
                this_systs.append(s_ps)
                this_systs.append(s_ms)

        elif wgt_types[isyst] == 0: # multisim
            this_wgts = wgts[wgts.isyst == isyst].wgt.groupby(level=[0,1]).head(250) # limit to 250 universes
            this_wgts = this_wgts.reset_index(level=2)
            this_wgts = this_wgts.pivot_table(values="wgt", index=["entry", "inu"], columns="iwgt")
            this_wgts.columns = pd.MultiIndex.from_tuples([(s, "univ_%i"% i) for i in range(len(this_wgts.columns))])

            for c in this_wgts.columns:
                this_systs.append(this_wgts[c])

        else:
            raise Exception("Cannot decode systematic uncertainty: %s" % s)

        for syst in this_systs:
            systs.append(syst)

    systs = pd.DataFrame(systs).T

    s_idx = systs.index.get_indexer(nuidx)
    systs_match = systs.iloc[s_idx]
    systs_match.loc[s_idx < 0, :] = 1.
    systs_match.index = nuind.index

    return systs_match