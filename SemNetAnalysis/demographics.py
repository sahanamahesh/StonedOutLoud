

import glob
import pandas as pd
import numpy as np

WDIR = '/Users/labmanager/Library/CloudStorage/OneDrive-UBC/StonedOutLoud/Derivatives/all_subs/'
demo = glob.glob(WDIR+"/Demographics/*.csv")
semnet = pd.read_excel(WDIR+'SemNetAnalysis/network_stats.xlsx')
for group in demo:
    new = pd.DataFrame()
    grp = group.rsplit('/',1)[1].split('.')[0]
    
    data = pd.read_csv(group)
    new = data[data['Subject'].isin(np.unique(np.array(semnet['Subject'])))]
    if len(str(data['Subject'][0])) !=3:
        for i, sub in data['Subject'].items():
            og_sub = sub
            if len(str(sub)) ==1:
                sub = '00'+str(sub)
            if len(str(sub)) ==2:
                sub = '0'+str(sub)
            if str(sub) in list(semnet['Subject']):
                spot = len(new)
                new.loc[spot] = data.loc[i]
                new['Subject'] = new['Subject'].replace(og_sub, sub)
    new.to_csv(WDIR+f"SemNetAnalysis/Demographics/{grp}.csv", index=False)
    
    print(new)

