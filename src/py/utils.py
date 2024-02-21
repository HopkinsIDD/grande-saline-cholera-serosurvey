import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import datetime
import logging

codedir = os.path.dirname(os.path.abspath(__file__))
rootdir = os.path.abspath(codedir+"/../")
serosurvey_dfdir = os.path.abspath(rootdir+"/from email FF/data")


casefile = serosurvey_dfdir + "/incidenceCDC/Surveillance Data Cholera Grande Saline_ 2010 - 2013_Mars.xlsx"
serofile = serosurvey_dfdir + "/serologyCDC/haiti_sero_chol.csv"


def get_git_revision_short_hash() -> str:
    import subprocess
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()

def create_folders(path):
    from pathlib import Path
    Path(path).mkdir(parents=True, exist_ok=True)

def combine_dims(a, start=0, count=2):
    """ Reshapes numpy array a by combining count dimensions, 
        starting at dimension index start """
    s = a.shape
    return np.reshape(a, s[:start] + (-1,) + s[start+count:])


def load_cases_and_serosurvey(casefile = casefile, serofile = serofile):
    # ... cases
    cn = ["date", "semaine_epi", "cas_vus_5-", "cas_vus_5+", "cas_vus", "cas_hosp_5-", "cas_hosp_5+", "cas_hosp"]
    cases_df = pd.read_excel(
        io = casefile,
        usecols = "A:H",
        parse_dates=["date"],
        skiprows = 10,
        nrows = 1096,
        names = cn
    )
    cases_df = cases_df.drop("semaine_epi", axis = 1)

    # ... serosurvey
    serosurvey_df = pd.read_csv(
        serofile,
        parse_dates = ["DateDiar", "Dateadmit", "date"]) #ints and bools come out as floats because there is NAs
    serosurvey_df["titer"] = np.round(serosurvey_df["titer"], decimals = 3)
    serosurvey_df["titer"] = serosurvey_df["titer"].replace(1.414,1.)
    serosurvey_df["titerinab"] = serosurvey_df["titerinab"].replace(1.414,1.)
    # error in serosurvey_df file:
    serosurvey_df["titer"] = serosurvey_df["titer"].replace(40560,40960) 
    serosurvey_df["titerinab"] = serosurvey_df["titerinab"].replace(40560,40960)

    serosurvey_df["titer_log2"] = np.log2(serosurvey_df["titer"])
    serosurvey_df["titerinab_log2"] = np.log2(serosurvey_df["titerinab"])

    # take the max, fill with NA when one is NA:
    serosurvey_df["titermax_log2"] = serosurvey_df[["titerinab_log2", "titer_log2"]].max(axis=1, skipna=False)
    # we only input the NA's from titerlog2, as measurement with only Inaba titers are not the safest. If you
    # believe otherwise you may uncomment the line to do this, that adds four individuals
    serosurvey_df["titermax_log2"] = serosurvey_df["titermax_log2"].fillna(serosurvey_df["titer_log2"]).fillna(serosurvey_df["titerinab_log2"])
    # date_idx makes the correspondance between serosurvey_df and cases_all_df dates
    # (ie same dates gives same date_model)
    ti = cases_df["date"][0].to_pydatetime().date()           # first cases day
    tf = serosurvey_df["date"].max().to_pydatetime().date() # last serosurvey day

    cases_df = cases_df.loc[cases_df["date"].dt.date <= tf]

    cases_df["date_idx"] = (cases_df["date"].dt.date - ti).apply(lambda x : x.days)
    serosurvey_df["date_idx"] = (serosurvey_df["date"].dt.date - ti).apply(lambda x : x.days)
    
    cases_df.set_index("date", drop = True, inplace = True)

    print(f">>> ti: {ti}, tf: {tf}")
    print(f">>> cases_df: {cases_df.shape}")
    print(f">>> serosurvey_df: {serosurvey_df.shape}")

    assert (cases_df["cas_vus_5-"] + cases_df["cas_vus_5+"] == cases_df["cas_vus"]).all()

    return cases_df, serosurvey_df, ti, tf
