import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import datetime
import logging

codedir = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.abspath(codedir+"/../../data")


casefile = datadir + "/incidence_clean.csv"
serofile = datadir + "/serology_clean.csv"


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


def load_cases_and_serosurvey_clean(casefile = casefile, serofile = serofile):
    # ... cases
    cases_df = pd.read_csv(casefile, parse_dates=["date"]) # cas_vus,cas_vus_lt5,cas_vus_geq5

    # ... serosurvey
    serosurvey_df = pd.read_csv(serofile, parse_dates = ["date"]) # date,age,vibriocidalMAX_titer,n

    # clean error in the serosurvey
    serosurvey_df["vibriocidalMAX_titer"] = serosurvey_df["vibriocidalMAX_titer"].replace(1.414,1.)
    serosurvey_df["vibriocidalMAX_titer"] = serosurvey_df["vibriocidalMAX_titer"].replace(40560,40960) 
    serosurvey_df["titermax_log2"] = np.log2(serosurvey_df["vibriocidalMAX_titer"])
    ti = cases_df["date"][0].to_pydatetime().date()           # first cases day
    tf = serosurvey_df["date"].max().to_pydatetime().date()   # last serosurvey day

    # take only cases before the last serosurvey day
    cases_df = cases_df.loc[cases_df["date"].dt.date <= tf]

    # Create a date idx column with a date index from 0 (ti) to n (tf)
    cases_df["date_idx"] = (cases_df["date"].dt.date - ti).apply(lambda x : x.days)
    serosurvey_df["date_idx"] = (serosurvey_df["date"].dt.date - ti).apply(lambda x : x.days)
    
    cases_df.set_index("date", drop = True, inplace = True)

    print(f">>> ti: {ti}, tf: {tf}")
    print(f">>> cases_df: {cases_df.shape}")
    print(f">>> serosurvey_df: {serosurvey_df.shape}")

    assert (cases_df["cas_vus_lt5"] + cases_df["cas_vus_geq5"] == cases_df["cas_vus"]).all()

    return cases_df, serosurvey_df, ti, tf
