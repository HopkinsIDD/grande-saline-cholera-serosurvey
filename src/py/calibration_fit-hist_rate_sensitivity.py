import utils
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as plt
import arviz as az
import pymc as pm
print(f"Running on PyMC3 v{pm.__version__}")
print(f"Running on ArviZ v{az.__version__}")
#import pytensor as pt
#import pytensor.tensor as ptt
import aesara.tensor as ptt
import aesara as pt
import matplotlib.dates as mdates
import xarray as xr

import os, datetime, pickle, shutil, bz2, scipy, math, tqdm, sys

def combine_dims(a, start=0, count=2):
    """ Reshapes numpy array a by combining count dimensions, 
        starting at dimension index start """
    s = a.shape
    return np.reshape(a, s[:start] + (-1,) + s[start+count:])

# WARNING NEEDS TO BE INSIDE THE if __name__ == "__main__" BLOCK !! for it to work in each thread

spec_id = int(sys.argv[1])

multipliers = np.concatenate((1/np.arange(2,11)[::-1], np.arange(1,11)))  # len is 19, from 1/10th to 10 times
this_mult=multipliers[spec_id]

## First
interactive = False
save = False
test_value = False  # wether to compute aesera test value
prefix = "cc-Fit-hist"

dynamics_specs = ["mvDynamics"]  # "inDynamics",
age_spec = ["U5", "adult", "all"]
model_specs = [f"{a}-{d}" for a in age_spec for d in dynamics_specs]
model_spec = model_specs[1] # hardcoded (0: children, 1: adult)
model_id = f"{prefix}_{model_spec}_{spec_id}_{utils.get_git_revision_short_hash()}_{datetime.date.today()}"
model_folder = f"/work/users/c/h/chadi/calib_Fit-hist/{model_id}"
print(f">>> {spec_id} -> using model spec {model_spec} \n          >>> {model_specs}")
print(f">>> saving to {model_id}")

# can be modified 
titer_var = "titermax_log2" # var to model

# Load data
cases_df, serosurvey_df, ti, tf = utils.load_cases_and_serosurvey()
age_spec = model_spec.split("-")[0]

population = 21131 # population from IHSI, 2009
# https://en.wikipedia.org/wiki/Demographics_of_Haiti for % under5 
# and age 2,3,4 should be roughly 3/5 of age 0,1,2,3,4
population2_4 = int(population * 12.53/100 * 3/5) 
population5_99 = int(population * (100-12.53)/100)
population2_99 = population2_4 + population5_99
probcase2_4 = 0.7394 # probability of a case being 2_4 the cases with MSF data Gonaives, Nov and Dec 2010

def get_gt_data(age_spec):
    if age_spec == "all":
        n_population = population2_99
        cases_gt = cases_df["cas_vus_5-"] * probcase2_4 + cases_df["cas_vus_5+"]
        serosurvey_gt = serosurvey_df
    elif age_spec == "U5":
        n_population = population2_4
        cases_gt = cases_df["cas_vus_5-"] * probcase2_4
        serosurvey_gt = serosurvey_df[serosurvey_df["Age"] < 5]
    elif age_spec == "adult":
        n_population = population5_99
        cases_gt = cases_df["cas_vus_5+"]
        serosurvey_gt = serosurvey_df[serosurvey_df["Age"] >= 5]

    cases_arr = np.array(cases_gt)
    n_sampled = len(serosurvey_gt)
    cases_total = np.sum(cases_arr) 
    n_days = len(cases_arr)

    # get all bins, even if zero for this age
    observed_count_index = serosurvey_df.groupby(by=titer_var)[titer_var].agg("count")
    observed_count = serosurvey_gt.groupby(by=titer_var)[titer_var].agg("count")
    observed_count = observed_count.reindex_like(observed_count_index).fillna(0)
    observed_count = observed_count.astype(int) # important for pymc bug
    
    return cases_gt, serosurvey_gt, n_population, cases_arr, n_sampled, cases_total, n_days, observed_count


cases_gt, serosurvey_gt, n_population,\
    cases_arr, n_sampled, cases_total,\
    n_days, observed_count = get_gt_data(age_spec)

print(f">>> selected age_spec {age_spec} with population: {n_population}, n_case: {cases_total}, n_sampled: {n_sampled}")


# data for model:
bin_lower_lim = observed_count.index.values

# some constants for plots
bins_pretty = [f"1:{int(np.round(2**x))}" for x in bin_lower_lim]
bins_pretty[0] = "0"
bins_pretty[7:] = [f"{b[:-3]}'{b[-3:]}" for b in bins_pretty[7:]]
bins_pretty_wlog = [f"{t} ({lt:4.2f})" for t, lt in zip(bins_pretty, bin_lower_lim)]
to_plot_vars = ["rho", "sigma", "delta", "omega_p", "lambda"]
sns.set_style('ticks')
#sns.set_context("paper", font_scale=2)
plt.rcParams['figure.facecolor'] = 'white'
sns.set_palette("dark")


# solve the nested braket bug
#pt.config.gcc__cxxflags = pt.config.gcc__cxxflags + " -fbracket-depth=800"

if __name__=="__main__":
    ## START COPY PASTED CODE
    if test_value:
        pt.config.compute_test_value = 'warn'
    else:
        pt.config.compute_test_value = 'off'

    sampling_time_mean = np.round(serosurvey_gt["date_idx"].mean()) # mean sampling time.
    sampling_time = np.array(np.round(serosurvey_gt["date_idx"]))
    # cases_df[-10:] is zeros, so no need to model the increase before exponential titer decay.
    delay_report2peaktiter = 9 # days

    # lower bins
    lbins_tiled = np.tile(bin_lower_lim, [n_sampled,1]).T 

    # remove the deterministic variable before run.
    with pm.Model() as titer_indivdecay_model4:
        # reporting rate:
        # reporting rate is allowed of overreporting (above 1) and underreporting (below 1)
        rho =  pm.TruncatedNormal("rho", 
                mu = (cases_total/n_population+1.5)/2, 
                sigma = 1, 
                lower = cases_total/n_population, 
                upper = 1.5, 
                shape = 1)
        
        # decay parameters:
        dynamics_spec = model_spec.split("-")[1]
        if dynamics_spec == "inDynamics":
            print("using independent dynamics")
            omega_p = pm.TruncatedNormal("omega_p", 
                    mu = 2+np.log2(5),#2+np.log2(5), 
                    sigma = 1.75, 
                    shape = n_sampled, 
                    lower = 0,
                    #upper=5+np.log2(5)
                    )
            llambda = pm.TruncatedNormal("lambda", 
                    mu = np.log2(45), # from Jones et al, table S3
                    sigma = np.log2(40)/2,#np.log2((40+150)/100)/1.96, , 
                    shape = n_sampled, 
                    lower = np.log2(20)#,upper = np.log2(250)
                    )
        elif dynamics_spec == "mvDynamics":
            print("Using multivariate dynamics")
            titer_params = pm.MvNormal("titer_params", 
                    mu = [2.98825782+np.log2(5), np.log2(45)], 
                    cov= np.array([[3.9655, -4.0231],
                                    [-4.0231, 8.1027]]), 
                    shape=(n_sampled, 2))
            omega_p_t = pm.Deterministic("omega_p_t", titer_params[:,0])
            llambda_t = pm.Deterministic("lambda_t", titer_params[:,1])
            # force at the mean in case it is < 0
            omega_p =  pm.Deterministic("omega_p", ptt.switch(omega_p_t>0, omega_p_t, 2.98825782+np.log2(5)))
            llambda =  pm.Deterministic("lambda", ptt.switch(llambda_t>0, llambda_t, np.log2(45)))
        else:
            raise ValueError(f"dynamic_spec {dynamics_spec} not recognized for dynamics")
        
        delta = pm.TruncatedNormal("delta", 
                mu = this_mult*np.log(2)/123,  # half life of 123
                sigma = .0005,
                shape = 1,
                lower = np.log(2)/250,
                upper = np.log(2)/20)

        # titer sd around mean
        sigma = pm.HalfNormal("sigma", 
                sigma = 1, 
                shape = 1)


        #sigma_p = pt.printing.Print('sigma')(sigma)
        #rho_p = pt.printing.Print('rho')(rho)
        #omega_p_p = pt.printing.Print('omega_p')(omega_p)
        #llambda_p = pt.printing.Print('lambda')(llambda)
        #delta_p = pt.printing.Print('delta')(delta)

        day_inf_reported = pm.Categorical("day_inf_reported", 
                p = cases_arr/cases_total, 
                shape = n_sampled)

        # TODO: this variable should be named if_infected only
        if_inf_reported = pm.Bernoulli("if_inf_reported",  
                p = (cases_total/rho)/n_population,
                shape = n_sampled)

        #if_inf_reported_p = pt.printing.Print('if_inf_reported')(if_inf_reported)
        #day_inf_reported_p = pt.printing.Print('day_inf_reported')(day_inf_reported)  # Warning: if I do sum or mean, the value is sent !!

        mean_titer_all = omega_p + if_inf_reported * llambda * ptt.exp(-delta * (sampling_time - (day_inf_reported + delay_report2peaktiter)))

        # 3. Measurement model: (exponential of the log) cumulative densitify function 
        # of the titer distribution around the mean titer evaluated on the measurement bins.
        titer_bins_per_event1 = ptt.exp(
            pm.logcdf(rv=pm.Normal.dist(
                            mu = mean_titer_all,
                            sigma = sigma, 
                            shape = n_sampled
                            ), 
                    value=lbins_tiled[1:,:]))  # here important to exclude the first zero
        # bounds of the bined probabilities, and then diff to get the probability to 
        # be in each bin for each event.
        titer_bins_per_event2 = ptt.concatenate([np.zeros((1,n_sampled)), titer_bins_per_event1, np.ones((1,n_sampled))], axis=0)
        titer_bins_per_event = pm.Deterministic("titer_bins_per_event", ptt.extra_ops.diff(titer_bins_per_event2, axis=0)) # 13 x 2622

        # bin probabilities weighted by the event probabilities
        #titer_bins_probs_weighted = pm.Deterministic("titer_bins_probs_weighted", event_prob_all * titer_bins_per_event)  # 13 x 2622
        
        titer_bin_probs = titer_bins_per_event.sum(axis=1)/n_sampled
        # titer_bin_probs_p = pt.printing.Print('titer_bin_probs')(titer_bin_probs)
        #titer_bin_probs_p = pt.printing.Print('titer_bin_probs')(titer_bin_probs)

        # observed counts as a Multinomial of the sum of all events prob_per_bin
        modeled_count = pm.Multinomial("modeled_count", 
                    p = titer_bin_probs, 
                    n = observed_count.sum(), 
                    observed = observed_count.to_numpy())

    model_catcount = titer_indivdecay_model4
    ## END COPY PASTED CODE

    if True:
        #if dynamics_spec == "mvDynamics":
        #    initvals = {}
        if True:
            initvals = {
                    "lambda":np.random.uniform(4.5, 16, n_sampled), 
                    "if_inf_reported": np.random.binomial(n=1, p=.2, size=n_sampled),
                    "day_inf_reported": pm.Categorical.dist(p=cases_arr/cases_total, size=n_sampled).eval()
                    }
        chainlength = 2000
        tunelength = 10000
        nchains = 4
        with model_catcount:
            trace =  pm.sample( #pm.sampling_jax.sample_numpyro_nuts(
            chainlength,
            tune = tunelength,
            chains = nchains,
    #        init = "adapt_full",
            initvals=initvals,
            #nuts_kwargs = dict(target_accept = 0.99)
            )
        
        utils.create_folders(model_folder)

        trace.to_netcdf(f'{model_folder}/fit_{model_id}.nc')

        trace_raw = az.from_netcdf(f'{model_folder}/fit_{model_id}.nc')
        trace_thinned = trace_raw

        with model_catcount:
            posterior_samples = pm.sample_posterior_predictive(trace_thinned,
                progressbar=True,
                var_names=['if_inf_reported', 'day_inf_reported', "modeled_count", "titer_bins_per_event", ]
                )

        trace_raw.extend(posterior_samples)

        n_sample=10000
        with model_catcount:
                prior_samples = pm.sample_prior_predictive(samples=n_sample)

        trace_raw.extend(prior_samples)

        trace_raw.to_netcdf(f'{model_folder}/full_{model_id}.nc')

