from cmdstanpy import cmdstan_path, CmdStanModel, set_cmdstan_path
import numpy as np
import pandas as pd
from scipy import stats
import cmdstanpy

model = "independent"  # "independent" or "dynamical"

if model == "independent":
    model = CmdStanModel(
        stan_file='/n/gershman_lab/users/Stan/1_code/publish_code/stan/final_independent_models/lt_hierarchical_independent_with_gq_final.stan',
                         cpp_options={'STAN_THREADS': True})

elif model == "dynamical":
        model = CmdStanModel(
            stan_file='/n/gershman_lab/users/Stan/1_code/publish_code/stan/final_dynamical_models/lt_hierarchical_super_model_noWeights_noPrev_0effects_N01Priors_GQ_varExplained_PD_exo2.stan',
            cpp_options={'STAN_THREADS': True})

# make data ITC
data_itc_all_weeks = pd.read_csv('/n/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/itc_data_for_stan_90s.csv')

# make data lt
data_lt_all_weeks = pd.read_csv('/n/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/lt_data_for_stan_90s.csv')

all_subjects = data_lt_all_weeks.subjectId
# make data frames
subjects = all_subjects.unique()

sub_exclude = [1, 4, 7, 8, 9, 11, 14, 15, 19, 23, 24, 25, 28, 32, 35, 36, 37, 47, 50, 53, 58, 60, 61, 65, 66, 68, 70,
               74, 77, 83, 84, 88];

# Remove the subjects names, so that we will not include these subjects at all in Stan (otherwise we will have parameters for these subjects, but these won't be used in the calculation of the log-likelihood, and weill therefore be harmful unused parameters)
subjects = np.delete(subjects, [x - 1 for x in sub_exclude])

N = subjects.shape[0]

N = N
W = 12

Tr_max_itc_all = 27
P_itc = 2
x1_itc_all = np.zeros([N, W, Tr_max_itc_all])
x2_itc_all = np.zeros([N, W, Tr_max_itc_all])
x3_itc_all = np.zeros([N, W, Tr_max_itc_all])
y1_itc_all = np.zeros([N, W, Tr_max_itc_all])
Tr_itc_all = np.zeros([N, W])
W_itc_obs_all = np.zeros([N])
W_itc_mis_all = np.zeros([N])
idx_itc_obs_all = np.zeros([N, W])
idx_itc_mis_all = np.zeros([N, W])

Tr_max_lt_all = 30
P_lt = 2
W_lt_obs_all = np.zeros([N])
W_lt_mis_all = np.zeros([N])
idx_lt_obs_all = np.zeros([N, W])
idx_lt_mis_all = np.zeros([N, W])
Tr_lt_all = np.zeros([N, W])
hi_p_all = np.zeros([N, W, Tr_max_lt_all])
hi_narr_all = np.zeros([N, W, Tr_max_lt_all])
hi_wide_all = np.zeros([N, W, Tr_max_lt_all])
lo_narr_all = np.zeros([N, W, Tr_max_lt_all])
lo_wide_all = np.zeros([N, W, Tr_max_lt_all])
choice_lt_all = np.zeros([N, W, Tr_max_lt_all])

sub = 0
for s in subjects[0:N]:
    # for s in sub_no_miss[0:N]:

    # ITC
    ########
    data_itc_tmp = data_itc_all_weeks.loc[data_itc_all_weeks['subjectId'] == s]

    W_itc_obs = data_itc_tmp['weekId'].unique().shape[0]
    W_itc_mis = W - W_itc_obs

    idx_itc_obs = data_itc_tmp['weekId'].unique()
    idx_itc_mis = np.setdiff1d(np.array(range(1, W + 1)), idx_itc_obs)

    tr_tmp = np.zeros([W])
    c = 0
    for w in data_itc_tmp['weekId'].unique():
        tr_tmp[idx_itc_obs[c] - 1] = np.shape(data_itc_tmp.loc[data_itc_tmp['weekId'] == w])[0]
        c = c + 1

    Tr_itc = np.transpose(tr_tmp.astype(int)).reshape(W)

    a = np.zeros(shape=(W, Tr_max_itc_all))
    x1_itc = np.zeros(a.shape, dtype=object)
    x2_itc = np.zeros(a.shape, dtype=object)
    x3_itc = np.zeros(a.shape, dtype=object)
    y1_itc = np.zeros(a.shape, dtype=object)

    # data_itc_tmp['choice'] = data_itc_tmp['choice'].replace({'now': 0, 'later': 1})
    tr_tmp = np.zeros([W])
    c = 0
    for w in data_itc_tmp['weekId'].unique():
        tr_tmp = np.shape(data_itc_tmp.loc[data_itc_tmp['weekId'] == w])[0]
        x1_tmp = np.array(data_itc_tmp['later_delay'].loc[data_itc_tmp['weekId'] == w])
        x2_tmp = np.array(data_itc_tmp['large_amount'].loc[data_itc_tmp['weekId'] == w])
        x3_tmp = np.array(data_itc_tmp['small_amount'].loc[data_itc_tmp['weekId'] == w])
        tmp = data_itc_tmp['choice'].loc[data_itc_tmp['weekId'] == w]

        if tr_tmp < Tr_max_itc_all:
            padd = int(Tr_max_itc_all - tr_tmp)
            z = np.zeros([padd, 1])
            x1_tmp = np.append(x1_tmp, z)
            x2_tmp = np.append(x2_tmp, z)
            x3_tmp = np.append(x3_tmp, z)
            tmp = np.append(tmp, z)

        x1_itc[idx_itc_obs[c] - 1, :] = x1_tmp
        x2_itc[idx_itc_obs[c] - 1, :] = x2_tmp
        x3_itc[idx_itc_obs[c] - 1, :] = x3_tmp
        y1_itc[idx_itc_obs[c] - 1, :] = tmp
        c = c + 1

    if W_itc_obs < W:
        padd = int(W - W_itc_obs)
        for i in range(0, idx_itc_mis.shape[0]):
            idx_itc_obs = np.insert(idx_itc_obs, idx_itc_mis[i] - 1, 0)
        # idx_itc_mis = np.append(idx_itc_mis, np.zeros([W - padd, 1])).astype(int)
    else:
        idx_itc_mis = np.zeros(W).astype(int)

    x1_itc_all[sub, :, :] = x1_itc
    x2_itc_all[sub, :, :] = x2_itc
    x3_itc_all[sub, :, :] = x3_itc
    y1_itc_all[sub, :, :] = y1_itc
    Tr_itc_all[sub, :] = Tr_itc
    W_itc_obs_all[sub] = W_itc_obs
    W_itc_mis_all[sub] = W_itc_mis
    idx_itc_obs_all[sub, :] = idx_itc_obs

    # lt
    ######
    data_lt_tmp = data_lt_all_weeks.loc[(data_lt_all_weeks['subjectId'] == s) & (data_lt_all_weeks['low_p'] > 0) & (
            data_lt_all_weeks[
                'top_p'] > 0)]  # Here we get rid of "certain" trials where there is not stochasticity, as these will cause "initialization errors" in Stan due to the variance standardization

    W_lt_obs = data_lt_tmp['weekId'].unique().shape[0]
    W_lt_mis = W - W_lt_obs

    idx_lt_obs = data_lt_tmp['weekId'].unique()
    idx_lt_mis = np.setdiff1d(np.array(range(1, W + 1)), idx_lt_obs)

    tr_tmp = np.zeros([W])
    c = 0
    for w in data_lt_tmp['weekId'].unique():
        tr_tmp[idx_lt_obs[c] - 1] = np.shape(data_lt_tmp.loc[data_lt_tmp['weekId'] == w])[0]
        c = c + 1

    Tr_lt = np.transpose(tr_tmp.astype(int)).reshape(W)

    hi_p = np.zeros([W, Tr_max_lt_all])
    hi_narr = np.zeros([W, Tr_max_lt_all])
    hi_wide = np.zeros([W, Tr_max_lt_all])
    lo_narr = np.zeros([W, Tr_max_lt_all])
    lo_wide = np.zeros([W, Tr_max_lt_all])
    choice_lt = np.zeros([W, Tr_max_lt_all])

    c = 0
    for w in data_lt_tmp['weekId'].unique():
        tr_tmp = np.shape(data_lt_tmp.loc[data_lt_tmp['weekId'] == w])[0]
        hi_p_tmp = data_lt_tmp['top_p'].loc[data_lt_tmp['weekId'] == w]
        hi_narr_tmp = data_lt_tmp['hi_narr'].loc[data_lt_tmp['weekId'] == w]
        hi_wide_tmp = data_lt_tmp['hi_wide'].loc[data_lt_tmp['weekId'] == w]
        lo_narr_tmp = data_lt_tmp['lo_narr'].loc[data_lt_tmp['weekId'] == w]
        lo_wide_tmp = data_lt_tmp['lo_wide'].loc[data_lt_tmp['weekId'] == w]
        choice_tmp = data_lt_tmp['choice'].loc[data_lt_tmp['weekId'] == w]

        if tr_tmp < Tr_max_lt_all:
            padd = int(Tr_max_lt_all - tr_tmp)
            z = np.zeros([padd, 1])
            hi_p_tmp = np.append(hi_p_tmp, z)
            hi_narr_tmp = np.append(hi_narr_tmp, z)
            hi_wide_tmp = np.append(hi_wide_tmp, z)
            lo_narr_tmp = np.append(lo_narr_tmp, z)
            lo_wide_tmp = np.append(lo_wide_tmp, z)
            choice_tmp = np.append(choice_tmp, z)

        hi_p[idx_lt_obs[c] - 1, :] = hi_p_tmp
        hi_narr[idx_lt_obs[c] - 1, :] = hi_narr_tmp
        hi_wide[idx_lt_obs[c] - 1, :] = hi_wide_tmp
        lo_narr[idx_lt_obs[c] - 1, :] = lo_narr_tmp
        lo_wide[idx_lt_obs[c] - 1, :] = lo_wide_tmp
        choice_lt[idx_lt_obs[c] - 1, :] = choice_tmp
        c = c + 1

    if W_lt_obs < W:
        padd = int(W - W_lt_obs)
        for i in range(0, idx_lt_mis.shape[0]):
            idx_lt_obs = np.insert(idx_lt_obs, idx_lt_mis[i] - 1, 0)
        # idx_lt_mis = np.append(idx_lt_mis, np.zeros([W - padd, 1])).astype(int)
    else:
        idx_lt_mis = np.zeros(W).astype(int)

    # Add subject's (zero-padded) data to the data of all subjects
    hi_p_all[sub, :, :] = hi_p
    hi_narr_all[sub, :, :] = hi_narr
    hi_wide_all[sub, :, :] = hi_wide
    lo_narr_all[sub, :, :] = lo_narr
    lo_wide_all[sub, :, :] = lo_wide
    choice_lt_all[sub, :, :] = choice_lt
    Tr_lt_all[sub, :] = Tr_lt
    W_lt_obs_all[sub] = W_lt_obs
    W_lt_mis_all[sub] = W_lt_mis
    idx_lt_obs_all[sub, :] = idx_lt_obs
    # idx_lt_mis_all[sub, :] = idx_lt_mis

    sub = sub + 1

lg_data_stan = {'W': W,
                'N': N,

                'P_itc': P_itc,
                'idx_itc_obs': idx_itc_obs_all.astype(int),
                'Tr_max_itc': Tr_max_itc_all,
                'Tr_itc': Tr_itc_all.astype(int),
                'amount_later': x2_itc_all.tolist(),
                'amount_sooner': x3_itc_all.tolist(),
                'delay_later': x1_itc_all.tolist(),
                'choice_itc': y1_itc_all.astype(int).tolist(),

                'P_lt': P_lt,
                'idx_lt_obs': idx_lt_obs_all.astype(int),
                'Tr_max_lt': Tr_max_lt_all,
                'Tr_lt': Tr_lt_all.astype(int),
                'hi_p_lt': hi_p_all.astype(float) / 100,
                'hi_narr_lt': hi_narr_all.astype(float) / 1000,
                'hi_wide_lt': hi_wide_all.astype(float) / 1000,
                'lo_narr_lt': lo_narr_all.astype(float) / 1000,
                'lo_wide_lt': lo_wide_all.astype(float) / 1000,
                'choice_lt': choice_lt_all.astype(int),
                }

fit = model.sample(data=lg_data_stan, inits=0, chains=4, parallel_chains=2, threads_per_chain=16,
                   max_treedepth=10, adapt_delta=0.95, step_size=0.05, iter_warmup=1000, save_warmup=False,
                   iter_sampling=1000,
                   thin=1, show_progress=True,
                   output_dir='/n/gershman_lab/users/Stan/4_hierarchical/lt_hierarchical_parallel/with_default_sampler/',
                   save_diagnostics=False)