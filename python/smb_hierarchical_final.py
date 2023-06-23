# Stan code was written by Daniel Reznik (reznikda@gmail.com) and Rahul Bhui (rbhui@mit.edu)
# Change detection model was originally implemented in Python by Hanna Hillman (hhillman231@gmail.com)
# Numerosity comparison model was originally written in Matlab by John Michael ()
# July 2021, Israel, Massachusetts, Connecticut, Saxony
##############################################################################

from cmdstanpy import cmdstan_path, CmdStanModel, set_cmdstan_path
import numpy as np
import pandas as pd
from scipy import stats
import cmdstanpy

# set_cmdstan_path('/ncf/gershman/Lab/cmdstan-2.24.0')

model = "independent"  # "independent" or "dynamical"

if model == "independent":
    model = CmdStanModel(
        stan_file='/n/gershman_lab/users/Stan/1_code/publish_code/stan/final_independent_models/smb_hierarchical_independent_with_gq_final.stan',
                         cpp_options={'STAN_THREADS': True})

elif model == "dynamical":
        model = CmdStanModel(
            stan_file='/n/gershman_lab/users/Stan/1_code/publish_code/stan/final_dynamical_models/smb_hierarchical_super_model_noWeights_noPrev_0effects_N01Priors_GQ_varExplained_PD_exo2.stan',
            cpp_options={'STAN_THREADS': True})

# make data ITC
data_itc_all_weeks = pd.read_csv('/n/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/itc_data_for_stan_90s.csv')

# make data SMB
data_smb_all_weeks = pd.read_csv('/n/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/smb_data_for_stan_90s_updatedRegressorsNov22_rescaled_signVoverTU.csv')


all_subjects = data_smb_all_weeks.subjectId
# make data frames
subjects = all_subjects.unique()

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

Tr_max_smb_all = 300
P_smb = 3
W_smb_obs_all = np.zeros([N])
W_smb_mis_all = np.zeros([N])
idx_smb_obs_all = np.zeros([N, W])
idx_smb_mis_all = np.zeros([N, W])
x_smb_all = np.zeros([N, W, Tr_max_smb_all, 3])
y_smb_all = np.zeros([N, W, Tr_max_smb_all])
Tr_smb_all = np.zeros([N, W])

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

    # SMB
    ######
    data_smb_tmp = data_smb_all_weeks.loc[data_smb_all_weeks['subjectId'] == s]

    W_smb_obs = data_smb_tmp['weekId'].unique().shape[0]
    W_smb_mis = W - W_smb_obs

    idx_smb_obs = data_smb_tmp['weekId'].unique()
    idx_smb_mis = np.setdiff1d(np.array(range(1, W + 1)), idx_smb_obs)

    tr_tmp = np.zeros([W])
    c = 0
    for w in data_smb_tmp['weekId'].unique():
        tr_tmp[idx_smb_obs[c] - 1] = np.shape(data_smb_tmp.loc[data_smb_tmp['weekId'] == w])[0]
        c = c + 1

    Tr_smb = np.transpose(tr_tmp.astype(int)).reshape(W)

    x_smb = np.zeros([W, Tr_max_smb_all, P_smb])
    y_smb = np.zeros([W, Tr_max_smb_all])

    c = 0
    for w in data_smb_tmp['weekId'].unique():
        tr_tmp = np.shape(data_smb_tmp['V'].loc[data_smb_tmp['weekId'] == w])[0]

        v_tmp = np.array(data_smb_tmp['V'].loc[data_smb_tmp['weekId'] == w])
        vtu_tmp = np.array(data_smb_tmp['VTU'].loc[data_smb_tmp['weekId'] == w])
        ru_tmp = np.array(data_smb_tmp['RU'].loc[data_smb_tmp['weekId'] == w])
        choice_tmp = data_smb_tmp['chosen_machine'].loc[data_smb_tmp['weekId'] == w]

        if tr_tmp < Tr_max_smb_all:
            padd = int(Tr_max_smb_all - tr_tmp)
            z = np.zeros([padd, 1])
            v_tmp = np.append(v_tmp, z)
            vtu_tmp = np.append(vtu_tmp, z)
            ru_tmp = np.append(ru_tmp, z)
            choice_tmp = np.append(choice_tmp, z)

        tmp = np.transpose(np.array([v_tmp, vtu_tmp, ru_tmp]))
        x_smb[idx_smb_obs[c] - 1, :, :] = tmp
        y_smb[idx_smb_obs[c] - 1, :] = choice_tmp
        c = c + 1

    if W_smb_obs < W:
        padd = int(W - W_smb_obs)
        for i in range(0, idx_smb_mis.shape[0]):
            idx_smb_obs = np.insert(idx_smb_obs, idx_smb_mis[i] - 1, 0)
        #  idx_smb_mis = np.append(idx_smb_mis, np.zeros([W - padd, 1])).astype(int)
    else:
        idx_smb_mis = np.zeros(W).astype(int)


    x_smb[:, :, 1] = x_smb[:, :, 1] / 100

    x_smb_all[sub, :, :, :] = x_smb
    y_smb_all[sub, :, :] = y_smb
    Tr_smb_all[sub, :] = Tr_smb
    W_smb_obs_all[sub] = W_smb_obs
    W_smb_mis_all[sub] = W_smb_mis
    idx_smb_obs_all[sub, :] = idx_smb_obs
    # idx_smb_mis_all[sub, :] = idx_smb_mis

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

                'P_smb': P_smb,
                'idx_smb_obs': idx_smb_obs_all.astype(int),
                'Tr_max_smb': Tr_max_smb_all,
                'Tr_smb': Tr_smb_all.astype(int),
                'x_smb': x_smb_all,
                'y_smb': y_smb_all.astype(int),
                }

fit = model.sample(data=lg_data_stan, inits=0, chains=4, parallel_chains=2, threads_per_chain=16,
                   max_treedepth=10, adapt_delta=0.95, step_size=0.05, iter_warmup=1000, save_warmup=False,
                   iter_sampling=1000,
                   thin=1, show_progress=True,
                   output_dir='/n/gershman_lab/users/Stan/4_hierarchical/smb_hierarchical_parallel/with_default_sampler',
                   save_diagnostics=False)