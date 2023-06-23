from cmdstanpy import cmdstan_path, CmdStanModel, set_cmdstan_path
import numpy as np
import pandas as pd
from scipy import stats
import cmdstanpy

# set_cmdstan_path('/ncf/gershman/Lab/cmdstan-2.24.0')

model = "independent"  # "independent" or "dynamical"
normalize_exogenous_variables = True

if model == "independent":
    model = CmdStanModel(
        stan_file='/n/gershman_lab/users/Stan/1_code/publish_code/stan/final_independent_models/nc_hierarchical_independent_with_gq_final.stan',
                         cpp_options={'STAN_THREADS': True})

elif model == "dynamical":
        model = CmdStanModel(
            stan_file='/n/gershman_lab/users/Stan/1_code/publish_code/stan/final_dynamical_models/nc_hierarchical_super_model_noWeights_noPrev_0effects_N01Priors_GQ_varExplained_PD_exo2.stan',
            cpp_options={'STAN_THREADS': True})

# make data ITC
data_itc_all_weeks = pd.read_csv('/n/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/itc_data_for_stan_90s.csv')

# make data NC
data_nc_all_weeks = pd.read_csv('/n/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/nc_data_for_stan_90s.csv')

# load exo data for the dynamic model
data_exo_all_weeks1 = np.array(pd.read_csv('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/2_task_data_files/'
                                  'data_90_subs/exogenous_variables/exo_day2_lt_nc_smb_PC1_valence.csv', header=None))
data_exo_all_weeks2 = np.array(pd.read_csv('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/2_task_data_files/'
                                  'data_90_subs/exogenous_variables/exo_day2_lt_nc_smb_PC2_arousal.csv', header=None))

all_subjects = data_nc_all_weeks.subjectId

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

Tr_max_nc_all = 160
P_nc = 1
mu_nc_all = np.zeros([N, W, Tr_max_nc_all])
sigma_nc_all = np.zeros([N, W, Tr_max_nc_all])
choice_nc_all = np.zeros([N, W, Tr_max_nc_all])
Tr_nc_all = np.zeros([N, W])
W_nc_obs_all = np.zeros([N])
W_nc_mis_all = np.zeros([N])
idx_nc_obs_all = np.zeros([N, W])
idx_nc_mis_all = np.zeros([N, W])

Num_exo = 2  # 2 PCs
exo_data_all = np.zeros([N, W, Num_exo])

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

    # NC
    #######
    data_nc_tmp = data_nc_all_weeks.loc[data_nc_all_weeks['subjectId'] == s]
    W_nc_obs = data_nc_tmp['weekId'].unique().shape[0]
    W_nc_mis = W - W_nc_obs

    idx_nc_obs = data_nc_tmp['weekId'].unique()
    idx_nc_mis = np.setdiff1d(np.array(range(1, W + 1)), idx_nc_obs)

    tr_tmp = np.zeros([W])
    Tr_nc = np.zeros([W])
    c = 0
    for w in data_nc_tmp['weekId'].unique():
        tr_tmp[idx_nc_obs[c] - 1] = np.shape(data_nc_tmp.loc[data_nc_tmp['weekId'] == w])[0]
        c = c + 1

    T_max_nc = np.max(tr_tmp).astype(int)
    Tr_nc = np.transpose(tr_tmp.astype(int)).reshape(W)

    # data_nc_tmp['key_choice'] = data_nc_tmp['key_choice'].replace({'left': -1, 'right': 1})
    mu = np.zeros([W, Tr_max_nc_all])
    sigma = np.zeros([W, Tr_max_nc_all])
    choice_nc = np.zeros([W, Tr_max_nc_all])

    tr_tmp = np.zeros([W])
    c = 0
    for w in data_nc_tmp['weekId'].unique():
        tr_tmp = np.shape(data_nc_tmp.loc[data_nc_tmp['weekId'] == w])[0]
        mu_tmp = data_nc_tmp['left_circles'].loc[data_nc_tmp['weekId'] == w] - \
                 data_nc_tmp['right_circles'].loc[data_nc_tmp['weekId'] == w]
        sigma_tmp = np.array(np.sqrt(data_nc_tmp['left_circles'].loc[data_nc_tmp['weekId'] == w] ** 2 +
                                     data_nc_tmp['right_circles'].loc[data_nc_tmp['weekId'] == w] ** 2))
        choice_tmp = data_nc_tmp['key_choice'].loc[data_nc_tmp['weekId'] == w]

        # mu_tmp[choice_tmp == -1] = -mu_tmp[choice_tmp == -1]

        if tr_tmp < Tr_max_nc_all:
            padd = int(Tr_max_nc_all - tr_tmp)
            z = np.zeros([padd, 1])
            mu_tmp = np.append(mu_tmp, z)
            sigma_tmp = np.append(sigma_tmp, z)
            choice_tmp = np.append(choice_tmp, z)

        mu[idx_nc_obs[c] - 1, :] = mu_tmp
        sigma[idx_nc_obs[c] - 1, :] = sigma_tmp
        choice_nc[idx_nc_obs[c] - 1, :] = choice_tmp
        c = c + 1

    if W_nc_obs < W:
        padd = int(W - W_nc_obs)
        for i in range(0, idx_nc_mis.shape[0]):
            idx_nc_obs = np.insert(idx_nc_obs, idx_nc_mis[i] - 1, 0)
        # idx_nc_mis = np.append(idx_nc_mis, np.zeros([W - padd, 1])).astype(int)
    else:
        idx_nc_mis = np.zeros(W).astype(int)

    mu_nc_all[sub, :, :] = mu
    sigma_nc_all[sub, :, :] = sigma
    choice_nc_all[sub, :, :] = choice_nc
    Tr_nc_all[sub, :] = Tr_nc
    W_nc_obs_all[sub] = W_nc_obs
    W_nc_mis_all[sub] = W_nc_mis
    idx_nc_obs_all[sub, :] = idx_nc_obs
    #  idx_nc_mis_all[sub, :] = idx_nc_mis

    # EXO
    ######
    data_exo_tmp1 = data_exo_all_weeks1[sub, :]
    data_exo_tmp2 = data_exo_all_weeks2[sub, :]

    idx_exo_obs1 = np.array([i for i, x in enumerate(~np.isnan(data_exo_tmp1)) if x]) + 1
    idx_exo_mis1 = np.array([i for i, x in enumerate(np.isnan(data_exo_tmp1)) if x]) + 1

    idx_exo_obs2 = np.array([i for i, x in enumerate(~np.isnan(data_exo_tmp2)) if x]) + 1
    idx_exo_mis2 = np.array([i for i, x in enumerate(np.isnan(data_exo_tmp2)) if x]) + 1

    exo_data_inter1 = np.append(data_exo_tmp1[idx_exo_obs1 - 1],
                                np.interp(idx_exo_mis1 - 1, idx_exo_obs1 - 1, data_exo_tmp1[idx_exo_obs1 - 1]))
    exo_data_inter2 = np.append(data_exo_tmp2[idx_exo_obs2 - 1],
                                np.interp(idx_exo_mis2 - 1, idx_exo_obs2 - 1, data_exo_tmp2[idx_exo_obs2 - 1]))

    exo_data_all[sub, :, 0] = exo_data_inter1
    exo_data_all[sub, :, 1] = exo_data_inter2

    sub = sub + 1


if normalize_exogenous_variables:
    for s in range(0, N):
        for d in range(0, Num_exo):
            if np.unique(exo_data_all[s, :, d]).shape[
                0] == 1:  # If the values are constant throughout sessions, convert to all zeros, to avoid getting nan's
                exo_data_all[s, :, d] = exo_data_all[s, :, d] * 0
            else:
                tmp_min = exo_data_all[s, :, d].min()
                tmp_max = exo_data_all[s, :, d].max()
                tmp = (exo_data_all[s, :, d] - tmp_min) / (tmp_max - tmp_min)  # This is normalized to [0,1]
                exo_data_all[s, :, d] = tmp * 2 - 1  # This is normalized to [-1,1]]


lg_data_stan = {'W': W,
                'N': N,

                'exo_q_num': Num_exo,
                'U': exo_data_all.tolist(),

                'P_itc': P_itc,
                'idx_itc_obs': idx_itc_obs_all.astype(int),
                'Tr_max_itc': Tr_max_itc_all,
                'Tr_itc': Tr_itc_all.astype(int),
                'amount_later': x2_itc_all.tolist(),
                'amount_sooner': x3_itc_all.tolist(),
                'delay_later': x1_itc_all.tolist(),
                'choice_itc': y1_itc_all.astype(int).tolist(),

                'P_nc': P_nc,
                'W_nc_obs': W_nc_obs_all.astype(int),
                'idx_nc_obs': idx_nc_obs_all.astype(int),
                'Tr_max_nc': Tr_max_nc_all,
                'Tr_nc': Tr_nc_all.astype(int),
                'deltaM': mu_nc_all.astype(int).tolist(),
                'TotalS': sigma_nc_all.tolist(),
                'choice_nc': choice_nc_all.astype(int).tolist(),
                }

fit = model.sample(data=lg_data_stan, inits=0, chains=4, parallel_chains=2, threads_per_chain=16,
                   max_treedepth=10, adapt_delta=0.95, step_size=0.05, iter_warmup=1000, save_warmup=False,
                   iter_sampling=1000,
                   thin=1, show_progress=True,
                   output_dir='OUTPUTDIR',
                   save_diagnostics=False)