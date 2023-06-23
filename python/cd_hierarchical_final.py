from cmdstanpy import cmdstan_path, CmdStanModel, set_cmdstan_path
import numpy as np
import pandas as pd
from scipy import stats
import cmdstanpy

model = "independent"  # "independent" or "dynamical"

if model == "independent":
    model = CmdStanModel(
        stan_file='/n/gershman_lab/users/Stan/1_code/publish_code/stan/final_independent_models/cd_hierarcahical_independent_with_gq_final.stan',
                         cpp_options={'STAN_THREADS': True})

elif model == "dynamical":
        model = CmdStanModel(
            stan_file='/n/gershman_lab/users/Stan/1_code/publish_code/stan/final_dynamical_models/cd_hierarchical_super_model_noWeights_noPrev_0effects_N01Priors_GQ_varExplained_PD_exo2.stan',
            cpp_options={'STAN_THREADS': True})


# make data ITC
data_itc_all_weeks = pd.read_csv(
    '/n/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/itc_data_for_stan_90s.csv')

# make data CD
data_cd_all_weeks = pd.read_csv(
    '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/cd_data_for_stan_90s_with_rt.csv')

all_subjects =  data_cd_all_weeks.subjectId

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

Tr_max_cd_all = 200

P_cd = 4  # This is 2x2 (separate criterion and sigma for each of the 4 blocks)

W_cd_obs_all = np.zeros([N])
W_cd_mis_all = np.zeros([N])
idx_cd_obs_all = np.zeros([N, W])
idx_cd_mis_all = np.zeros([N, W])
Tr_cd_all = np.zeros([N, W])
Tar_cd_all = np.zeros([N, W, Tr_max_cd_all])
Nb_cd_all = np.zeros([N, W, Tr_max_cd_all])
delta_cd_all = np.zeros([N, W, Tr_max_cd_all])
choice_cd_all = np.zeros([N, W, Tr_max_cd_all])

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

    # CD
    ######
    data_cd_tmp = data_cd_all_weeks.loc[data_cd_all_weeks['subjectId'] == s]

    W_cd_obs = data_cd_tmp['weekId'].unique().shape[0]
    W_cd_mis = W - W_cd_obs

    idx_cd_obs = data_cd_tmp['weekId'].unique()
    idx_cd_mis = np.setdiff1d(np.array(range(1, W + 1)), idx_cd_obs)

    tr_tmp = np.zeros([W])
    c = 0
    for w in data_cd_tmp['weekId'].unique():
        tr_tmp[idx_cd_obs[c] - 1] = np.shape(data_cd_tmp.loc[data_cd_tmp['weekId'] == w])[0]
        c = c + 1

    Tr_cd = np.transpose(tr_tmp.astype(int)).reshape(W)

    Nb_cd = np.zeros([W, Tr_max_cd_all])
    Tar_cd = np.zeros([W, Tr_max_cd_all])
    Delta_cd = np.zeros([W, Tr_max_cd_all])
    choice_cd = np.zeros([W, Tr_max_cd_all])

    c = 0
    for w in data_cd_tmp['weekId'].unique():
        tr_tmp = np.shape(data_cd_tmp.loc[data_cd_tmp['weekId'] == w])[0]
        nb_tmp = data_cd_tmp['n'].loc[data_cd_tmp['weekId'] == w]
        tar_tmp = data_cd_tmp['t'].loc[data_cd_tmp['weekId'] == w]
        delta_tmp = data_cd_tmp['d'].loc[data_cd_tmp['weekId'] == w]
        choice_tmp = data_cd_tmp['behav'].loc[data_cd_tmp['weekId'] == w]

        if tr_tmp < Tr_max_cd_all:
            padd = int(Tr_max_cd_all - tr_tmp)
            z = np.zeros([padd, 1])
            nb_tmp = np.append(nb_tmp, z)
            tar_tmp = np.append(tar_tmp, z)
            delta_tmp = np.append(delta_tmp, z)
            choice_tmp = np.append(choice_tmp, z)

        Nb_cd[idx_cd_obs[c] - 1, :] = nb_tmp
        Tar_cd[idx_cd_obs[c] - 1, :] = tar_tmp
        Delta_cd[idx_cd_obs[c] - 1, :] = delta_tmp
        choice_cd[idx_cd_obs[c] - 1, :] = choice_tmp
        c = c + 1

    if W_cd_obs < W:
        padd = int(W - W_cd_obs)
        for i in range(0, idx_cd_mis.shape[0]):
            idx_cd_obs = np.insert(idx_cd_obs, idx_cd_mis[i] - 1, 0)
        # idx_cd_mis = np.append(idx_cd_mis, np.zeros([W - padd, 1])).astype(int)
    else:
        idx_cd_mis = np.zeros(W).astype(int)

    Nb_cd_all[sub, :, :] = Nb_cd
    Tar_cd_all[sub, :, :] = Tar_cd
    delta_cd_all[sub, :, :] = Delta_cd
    choice_cd_all[sub, :, :] = choice_cd
    Tr_cd_all[sub, :] = Tr_cd
    W_cd_obs_all[sub] = W_cd_obs
    W_cd_mis_all[sub] = W_cd_mis
    idx_cd_obs_all[sub, :] = idx_cd_obs
    # idx_cd_mis_all[sub, :] = idx_cd_mis

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

                'P_cd': P_cd,
                'idx_cd_obs': idx_cd_obs_all.astype(int),
                'Tr_max_cd': Tr_max_cd_all,
                'Tr_cd': Tr_cd_all.astype(int),
                'Nb_cd': Nb_cd_all.astype(int),
                'Tar_cd': Tar_cd_all.astype(int),
                'delta_cd': delta_cd_all / 100,
                'choice_cd': choice_cd_all.astype(int),
                }


fit = model.sample(data=lg_data_stan, inits=0, chains=4, parallel_chains=2, threads_per_chain=16,
                   max_treedepth=10, adapt_delta=0.95, step_size=0.05, iter_warmup=1000, save_warmup=False,
                   iter_sampling=1000,
                   thin=1, show_progress=True,
                   output_dir='/n/gershman_lab/users/Stan/4_hierarchical/cd_hierarchical_parallel/with_default_sampler/',
                   save_diagnostics=False)