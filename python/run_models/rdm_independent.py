from cmdstanpy import cmdstan_path, CmdStanModel, set_cmdstan_path
import numpy as np
import pandas as pd
from scipy import stats
import cmdstanpy

model = "independent"  # "independent" or "dynamical"
normalize_exogenous_variables = True

if model == "independent":
    model = CmdStanModel(
        stan_file='/n/gershman_lab/users/Stan/1_code/publish_code/stan/final_independent_models/rdm_independent.stan',
                         cpp_options={'STAN_THREADS': True})

elif model == "dynamical":
        model = CmdStanModel(
            stan_file='/n/gershman_lab/users/Stan/1_code/publish_code/stan/final_dynamical_models/rdm_dynamic.stan',
            cpp_options={'STAN_THREADS': True})

# make data ITC
data_itc_all_weeks = pd.read_csv(
    '/n/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/itc_data_for_stan_90s.csv')

# make data RDM
data_rdm_all_weeks = pd.read_csv(
    '/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/rdm_data_for_stan_90s.csv')

# load exo data for the dynamic model
data_exo_all_weeks1 = np.array(pd.read_csv('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/2_task_data_files/'
                                  'data_90_subs/exogenous_variables/exo_day1_gng_rdm_PC1_valence.csv', header=None))
data_exo_all_weeks2 = np.array(pd.read_csv('/net/rcstorenfs02/ifs/rc_labs/gershman_lab/users/Stan/2_task_data_files/'
                                  'data_90_subs/exogenous_variables/exo_day1_gng_rdm_PC2_arousal.csv', header=None))

all_subjects =  data_rdm_all_weeks.subjectId
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

Nu_max_rdm_all = 356  # across all subjects
Nl_max_rdm_all = 208  # across all subjects
P_rdm = 3  # alpha (decision boundary), delta (drift rate), tau (non-decision time)
Nu_rdm_all = np.zeros([N, W])
Nl_rdm_all = np.zeros([N, W])
RTu_rdm_all = np.zeros([N, W, Nu_max_rdm_all])
RTl_rdm_all = np.zeros([N, W, Nl_max_rdm_all])
Cohu_rdm_all = np.zeros([N, W, Nu_max_rdm_all])
Cohl_rdm_all = np.zeros([N, W, Nl_max_rdm_all])
W_rdm_obs_all = np.zeros([N])
W_rdm_mis_all = np.zeros([N])
idx_rdm_obs_all = np.zeros([N, W])
idx_rdm_mis_all = np.zeros([N, W])
minRT_rdm = np.zeros([N, W])

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

    # RDM
    ######
    data_rdm_tmp = data_rdm_all_weeks.loc[data_rdm_all_weeks['subjectId'] == s]

    good_trials = data_rdm_tmp['key_choice'] != 'no_choice'
    data_rdm_tmp = data_rdm_tmp[good_trials]
    good_rt = data_rdm_tmp['rt'] > 200
    data_rdm_tmp = data_rdm_tmp[good_rt]

    good_rt = data_rdm_tmp['rt'] < 1500
    data_rdm_tmp = data_rdm_tmp[good_rt]

    W_rdm_obs = data_rdm_tmp['weekId'].unique().shape[0]
    W_rdm_mis = W - W_rdm_obs

    idx_rdm_obs = data_rdm_tmp['weekId'].unique()
    idx_rdm_mis = np.setdiff1d(np.array(range(1, W + 1)), idx_rdm_obs)

    c = 0
    nu = np.zeros([W])
    nl = np.zeros([W])
    min_rt = np.zeros([W])

    for w in data_rdm_tmp['weekId'].unique():
        tmp = data_rdm_tmp['correct'].loc[data_rdm_tmp['weekId'] == w]
        nu[idx_rdm_obs[c] - 1] = sum(tmp == 1)
        nl[idx_rdm_obs[c] - 1] = sum(tmp == 0)
        min_rt[idx_rdm_obs[c] - 1] = np.min(data_rdm_tmp['rt'].loc[data_rdm_tmp['weekId'] == w]) / 1000
        c = c + 1

    Nu_max = int(np.max(nu))
    Nl_max = int(np.max(nl))
    Nu = np.transpose(nu.astype(int)).reshape(W)
    Nl = np.transpose(nl.astype(int)).reshape(W)
    min_RT = np.transpose(min_rt).reshape(W)

    x1_rdm = np.zeros([W, Nu_max_rdm_all])
    x3_rdm = np.zeros([W, Nu_max_rdm_all])

    x2_rdm = np.zeros([W, Nl_max_rdm_all])
    x4_rdm = np.zeros([W, Nl_max_rdm_all])

    c = 0
    for w in data_rdm_tmp['weekId'].unique():
        tmp_d = data_rdm_tmp.loc[data_rdm_tmp['weekId'] == w]
        # tmp_d['correct'] = tmp_d['correct'].fillna(0)
        tmp_rtu = np.array((tmp_d['rt'].loc[tmp_d['correct'] == 1]) / 1000)
        tmp_rtl = np.array((tmp_d['rt'].loc[tmp_d['correct'] == 0]) / 1000)
        tmp_cohu = np.array((tmp_d['coh'].loc[tmp_d['correct'] == 1]))
        tmp_cohl = np.array((tmp_d['coh'].loc[tmp_d['correct'] == 0]))

        if len(tmp_rtu) < Nu_max_rdm_all:
            padd = Nu_max_rdm_all - len(tmp_rtu)
            z = np.zeros([padd, 1])
            tmp_rtu = np.append(tmp_rtu, z)
            tmp_cohu = np.append(tmp_cohu, z)

        if len(tmp_rtl) < Nl_max_rdm_all:
            padd = Nl_max_rdm_all - len(tmp_rtl)
            z = np.zeros([padd, 1])
            tmp_rtl = np.append(tmp_rtl, z)
            tmp_cohl = np.append(tmp_cohl, z)

        x1_rdm[idx_rdm_obs[c] - 1, :] = tmp_rtu
        x2_rdm[idx_rdm_obs[c] - 1, :] = tmp_rtl
        x3_rdm[idx_rdm_obs[c] - 1, :] = tmp_cohu
        x4_rdm[idx_rdm_obs[c] - 1, :] = tmp_cohl
        c = c + 1

    if W_rdm_obs < W:
        padd = int(W - W_rdm_obs)
        min_RT[idx_rdm_mis - 1] = np.min(min_RT[np.nonzero(min_RT)])
        for i in range(0, idx_rdm_mis.shape[0]):
            idx_rdm_obs = np.insert(idx_rdm_obs, idx_rdm_mis[i] - 1, 0)
        #  idx_rdm_mis = np.append(idx_rdm_mis, np.zeros([W - padd, 1])).astype(int)
    else:
        idx_rdm_mis = np.zeros(W).astype(int)

    Nu_rdm_all[sub, :] = Nu
    Nl_rdm_all[sub, :] = Nl
    RTu_rdm_all[sub, :, :] = x1_rdm
    RTl_rdm_all[sub, :, :] = x2_rdm
    Cohu_rdm_all[sub, :, :] = x3_rdm
    Cohl_rdm_all[sub, :, :] = x4_rdm
    minRT_rdm[sub, :] = min_RT
    W_rdm_obs_all[sub] = W_rdm_obs
    W_rdm_mis_all[sub] = W_rdm_mis
    idx_rdm_obs_all[sub, :] = idx_rdm_obs

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

                'P_rdm': P_rdm,
                'idx_rdm_obs': idx_rdm_obs_all.astype(int).tolist(),
                'Nu_max_rdm': Nu_max_rdm_all,
                'Nl_max_rdm': Nl_max_rdm_all,
                'Nu_rdm': Nu_rdm_all.astype(int).tolist(),
                'Nl_rdm': Nl_rdm_all.astype(int).tolist(),
                'RTu_rdm': RTu_rdm_all.tolist(),
                'RTl_rdm': RTl_rdm_all.tolist(),
                'Cohu_rdm': (Cohu_rdm_all).tolist(),
                'Cohl_rdm': (Cohl_rdm_all).tolist(),
                'minRT_rdm': minRT_rdm.tolist(),
                'RTbound_rdm': 0.1,
                }


fit = model.sample(data=lg_data_stan, inits=0, chains=4, parallel_chains=2, threads_per_chain=16,
                   max_treedepth=10, adapt_delta=0.95, step_size=0.05, iter_warmup=1000, save_warmup=False,
                   iter_sampling=1000,
                   thin=1, show_progress=True,
                   output_dir='OUTPUTDIR',
                   save_diagnostics=False)

