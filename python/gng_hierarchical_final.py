from cmdstanpy import cmdstan_path, CmdStanModel, set_cmdstan_path
import numpy as np
import pandas as pd
from scipy import stats
import cmdstanpy

model = "independent"  # "independent" or "dynamical"

if model == "independent":
    model = CmdStanModel(
        stan_file='/n/gershman_lab/users/Stan/1_code/publish_code/stan/final_independent_models/gng_hierarchical_independent_modified_model_with_gq_final.stan',
                         cpp_options={'STAN_THREADS': True})

elif model == "dynamical":
        model = CmdStanModel(
            stan_file='/n/gershman_lab/users/Stan/1_code/publish_code/stan/final_dynamical_models/gng_hierarchical_super_model_noWeights_noPrev_0effects_rho2_RewPun_Neut_GQ_varExplained_PD_exo2_piPos.stan',
            cpp_options={'STAN_THREADS': True})


# make data ITC - this is required only for reduce_sum sliced argument
data_itc_all_weeks = pd.read_csv('/n/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/itc_data_for_stan_90s.csv')
#data_itc_all_weeks['choice'] = data_itc_all_weeks['choice'].replace({'now': 0, 'later': 1})

# make data GNG
data_gng_all_weeks = pd.read_csv('/n/gershman_lab/users/Stan/2_task_data_files/data_90_subs/final/gng_data_for_stan_90s.csv')

all_subjects = data_gng_all_weeks.subjectId # make data frames
subjects = all_subjects.unique()
# exclude_gng_subjects_with_accuracy_below_055:
sub_exclude = [1, 2, 3, 6, 8, 9, 15, 16, 17, 28, 33, 34, 37, 39, 44, 59, 60, 68, 70, 71, 74, 76, 78, 86]

for idx in range(0, np.shape(sub_exclude)[0]):
    data_gng_all_weeks = data_gng_all_weeks.drop(
        data_gng_all_weeks[(data_gng_all_weeks['subjectId'] == subjects[sub_exclude[idx] - 1])].index)

all_subjects = data_gng_all_weeks.subjectId # make data frames
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

Tr_max_gng_all = 80
Bl = 3  # number of blocks in gng
P_gng = 5
W_gng_obs_all = np.zeros([N])
W_gng_mis_all = np.zeros([N])
idx_gng_obs_all = np.zeros([N, W])
idx_gng_mis_all = np.zeros([N, W])
Tr_gng_all = np.zeros([N, W, Bl])
cue_gng_all = np.zeros([N, W, Bl, Tr_max_gng_all])
pressed_gng_all = np.zeros([N, W, Bl, Tr_max_gng_all])
outcome_gng_all = np.zeros([N, W, Bl, Tr_max_gng_all])

sub = 0
for s in subjects[0:N]:

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

    # GNG
    ######
    data_gng_tmp = data_gng_all_weeks.loc[data_gng_all_weeks['subjectId'] == s]

    W_gng_obs = data_gng_tmp['weekId'].unique().shape[0]
    W_gng_mis = W - W_gng_obs

    idx_gng_obs = data_gng_tmp['weekId'].unique()
    idx_gng_mis = np.setdiff1d(np.array(range(1, W + 1)), idx_gng_obs)

    tr_tmp = np.zeros([W, Bl])
    Tr = np.zeros([W, Bl])
    c = 0
    for w in data_gng_tmp['weekId'].unique():
        c1 = 0
        for bl in data_gng_tmp['block'].unique():
            tr_tmp[idx_gng_obs[c] - 1, c1] = \
                np.shape(data_gng_tmp.loc[(data_gng_tmp['weekId'] == w) & (data_gng_tmp['block'] == bl)])[0]
            c1 = c1 + 1
        c = c + 1

    # Tr_gng = np.transpose(tr_tmp.astype(int)).reshape(W, Bl)
    Tr_gng = tr_tmp

    c = 0
    tr_tmp = np.zeros([W, Bl])
    cue_gng = np.zeros([W, Bl, Tr_max_gng_all])
    pressed_gng = np.zeros([W, Bl, Tr_max_gng_all])
    outcome_gng = np.zeros([W, Bl, Tr_max_gng_all])

    for w in data_gng_tmp['weekId'].unique():
        c1 = 0
        for bl in data_gng_tmp['block'].unique():
            tr_tmp = np.shape(data_gng_tmp['cond'].loc[(data_gng_tmp['weekId'] == w) & (data_gng_tmp['block'] == bl)])[
                0]

            cue_tmp = data_gng_tmp['cond'].loc[(data_gng_tmp['weekId'] == w) & (data_gng_tmp['block'] == bl)]
            outcome_tmp = data_gng_tmp['feedback_points'].loc[
                (data_gng_tmp['weekId'] == w) & (data_gng_tmp['block'] == bl)]
            pressed_tmp = data_gng_tmp['choice'].loc[(data_gng_tmp['weekId'] == w) & (data_gng_tmp['block'] == bl)]

            if tr_tmp < Tr_max_gng_all:
                padd = int(Tr_max_gng_all - tr_tmp)
                z = np.zeros([padd, 1])
                cue_tmp = np.append(cue_tmp, z)
                outcome_tmp = np.append(outcome_tmp, z)
                pressed_tmp = np.append(pressed_tmp, z)

            cue_gng[idx_gng_obs[c] - 1, c1, :] = cue_tmp
            outcome_gng[idx_gng_obs[c] - 1, c1, :] = outcome_tmp
            pressed_gng[idx_gng_obs[c] - 1, c1, :] = pressed_tmp
            c1 = c1 + 1
        c = c + 1

    if W_gng_obs < W:
        padd = int(W - W_gng_obs)
        for i in range(0, idx_gng_mis.shape[0]):
            idx_gng_obs = np.insert(idx_gng_obs, idx_gng_mis[i] - 1, 0)
        # idx_gng_mis = np.append(idx_gng_mis, np.zeros([W - padd, 1])).astype(int)
    else:
        idx_gng_mis = np.zeros(W).astype(int)

    Tr_gng_all[sub, :, :] = Tr_gng
    cue_gng_all[sub, :, :, :] = cue_gng
    pressed_gng_all[sub, :, :, :] = pressed_gng
    outcome_gng_all[sub, :, :, :] = outcome_gng
    W_gng_obs_all[sub] = W_gng_obs
    W_gng_mis_all[sub] = W_gng_mis
    idx_gng_obs_all[sub, :] = idx_gng_obs
    # idx_gng_mis_all[sub, :] = idx_gng_mis

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

                'P_gng': P_gng,
                'idx_gng_obs': idx_gng_obs_all.astype(int),
                'Bl': Bl,
                'Tr_gng': Tr_gng_all.astype(int),
                'Tr_max_gng': Tr_max_gng_all,
                'cue_gng': cue_gng_all.astype(int),
                'pressed_gng': pressed_gng_all.astype(int),
                'outcome_gng': outcome_gng_all
                }


fit = model.sample(data=lg_data_stan, inits=0, chains=4, parallel_chains=2, threads_per_chain=16,
                   max_treedepth=10, adapt_delta=0.95, step_size=0.05, iter_warmup=1000, save_warmup=False,
                   iter_sampling=1000,
                   thin=1, show_progress=True,
                   output_dir='/n/gershman_lab/users/Stan/4_hierarchical/gng_hierarchical_parallel/with_default_sampler/',
                   save_diagnostics=False)