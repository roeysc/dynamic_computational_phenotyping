import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# CD
cd_predicted = pd.read_csv('/Users/daniel/PycharmProjects/comp_pheno_cmdStan/results/hierarchical/tasks/cd_hierarchical/final/cd_predicted_independent_norm_priors.csv')
cd_real = pd.read_csv('/Users/daniel/PycharmProjects/comp_pheno_cmdStan/cd_data_for_stan_90s_with_rt.csv')

cd_predicted = cd_predicted[cd_predicted['0'] != -2]
subjects = np.unique(cd_real.subjectId)

cd_predicted_reset = cd_predicted.reset_index()

pred = []
for t in range(0, np.shape(cd_predicted)[0]):
    tmp = cd_predicted_reset['0'][t]
    pred = np.append(pred, tmp)

cd_real['pred'] = pred

## ALL WEEKS - STD ACROSS ALL SUBJECTS
cd_real_s1 = cd_real.loc[cd_real['t'] != 0]
cd_real_s1 = cd_real_s1.sort_values(by=['n'])

cd_real_s3 = cd_real.loc[cd_real['t'] == 0]
cd_real_s3 = cd_real_s3.sort_values(by=['n'])

subjects = cd_real_s1.subjectId.unique()
predicted_cd_s2_n3 = np.zeros([np.shape(subjects)[0]])
real_cd_s2_n3 = np.zeros([np.shape(subjects)[0]])

predicted_cd_s2_n4 = np.zeros([np.shape(subjects)[0]])
real_cd_s2_n4 = np.zeros([np.shape(subjects)[0]])

predicted_cd_s2_n6 = np.zeros([np.shape(subjects)[0]])
real_cd_s2_n6 = np.zeros([np.shape(subjects)[0]])

predicted_cd_s2_n8 = np.zeros([np.shape(subjects)[0]])
real_cd_s2_n8 = np.zeros([np.shape(subjects)[0]])

sub = 0
for s in subjects:
    predicted_cd_s2_n3[sub] = np.mean(cd_real_s1['pred'].loc[(cd_real_s1['n'] == 3) & (cd_real_s1['subjectId'] == s)])
    real_cd_s2_n3[sub] = np.mean(cd_real_s1['behav'].loc[(cd_real_s1['n'] == 3) & (cd_real_s1['subjectId'] == s)])

    predicted_cd_s2_n4[sub] = np.mean(cd_real_s1['pred'].loc[(cd_real_s1['n'] == 4) & (cd_real_s1['subjectId'] == s)])
    real_cd_s2_n4[sub] = np.mean(cd_real_s1['behav'].loc[(cd_real_s1['n'] == 4) & (cd_real_s1['subjectId'] == s)])

    predicted_cd_s2_n6[sub] = np.mean(cd_real_s1['pred'].loc[(cd_real_s1['n'] == 6) & (cd_real_s1['subjectId'] == s)])
    real_cd_s2_n6[sub] = np.mean(cd_real_s1['behav'].loc[(cd_real_s1['n'] == 6) & (cd_real_s1['subjectId'] == s)])

    predicted_cd_s2_n8[sub] = np.mean(cd_real_s1['pred'].loc[(cd_real_s1['n'] == 8) & (cd_real_s1['subjectId'] == s) & (cd_real_s1['t'] == 1)])
    real_cd_s2_n8[sub] = np.mean(cd_real_s1['behav'].loc[(cd_real_s1['n'] == 8) & (cd_real_s1['subjectId'] == s) & (cd_real_s1['t'] == 1)])

    sub = sub + 1

pred_n = [predicted_cd_s2_n3, predicted_cd_s2_n4, predicted_cd_s2_n6, predicted_cd_s2_n8]
real_n = [real_cd_s2_n3, real_cd_s2_n4, real_cd_s2_n6, real_cd_s2_n8]


subjects = cd_real_s3.subjectId.unique()
predicted_cd_s2_n3_nt = np.zeros([np.shape(subjects)[0]])
real_cd_s2_n3_nt = np.zeros([np.shape(subjects)[0]])

predicted_cd_s2_n4_nt = np.zeros([np.shape(subjects)[0]])
real_cd_s2_n4_nt = np.zeros([np.shape(subjects)[0]])

predicted_cd_s2_n6_nt = np.zeros([np.shape(subjects)[0]])
real_cd_s2_n6_nt = np.zeros([np.shape(subjects)[0]])

predicted_cd_s2_n8_nt = np.zeros([np.shape(subjects)[0]])
real_cd_s2_n8_nt = np.zeros([np.shape(subjects)[0]])

sub = 0
for s in subjects:
    predicted_cd_s2_n3_nt[sub] = np.mean(cd_real_s3['pred'].loc[(cd_real_s3['n'] == 3) & (cd_real_s3['subjectId'] == s)])
    real_cd_s2_n3_nt[sub] = np.mean(cd_real_s3['behav'].loc[(cd_real_s3['n'] == 3) & (cd_real_s3['subjectId'] == s)])

    predicted_cd_s2_n4_nt[sub] = np.mean(cd_real_s3['pred'].loc[(cd_real_s3['n'] == 4) & (cd_real_s3['subjectId'] == s)])
    real_cd_s2_n4_nt[sub] = np.mean(cd_real_s3['behav'].loc[(cd_real_s3['n'] == 4) & (cd_real_s3['subjectId'] == s)])

    predicted_cd_s2_n6_nt[sub] = np.mean(cd_real_s3['pred'].loc[(cd_real_s3['n'] == 6) & (cd_real_s3['subjectId'] == s)])
    real_cd_s2_n6_nt[sub] = np.mean(cd_real_s3['behav'].loc[(cd_real_s3['n'] == 6) & (cd_real_s3['subjectId'] == s)])

    predicted_cd_s2_n8_nt[sub] = np.mean(cd_real_s3['pred'].loc[(cd_real_s3['n'] == 8) & (cd_real_s3['subjectId'] == s)])
    real_cd_s2_n8_nt[sub] = np.mean(cd_real_s3['behav'].loc[(cd_real_s3['n'] == 8) & (cd_real_s3['subjectId'] == s)])

    sub = sub + 1


pred_n_nt = [predicted_cd_s2_n3_nt, predicted_cd_s2_n4_nt, predicted_cd_s2_n6_nt, predicted_cd_s2_n8_nt]
real_n_nt = [real_cd_s2_n3_nt, real_cd_s2_n4_nt, real_cd_s2_n6_nt, real_cd_s2_n8_nt]


cd_real_s2 = cd_real.loc[cd_real['n'] == 8]
cd_real_s2 = cd_real_s2.sort_values(by=['t'])

subjects = cd_real_s2.subjectId.unique()
predicted_cd_s2_t0 = np.zeros([np.shape(subjects)[0]])
real_cd_s2_t0 = np.zeros([np.shape(subjects)[0]])

predicted_cd_s2_t1 = np.zeros([np.shape(subjects)[0]])
real_cd_s2_t1 = np.zeros([np.shape(subjects)[0]])

predicted_cd_s2_t2 = np.zeros([np.shape(subjects)[0]])
real_cd_s2_t2 = np.zeros([np.shape(subjects)[0]])

predicted_cd_s2_t3 = np.zeros([np.shape(subjects)[0]])
real_cd_s2_t3 = np.zeros([np.shape(subjects)[0]])

predicted_cd_s2_t4 = np.zeros([np.shape(subjects)[0]])
real_cd_s2_t4 = np.zeros([np.shape(subjects)[0]])

sub = 0
for s in subjects:
    predicted_cd_s2_t0[sub] = np.mean(cd_real_s2['pred'].loc[(cd_real_s2['t'] == 0) & (cd_real_s2['subjectId'] == s)])
    real_cd_s2_t0[sub] = np.mean(cd_real_s2['behav'].loc[(cd_real_s2['t'] == 0) & (cd_real_s2['subjectId'] == s)])

    predicted_cd_s2_t1[sub] = np.mean(cd_real_s2['pred'].loc[(cd_real_s2['t'] == 1) & (cd_real_s2['subjectId'] == s)])
    real_cd_s2_t1[sub] = np.mean(cd_real_s2['behav'].loc[(cd_real_s2['t'] == 1) & (cd_real_s2['subjectId'] == s)])

    predicted_cd_s2_t2[sub] = np.mean(cd_real_s2['pred'].loc[(cd_real_s2['t'] == 2) & (cd_real_s2['subjectId'] == s)])
    real_cd_s2_t2[sub] = np.mean(cd_real_s2['behav'].loc[(cd_real_s2['t'] == 2) & (cd_real_s2['subjectId'] == s)])

    predicted_cd_s2_t3[sub] = np.mean(cd_real_s2['pred'].loc[(cd_real_s2['t'] == 3) & (cd_real_s2['subjectId'] == s)])
    real_cd_s2_t3[sub] = np.mean(cd_real_s2['behav'].loc[(cd_real_s2['t'] == 3) & (cd_real_s3['subjectId'] == s)])

    predicted_cd_s2_t4[sub] = np.mean(cd_real_s2['pred'].loc[(cd_real_s2['t'] == 4) & (cd_real_s2['subjectId'] == s)])
    real_cd_s2_t4[sub] = np.mean(cd_real_s2['behav'].loc[(cd_real_s2['t'] == 4) & (cd_real_s2['subjectId'] == s)])

    sub = sub + 1


pred_t = [predicted_cd_s2_t0, predicted_cd_s2_t1, predicted_cd_s2_t2, predicted_cd_s2_t3, predicted_cd_s2_t4]
real_t = [real_cd_s2_t0, real_cd_s2_t1, real_cd_s2_t2, real_cd_s2_t3, real_cd_s2_t4]


fig, axs = plt.subplots(1, 2)
plt.suptitle('Change detection')
axs[0].plot([3,4,6,8], np.mean(real_n, axis=1), 'ko-', label='data', markersize=4)
axs[0].fill_between([3,4,6,8], np.mean(real_n, axis=1) - np.nanstd(real_n, axis=1)/np.sqrt(90) , np.mean(real_n, axis=1) +
                    np.nanstd(real_n, axis=1)/np.sqrt(90), color='black', alpha=0.3)

axs[0].plot([3,4,6,8], np.mean(pred_n, axis=1), 'ro-', label='data', markersize=4)
axs[0].fill_between([3,4,6,8], np.mean(pred_n, axis=1) - np.nanstd(pred_n, axis=1)/np.sqrt(90) , np.mean(pred_n, axis=1) +
                    np.nanstd(pred_n, axis=1)/np.sqrt(90), color='red', alpha=0.3)

axs[0].plot([3,4,6,8], np.mean(real_n_nt, axis=1), 'ko-', label='data', markersize=4)
axs[0].fill_between([3,4,6,8], np.mean(real_n_nt, axis=1) - np.nanstd(real_n_nt, axis=1)/np.sqrt(90) , np.mean(real_n_nt, axis=1) +
                    np.nanstd(real_n_nt, axis=1)/np.sqrt(90), color='black', alpha=0.3)

axs[0].plot([3,4,6,8], np.mean(pred_n_nt, axis=1), 'ro-', label='data', markersize=4)
axs[0].fill_between([3,4,6,8], np.mean(pred_n_nt, axis=1) - np.nanstd(pred_n_nt, axis=1)/np.sqrt(90) , np.mean(pred_n_nt, axis=1) +
                    np.nanstd(pred_n_nt, axis=1)/np.sqrt(90), color='red', alpha=0.3)
axs[0].set_title("T = 1")
axs[0].set_ylabel('% correct')
axs[0].set_xlabel('Set size')
plt.ylim(0, 1.1)

axs[1].plot(list(range(0,5)), np.mean(real_t, axis=1), 'ko-', label='data', markersize=4)
axs[1].fill_between(list(range(0,5)), np.mean(real_t, axis=1) - np.nanstd(real_t, axis=1)/np.sqrt(90) , np.mean(real_t, axis=1) +
                    np.nanstd(real_t, axis=1)/np.sqrt(90), color='black', alpha=0.3)

axs[1].plot(list(range(0,5)), np.mean(pred_t, axis=1), 'ro-', label='data', markersize=4)
axs[1].fill_between(list(range(0,5)), np.mean(pred_t, axis=1) - np.nanstd(pred_t, axis=1)/np.sqrt(90) , np.mean(pred_t, axis=1) +
                    np.nanstd(pred_t, axis=1)/np.sqrt(90), color='red', alpha=0.3)
axs[1].set_title("N = 8")
axs[1].set_xlabel('Targets')
plt.ylim(0, 1.1)
plt.tight_layout()
plt.savefig('cd_ppc_all_weeks_er_final.pdf', dpi=600)
plt.close()
