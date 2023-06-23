import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# SMB
smb_predicted = pd.read_csv('/Users/daniel/PycharmProjects/comp_pheno_cmdStan/results/hierarchical/tasks/smb_hierarchical/final//smb_predicted_independent_norm_priors_new.csv')
smb_real = pd.read_csv('/Users/daniel/PycharmProjects/comp_pheno_cmdStan/smb_data_for_stan_90s_updatedRegressorsNov22_rescaled_signVoverTU.csv')

smb_predicted = smb_predicted[smb_predicted['0'] != -2]
smb_predicted_reset = smb_predicted.reset_index()

pred = []
for t in range(0, np.shape(smb_predicted)[0]):
    tmp = smb_predicted_reset['0'][t]
    pred = np.append(pred, tmp)

smb_real['pred'] = pred
smb_real['ev_diff'] = smb_real.V

## ALL WEEKS - STD ACROSS ALL SUBJECTS
jumps = 15
smb_real = smb_real.sort_values(by=['ev_diff'])
bins = np.linspace(smb_real.ev_diff.min() - 0.01, smb_real.ev_diff.max() + 0.01, jumps)
smb_real['bin'] = pd.cut(smb_real['ev_diff'], bins)
subjects = smb_real.subjectId.unique()

bins_l = smb_real['bin'].unique()
c = 0
predicted_smb = np.zeros([bins_l.shape[0], np.shape(subjects)[0]])
real_smb = np.zeros([bins_l.shape[0], np.shape(subjects)[0]])
for d in bins_l:
    sub = 0
    for s in subjects:
        predicted_smb[c, sub] = np.mean(smb_real['pred'].loc[(smb_real['bin'] == d) & (smb_real['subjectId'] == s)])
        real_smb[c, sub] = np.mean(smb_real['chosen_machine'].loc[(smb_real['bin'] == d) & (smb_real['subjectId'] == s)])
        sub = sub + 1
    c = c + 1

mid = [(a.left + a.right) / 2 for a in bins_l]

ax = plt.gca()
plt.plot(mid, np.nanmean(real_smb, axis=1), 'ko-', markersize=3, linewidth=1)
plt.fill_between(mid, np.nanmean(real_smb, axis=1) - np.nanstd(real_smb, axis=1)/np.sqrt(90), np.nanmean(real_smb, axis=1) +
                    np.nanstd(real_smb, axis=1)/np.sqrt(90), color='black', alpha=0.3)
plt.plot(mid, np.nanmean(predicted_smb, axis=1), 'ro-', markersize=3, linewidth=1)
plt.fill_between(mid, np.nanmean(predicted_smb, axis=1) - np.nanstd(predicted_smb, axis=1)/np.sqrt(90), np.nanmean(predicted_smb, axis=1) +
                    np.nanstd(predicted_smb, axis=1)/np.sqrt(90), color='red', alpha=0.3)
plt.ylim(-0.1, 1.1)
plt.ylabel('Choice Probability')
plt.xlabel('EV difference')
plt.title('Two-armed bandit')
plt.tight_layout()
plt.savefig('smb_ppc_all_weeks_er_final.pdf', dpi=600)
plt.close()
