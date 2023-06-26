import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# NC
nc_predicted = pd.read_csv('/Users/daniel/PycharmProjects/comp_pheno_cmdStan/results/hierarchical/tasks/nc_hierarchical/final//nc_predicted_independent_norm_priors.csv')
nc_real = pd.read_csv('/Users/daniel/PycharmProjects/comp_pheno_cmdStan/nc_data_for_stan_90s.csv')

nc_predicted = nc_predicted[nc_predicted['0'] != -2]
nc_predicted_reset = nc_predicted.reset_index()

pred = []
for t in range(0, np.shape(nc_predicted)[0]):
    tmp = nc_predicted_reset['0'][t]
    pred = np.append(pred, tmp)

nc_real['pred'] = pred
nc_real['mu_diff'] = nc_real['right_circles'] - nc_real['left_circles']

## ALL WEEKS - STD ACROSS ALL SUBJECTS
jumps = 20
nc_real = nc_real.sort_values(by=['mu_diff'])
bins = np.linspace(nc_real.mu_diff.min() - 0.01, nc_real.mu_diff.max() + 0.01, jumps)
nc_real['bin'] = pd.cut(nc_real['mu_diff'], bins)

bins_l = nc_real['bin'].unique()
subjects = nc_real.subjectId.unique()
c = 0
predicted = np.zeros([bins_l.shape[0], np.shape(subjects)[0]])
real = np.zeros([bins_l.shape[0], np.shape(subjects)[0]])
for d in bins_l:
    sub = 0
    for s in subjects:
        predicted[c, sub] = np.mean(nc_real['pred'].loc[(nc_real['bin'] == d) & (nc_real['subjectId'] == s)])
        real[c, sub] = np.mean(nc_real['key_choice'].loc[(nc_real['bin'] == d) & (nc_real['subjectId'] == s)])
        sub = sub + 1
    c = c + 1

mid = [(a.left + a.right) / 2 for a in bins_l]
mid = np.round(np.unique(mid))

ax = plt.gca()
plt.plot(mid, np.nanmean(real, axis=1), 'ko-', markersize=3, linewidth=1)
plt.fill_between(mid, np.nanmean(real, axis=1) - np.nanstd(real, axis=1)/np.sqrt(90), np.nanmean(real, axis=1) +
                    np.nanstd(real, axis=1)/np.sqrt(90), color='black', alpha=0.3)
plt.plot(mid, np.nanmean(predicted, axis=1), 'ro-', markersize=3, linewidth=1)
plt.fill_between(mid, np.nanmean(predicted, axis=1) - np.nanstd(predicted, axis=1)/np.sqrt(90), np.nanmean(predicted, axis=1) +
                    np.nanstd(predicted, axis=1)/np.sqrt(90), color='red', alpha=0.3)
plt.ylim(-0.1, 1.1)
plt.ylabel('P(right)')
plt.xlabel('circle difference (right-left)')
plt.title('Numerosity Comparison')
plt.tight_layout()
plt.savefig('nc_ppc_all_weeks_er_final.pdf', dpi=600)
plt.close()