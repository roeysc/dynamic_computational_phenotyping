import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# ITC
itc_predicted = pd.read_csv('/Users/daniel/PycharmProjects/comp_pheno_cmdStan/results/hierarchical/tasks/itc_hierarchical/final//itc_predicted_independent_norm_priors.csv')
itc_real = pd.read_csv('/Users/daniel/PycharmProjects/comp_pheno_cmdStan/itc_data_for_stan_90s.csv')

itc_predicted = itc_predicted[itc_predicted['0'] != -2]
itc_predicted_reset = itc_predicted.reset_index()

pred = []
for t in range(0, np.shape(itc_predicted)[0]):
    tmp = itc_predicted_reset['0'][t]
    pred = np.append(pred, tmp)

itc_real['pred'] = pred

## ALL WEEKS - STD ACROSS ALL SUBJECTS
bins = 5
# short delay
itc_real_short = itc_real.loc[(itc_real['later_delay'] >= 7) & (itc_real['later_delay'] < 100)]
itc_real_short['small/large'] = np.array(itc_real_short.small_amount/itc_real_short.large_amount)
itc_real_short = itc_real_short.sort_values(by=['small/large'])

# medium delay
itc_real_medium = itc_real.loc[(itc_real['later_delay'] >= 100) & (itc_real['later_delay'] < 171)]
itc_real_medium['small/large'] = itc_real_medium.small_amount/itc_real_medium.large_amount
itc_real_medium = itc_real_medium.sort_values(by=['small/large'])

# long delay
itc_real_long = itc_real.loc[(itc_real['later_delay'] >= 171) & (itc_real['later_delay'] <= 260)]
itc_real_long['small/large'] = itc_real_long.small_amount/itc_real_long.large_amount
itc_real_long = itc_real_long.sort_values(by=['small/large'])

bins_s = np.linspace(itc_real_short['small/large'].min() - 0.01, itc_real_short['small/large'].max() + 0.01, bins)
itc_real_short['bin'] = pd.cut(itc_real_short['small/large'], bins_s)

bins_m = np.linspace(itc_real_medium['small/large'].min() - 0.01, itc_real_medium['small/large'].max() + 0.01, bins)
itc_real_medium['bin'] = pd.cut(itc_real_medium['small/large'], bins_m)

bins_l = np.linspace(itc_real_long['small/large'].min() - 0.01, itc_real_long['small/large'].max() + 0.01, bins)
itc_real_long['bin'] = pd.cut(itc_real_long['small/large'], bins_l)

subjects = itc_real_short.subjectId.unique()
bins_ls = itc_real_short['bin'].unique()
c = 0
predicted_itc_s = np.zeros([bins_ls.shape[0], np.shape(subjects)[0]])
real_itc_s = np.zeros([bins_ls.shape[0], np.shape(subjects)[0]])
for d in bins_ls:
    sub = 0
    for s in subjects:
        predicted_itc_s[c, sub] = np.mean(itc_real_short['pred'].loc[(itc_real_short['bin'] == d) & (itc_real_short['subjectId'] == s)])
        real_itc_s[c, sub] = np.mean(itc_real_short['choice'].loc[(itc_real_short['bin'] == d) & (itc_real_short['subjectId'] == s)])
        sub = sub + 1
    c = c + 1

mid_s = [(a.left + a.right) / 2 for a in bins_ls]

bins_lm = itc_real_medium['bin'].unique()
c = 0
predicted_itc_m = np.zeros([bins_lm.shape[0], np.shape(subjects)[0]])
real_itc_m = np.zeros([bins_lm.shape[0], np.shape(subjects)[0]])
for d in bins_lm:
    sub = 0
    for s in subjects:
        predicted_itc_m[c, sub] = np.mean(itc_real_medium['pred'].loc[(itc_real_medium['bin'] == d) & (itc_real_medium['subjectId'] == s)])
        real_itc_m[c, sub] = np.mean(itc_real_medium['choice'].loc[(itc_real_medium['bin'] == d) & (itc_real_medium['subjectId'] == s)])
        sub = sub + 1
    c = c + 1

mid_m = [(a.left + a.right) / 2 for a in bins_lm]

bins_ll = itc_real_long['bin'].unique()
c = 0
predicted_itc_l = np.zeros([bins_ll.shape[0], np.shape(subjects)[0]])
real_itc_l = np.zeros([bins_ll.shape[0], np.shape(subjects)[0]])
for d in bins_ll:
    sub = 0
    for s in subjects:
        predicted_itc_l[c, sub] = np.mean(itc_real_long['pred'].loc[(itc_real_long['bin'] == d) & (itc_real_long['subjectId'] == s)])
        real_itc_l[c, sub] = np.mean(itc_real_long['choice'].loc[(itc_real_long['bin'] == d) & (itc_real_long['subjectId'] == s)])
        sub = sub + 1
    c = c + 1

mid_l = [(a.left + a.right) / 2 for a in bins_ll]

fig, axs = plt.subplots(1,3)
plt.suptitle('Intertemporal Choice')
axs[0].plot(mid_s, np.nanmean(real_itc_s, axis=1), 'ko-', label='data', markersize=3, linewidth=1)
axs[0].fill_between(mid_s, np.nanmean(real_itc_s, axis=1) - np.nanstd(real_itc_s, axis=1)/np.sqrt(90), np.nanmean(real_itc_s, axis=1) +
                    np.nanstd(real_itc_s, axis=1)/np.sqrt(90), color='black', alpha=0.3)

axs[0].plot(mid_s, np.nanmean(predicted_itc_s, axis=1), 'ro-', label='model', linewidth=1, markersize=3)
axs[0].fill_between(mid_s, np.nanmean(predicted_itc_s, axis=1) - np.nanstd(predicted_itc_s, axis=1)/np.sqrt(90), np.nanmean(predicted_itc_s, axis=1) +
                    np.nanstd(predicted_itc_s, axis=1)/np.sqrt(90), color='red', alpha=0.3)
axs[0].set_title("Short Delay - 7 to 99 days")
#axs[0].set_ylabel('Probability to choose "later"')
axs[0].set_xlabel('Small/large ratio')
axs[0].set_ylim(-0.1, 1.1)
axs[0].set_xlim(0, 1.2)
plt.tight_layout()

axs[1].plot(mid_m, np.nanmean(real_itc_m, axis=1), 'ko-', label='data', markersize=3, linewidth=1)
axs[1].fill_between(mid_m, np.nanmean(real_itc_m, axis=1) - np.nanstd(real_itc_m, axis=1)/np.sqrt(90), np.nanmean(real_itc_m, axis=1) +
                    np.nanstd(real_itc_m, axis=1)/np.sqrt(90), color='black', alpha=0.3)

axs[1].plot(mid_m, np.nanmean(predicted_itc_m, axis=1), 'ro-', label='model', linewidth=1, markersize=3)
axs[1].fill_between(mid_m, np.nanmean(predicted_itc_m, axis=1) - np.nanstd(predicted_itc_m, axis=1)/np.sqrt(90), np.nanmean(predicted_itc_m, axis=1) +
                    np.nanstd(predicted_itc_m, axis=1)/np.sqrt(90), color='red', alpha=0.3)
axs[1].set_title("Medium Delay - 100 to 170 days")
axs[1].set_xlabel('Small/large ratio')
axs[1].set_ylim(-0.1, 1.1)
axs[1].set_xlim(0, 1.2)
plt.tight_layout()

axs[2].plot(mid_l, np.nanmean(real_itc_l, axis=1), 'ko-', label='data', markersize=3, linewidth=1)
axs[2].fill_between(mid_l, np.nanmean(real_itc_l, axis=1) - np.nanstd(real_itc_l, axis=1)/np.sqrt(90), np.nanmean(real_itc_l, axis=1) +
                    np.nanstd(real_itc_l, axis=1)/np.sqrt(90), color='black', alpha=0.3)

axs[2].plot(mid_l, np.nanmean(predicted_itc_l, axis=1), 'ro-', label='model', linewidth=1, markersize=3)
axs[2].fill_between(mid_l, np.nanmean(predicted_itc_l, axis=1) - np.nanstd(predicted_itc_l, axis=1)/np.sqrt(90), np.nanmean(predicted_itc_l, axis=1) +
                    np.nanstd(predicted_itc_l, axis=1)/np.sqrt(90), color='red', alpha=0.3)
axs[2].set_title("Long Delay - 171 to 260 days ")
axs[2].set_xlabel('Small/large ratio')

axs[2].set_ylim(-0.1, 1.1)
axs[2].set_xlim(0, 1.2)
plt.suptitle('Intertemporal choice')
plt.tight_layout()
plt.savefig('itc_ppc_all_weeks_er_final.pdf', dpi=600)
plt.close()