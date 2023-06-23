import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# LT
lt_predicted = pd.read_csv('/Users/daniel/PycharmProjects/comp_pheno_cmdStan/results/hierarchical/tasks/lt_hierarchical/final//lt_predicted_independent_norm_priors_remove_80_per.csv')
lt_real = pd.read_csv('/Users/daniel/PycharmProjects/comp_pheno_cmdStan/lt_data_for_stan_90s.csv')

subjects = lt_real.subjectId.unique()

lt_sub2exclude = [1, 4, 7, 8, 9, 11, 14, 15, 19, 23, 24, 25, 28, 32, 35, 36, 37, 47, 50, 53, 58, 60, 61, 65, 66, 68, 70, 74, 77, 83, 84, 88]

for s in lt_sub2exclude:
    lt_real = lt_real.drop(lt_real[(lt_real['subjectId'] == subjects[s-1])].index)


lt_real = lt_real[lt_real.top_p != 100]

lt_predicted = lt_predicted[lt_predicted['0'] != -2]
lt_predicted_reset = lt_predicted.reset_index()

pred = []
for t in range(0, np.shape(lt_predicted)[0]):
    tmp = lt_predicted_reset['0'][t]
    pred = np.append(pred, tmp)


## ALL WEEKS - STD ACROSS ALL SUBJECTS
bins = 15

lt_real['pred'] = pred
lt_real['ev1'] = np.array(lt_real.top_p*lt_real.hi_narr/100 + lt_real.low_p*lt_real.lo_narr/100)
lt_real['ev2'] = np.array(lt_real.top_p*lt_real.hi_wide/100 + lt_real.low_p*lt_real.lo_wide/100)
lt_real['ev_diff'] = lt_real['ev1'] - lt_real['ev2']

lt_real_bl1 = lt_real.loc[lt_real['block'] == 'block_1']
lt_real_bl2 = lt_real.loc[lt_real['block'] == 'block_2']
lt_real_bl3 = lt_real.loc[lt_real['block'] == 'block_3']

lt_real_bl1 = lt_real_bl1.sort_values(by=['ev_diff'])
lt_real_bl2 = lt_real_bl2.sort_values(by=['ev_diff'])
lt_real_bl3 = lt_real_bl3.sort_values(by=['ev_diff'])

bins1 = np.linspace(lt_real_bl1.ev_diff.min() - 0.01, lt_real_bl1.ev_diff.max() + 0.01, bins)
lt_real_bl1['bin'] = pd.cut(lt_real_bl1['ev_diff'], bins1)

bins2 = np.linspace(lt_real_bl2.ev_diff.min() - 0.01, lt_real_bl2.ev_diff.max() + 0.01, bins)
lt_real_bl2['bin'] = pd.cut(lt_real_bl2['ev_diff'], bins2)

bins3 = np.linspace(lt_real_bl3.ev_diff.min() - 0.01, lt_real_bl3.ev_diff.max() + 0.01, bins)
lt_real_bl3['bin'] = pd.cut(lt_real_bl3['ev_diff'], bins3)

bins_l1 = lt_real_bl1['bin'].unique()
subjects = lt_real_bl1.subjectId.unique()
c = 0
predicted_lt1 = np.zeros([bins_l1.shape[0], np.shape(subjects)[0]])
real_lt1 = np.zeros([bins_l1.shape[0], np.shape(subjects)[0]])
for d in bins_l1:
    sub = 0
    for s in subjects:
        predicted_lt1[c, sub] = np.mean(lt_real_bl1['pred'].loc[(lt_real_bl1['bin'] == d) & (lt_real_bl1['subjectId'] == s)])
        real_lt1[c, sub] = np.mean(lt_real_bl1['choice'].loc[(lt_real_bl1['bin'] == d) & (lt_real_bl1['subjectId'] == s)])
        sub = sub + 1
    c = c + 1


bins_l2 = lt_real_bl2['bin'].unique()
subjects = lt_real_bl2.subjectId.unique()
c = 0
predicted_lt2 = np.zeros([bins_l2.shape[0], np.shape(subjects)[0]])
real_lt2 = np.zeros([bins_l2.shape[0], np.shape(subjects)[0]])
for d in bins_l2:
    sub = 0
    for s in subjects:
        predicted_lt2[c, sub] = np.mean(lt_real_bl2['pred'].loc[(lt_real_bl2['bin'] == d) & (lt_real_bl2['subjectId'] == s)])
        real_lt2[c, sub] = np.mean(lt_real_bl2['choice'].loc[(lt_real_bl2['bin'] == d) & (lt_real_bl2['subjectId'] == s)])
        sub = sub + 1
    c = c + 1


bins_l3 = lt_real_bl3['bin'].unique()
subjects = lt_real_bl3.subjectId.unique()
c = 0
predicted_lt3 = np.zeros([bins_l3.shape[0], np.shape(subjects)[0]])
real_lt3 = np.zeros([bins_l3.shape[0], np.shape(subjects)[0]])
for d in bins_l3:
    sub = 0
    for s in subjects:
        predicted_lt3[c, sub] = np.mean(lt_real_bl3['pred'].loc[(lt_real_bl3['bin'] == d) & (lt_real_bl3['subjectId'] == s)])
        real_lt3[c, sub] = np.mean(lt_real_bl3['choice'].loc[(lt_real_bl3['bin'] == d) & (lt_real_bl3['subjectId'] == s)])
        sub = sub + 1
    c = c + 1

mid1 = [(a.left + a.right) / 2 for a in bins_l1]
mid2 = [(a.left + a.right) / 2 for a in bins_l2]
mid3 = [(a.left + a.right) / 2 for a in bins_l3]

fig, axs = plt.subplots(1,3)
plt.suptitle('Lottery Ticket')
axs[0].plot(mid1, np.nanmean(real_lt1, axis=1), 'ko-', label='data', markersize=3, linewidth=1)
axs[0].fill_between(mid1, np.nanmean(real_lt1, axis=1) - np.nanstd(real_lt1, axis=1)/np.sqrt(90-32), np.nanmean(real_lt1, axis=1) +
                    np.nanstd(real_lt1, axis=1)/np.sqrt(90-32), color='black', alpha=0.3)

axs[0].plot(mid1, np.nanmean(predicted_lt1, axis=1), 'ro-', label='model', linewidth=1, markersize=3)
axs[0].fill_between(mid1, np.nanmean(predicted_lt1, axis=1) - np.nanstd(predicted_lt1, axis=1)/np.sqrt(90-32), np.nanmean(predicted_lt1, axis=1) +
                    np.nanstd(predicted_lt1, axis=1)/np.sqrt(90-32), color='red', alpha=0.3)
axs[0].set_title("Block 1")
axs[0].set_ylabel('Choice Probability')
axs[0].set_xlabel('EV difference')
axs[0].set_xlim([-2.5, 2.5])

axs[1].plot(mid2, np.nanmean(real_lt2, axis=1), 'ko-', label='data', markersize=3, linewidth=1)
axs[1].fill_between(mid2, np.nanmean(real_lt2, axis=1) - np.nanstd(real_lt2, axis=1)/np.sqrt(90-32), np.nanmean(real_lt2, axis=1) +
                    np.nanstd(real_lt2, axis=1)/np.sqrt(90-32), color='black', alpha=0.3)

axs[1].plot(mid2, np.nanmean(predicted_lt2, axis=1), 'ro-', label='model', linewidth=1, markersize=3)
axs[1].fill_between(mid2, np.nanmean(predicted_lt2, axis=1) - np.nanstd(predicted_lt2, axis=1)/np.sqrt(90-32), np.nanmean(predicted_lt2, axis=1) +
                    np.nanstd(predicted_lt2, axis=1)/np.sqrt(90-32), color='red', alpha=0.3)
axs[1].set_title("Block 2")
axs[1].set_xlabel('EV difference')
axs[1].set_xlim([-70, 70])


axs[2].plot(mid3, np.nanmean(real_lt3, axis=1), 'ko-', label='data', markersize=3, linewidth=1)
axs[2].fill_between(mid3, np.nanmean(real_lt3, axis=1) - np.nanstd(real_lt3, axis=1)/np.sqrt(90-32), np.nanmean(real_lt3, axis=1) +
                    np.nanstd(real_lt3, axis=1)/np.sqrt(90-32), color='black', alpha=0.3)

axs[2].plot(mid3, np.nanmean(predicted_lt3, axis=1), 'ro-', label='model', linewidth=1, markersize=3)
axs[2].fill_between(mid3, np.nanmean(predicted_lt3, axis=1) - np.nanstd(predicted_lt3, axis=1)/np.sqrt(90-32), np.nanmean(predicted_lt3, axis=1) +
                    np.nanstd(predicted_lt3, axis=1)/np.sqrt(90-32), color='red', alpha=0.3)

axs[2].set_xlim([-210, 210])

axs[2].set_title("Block 3")
axs[2].set_xlabel('EV difference')
plt.tight_layout()
plt.savefig('lottery_ppc_all_weeks_er_final_excl_80_per.pdf', dpi=600)
plt.close()