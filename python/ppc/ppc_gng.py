import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# GNG
gng_predicted = pd.read_csv('/Users/daniel/PycharmProjects/comp_pheno_cmdStan/results/hierarchical/tasks/gng_hierarchical/final//gng_predicted_independent_norm_priors.csv')
gng_real = pd.read_csv('/Users/daniel/PycharmProjects/comp_pheno_cmdStan/gng_data_for_stan_90s.csv')


all_subjects = gng_real.subjectId
subjects = all_subjects.unique()
# exclude_gng_subjects_with_accuracy_below_055:
sub_exclude = [1, 2, 3, 6, 8, 9, 15, 16, 17, 28, 33, 34, 37, 39, 44, 59, 60, 68, 70, 71, 74, 76, 78, 86]

for idx in range(0, np.shape(sub_exclude)[0]):
    gng_real = gng_real.drop(
        gng_real[(gng_real['subjectId'] == subjects[sub_exclude[idx] - 1])].index)

gng_predicted = gng_predicted[gng_predicted['0'] != -2]
gng_predicted_reset = gng_predicted.reset_index()


# way 2
pred = []
for t in range(0, np.shape(gng_predicted)[0]):
    tmp = gng_predicted_reset['0'][t]
    pred = np.append(pred, tmp)

gng_real['pred'] = pred

## ALL WEEKS - STD ACROSS ALL SUBJECTS
data1_real_s = np.zeros((gng_real.subjectId.unique().shape[0], 20))
data2_real_s = np.zeros((gng_real.subjectId.unique().shape[0], 20))
data3_real_s = np.zeros((gng_real.subjectId.unique().shape[0], 20))
data4_real_s = np.zeros((gng_real.subjectId.unique().shape[0], 20))

data1_pred_s = np.zeros((gng_real.subjectId.unique().shape[0], 20))
data2_pred_s = np.zeros((gng_real.subjectId.unique().shape[0], 20))
data3_pred_s = np.zeros((gng_real.subjectId.unique().shape[0], 20))
data4_pred_s = np.zeros((gng_real.subjectId.unique().shape[0], 20))

c=0
for s in gng_real.subjectId.unique():

    data1_real = np.zeros((gng_real.weekId[gng_real.subjectId == s].unique().shape[0], 20))
    data2_real = np.zeros((gng_real.weekId[gng_real.subjectId == s].unique().shape[0], 20))
    data3_real = np.zeros((gng_real.weekId[gng_real.subjectId == s].unique().shape[0], 20))
    data4_real = np.zeros((gng_real.weekId[gng_real.subjectId == s].unique().shape[0], 20))

    data1_pred = np.zeros((gng_real.weekId[gng_real.subjectId == s].unique().shape[0], 20))
    data2_pred = np.zeros((gng_real.weekId[gng_real.subjectId == s].unique().shape[0], 20))
    data3_pred = np.zeros((gng_real.weekId[gng_real.subjectId == s].unique().shape[0], 20))
    data4_pred = np.zeros((gng_real.weekId[gng_real.subjectId == s].unique().shape[0], 20))

    c1 = 0

    for w in gng_real.weekId[gng_real.subjectId == s].unique():
        data1b1_tmp = gng_real[(gng_real.subjectId == s) & (gng_real.block == 'block_1') & (gng_real.cond == 1)
        & (gng_real.weekId == w)]
        data1b2_tmp = gng_real[(gng_real.subjectId == s) & (gng_real.block == 'block_2') & (gng_real.cond == 1)
        & (gng_real.weekId == w)]
        data1b3_tmp = gng_real[(gng_real.subjectId == s) & (gng_real.block == 'block_3') & (gng_real.cond == 1)
        & (gng_real.weekId == w)]

        data2b1_tmp = gng_real[(gng_real.subjectId == s) & (gng_real.block == 'block_1') & (gng_real.cond == 2)
        & (gng_real.weekId == w)]
        data2b2_tmp = gng_real[(gng_real.subjectId == s) & (gng_real.block == 'block_2') & (gng_real.cond == 2)
        & (gng_real.weekId == w)]
        data2b3_tmp = gng_real[(gng_real.subjectId == s) & (gng_real.block == 'block_3') & (gng_real.cond == 2)
        & (gng_real.weekId == w)]

        data3b1_tmp = gng_real[(gng_real.subjectId == s) & (gng_real.block == 'block_1') & (gng_real.cond == 3)
        & (gng_real.weekId == w)]
        data3b2_tmp = gng_real[(gng_real.subjectId == s) & (gng_real.block == 'block_2') & (gng_real.cond == 3)
        & (gng_real.weekId == w)]
        data3b3_tmp = gng_real[(gng_real.subjectId == s) & (gng_real.block == 'block_3') & (gng_real.cond == 3)
        & (gng_real.weekId == w)]

        data4b1_tmp = gng_real[(gng_real.subjectId == s) & (gng_real.block == 'block_1') & (gng_real.cond == 4)
        & (gng_real.weekId == w)]
        data4b2_tmp = gng_real[(gng_real.subjectId == s) & (gng_real.block == 'block_2') & (gng_real.cond == 4)
        & (gng_real.weekId == w)]
        data4b3_tmp = gng_real[(gng_real.subjectId == s) & (gng_real.block == 'block_3') & (gng_real.cond == 4)
        & (gng_real.weekId == w)]

        data1_real[c1, :] = (np.array(data1b1_tmp.choice) + np.array(data1b2_tmp.choice) + np.array(data1b3_tmp.choice))/3
        data2_real[c1, :] = (np.array(data2b1_tmp.choice) + np.array(data2b2_tmp.choice) + np.array(data2b3_tmp.choice))/3
        data3_real[c1, :] = (np.array(data3b1_tmp.choice) + np.array(data3b2_tmp.choice) + np.array(data3b3_tmp.choice))/3
        data4_real[c1, :] = (np.array(data4b1_tmp.choice) + np.array(data4b2_tmp.choice) + np.array(data4b3_tmp.choice))/3

        data1_pred[c1, :] = (np.array(data1b1_tmp.pred) + np.array(data1b2_tmp.pred) + np.array(data1b3_tmp.pred))/3
        data2_pred[c1, :] = (np.array(data2b1_tmp.pred) + np.array(data2b2_tmp.pred) + np.array(data2b3_tmp.pred))/3
        data3_pred[c1, :] = (np.array(data3b1_tmp.pred) + np.array(data3b2_tmp.pred) + np.array(data3b3_tmp.pred))/3
        data4_pred[c1, :] = (np.array(data4b1_tmp.pred) + np.array(data4b2_tmp.pred) + np.array(data4b3_tmp.pred))/3

        c1 = c1 + 1

    data1_real_s[c, :] = np.mean(data1_real, axis=0)
    data2_real_s[c, :] = np.mean(data2_real, axis=0)
    data3_real_s[c, :] = np.mean(data3_real, axis=0)
    data4_real_s[c, :] = np.mean(data4_real, axis=0)

    data1_pred_s[c, :] = np.mean(data1_pred, axis=0)
    data2_pred_s[c, :] = np.mean(data2_pred, axis=0)
    data3_pred_s[c, :] = np.mean(data3_pred, axis=0)
    data4_pred_s[c, :] = np.mean(data4_pred, axis=0)

    c = c + 1


data1_real_final_mean = np.nanmean(data1_real_s, axis=0)
data2_real_final_mean = np.nanmean(data2_real_s, axis=0)
data3_real_final_mean = np.nanmean(data3_real_s, axis=0)
data4_real_final_mean = np.nanmean(data4_real_s, axis=0)

data1_real_final_std = np.nanstd(data1_real_s, axis=0)
data2_real_final_std = np.nanstd(data2_real_s, axis=0)
data3_real_final_std = np.nanstd(data3_real_s, axis=0)
data4_real_final_std = np.nanstd(data4_real_s, axis=0)

data1_pred_final_mean = np.nanmean(data1_pred_s, axis=0)
data2_pred_final_mean = np.nanmean(data2_pred_s, axis=0)
data3_pred_final_mean = np.nanmean(data3_pred_s, axis=0)
data4_pred_final_mean = np.nanmean(data4_pred_s, axis=0)

data1_pred_final_std = np.nanstd(data1_pred_s, axis=0)
data2_pred_final_std = np.nanstd(data2_pred_s, axis=0)
data3_pred_final_std = np.nanstd(data3_pred_s, axis=0)
data4_pred_final_std = np.nanstd(data4_pred_s, axis=0)

fig, ((ax1, ax3), (ax5, ax7)) = plt.subplots(2, 2)
ax2 = ax1.twiny()
ax1.plot(data1_pred_final_mean, 'ro-', linewidth=1, markersize=3)
ax1.fill_between(np.arange(0,20),data1_pred_final_mean - data1_pred_final_std/np.sqrt(66), data1_pred_final_mean + data1_pred_final_std/np.sqrt(66), color='red', alpha=0.3)
ax1.set_xlabel('Trial')
ax1.set_xlim(-1, 21)
ax2.plot(data1_real_final_mean, 'ko-', markersize=3, linewidth=1)
ax2.fill_between(np.arange(0,20),data1_real_final_mean - data1_real_final_std/np.sqrt(66), data1_real_final_mean + data1_real_final_std/np.sqrt(66), color='black', alpha=0.3)
ax2.set_xticks([])
ax2.set_xlim(-1, 21)
plt.ylim(0.35, 1.1)
ax1.set_ylabel('P(Go)')
ax1.set_title('Go to win')
plt.tight_layout()

ax4 = ax3.twiny()
ax3.plot(data2_pred_final_mean, 'ro-', linewidth=1, markersize=3)
ax3.fill_between(np.arange(0,20),data2_pred_final_mean - data2_pred_final_std/np.sqrt(66), data2_pred_final_mean + data2_pred_final_std/np.sqrt(66), color='red', alpha=0.3)
ax3.set_xlabel('Trial')
ax3.set_xlim(-1, 21)
ax4.plot(data2_real_final_mean, 'ko-', linewidth=1, markersize=3)
ax4.fill_between(np.arange(0,20),data2_real_final_mean - data2_real_final_std/np.sqrt(66), data2_real_final_mean + data2_real_final_std/np.sqrt(66), color='black', alpha=0.3)
ax4.set_xticks([])
ax4.set_xlim(-1, 21)
plt.ylim(-0.1, 0.65)
ax3.set_ylabel('P(Go)')
ax3.set_title('NoGo to win')
plt.tight_layout()

ax6 = ax5.twiny()
ax5.plot(data3_pred_final_mean, 'ro-', linewidth=1, markersize=3)
ax5.fill_between(np.arange(0,20),data3_pred_final_mean - data3_pred_final_std/np.sqrt(66), data3_pred_final_mean + data3_pred_final_std/np.sqrt(66), color='red', alpha=0.3)
ax5.set_xlabel('Trial')
ax5.set_xlim(-1, 21)
ax6.plot(data3_real_final_mean, 'ko-', linewidth=1, markersize=3)
ax6.fill_between(np.arange(0,20),data3_real_final_mean - data3_real_final_std/np.sqrt(66), data3_real_final_mean + data3_real_final_std/np.sqrt(66), color='black', alpha=0.3)
ax6.set_xticks([])
ax6.set_xlim(-1, 21)
plt.ylim(0.35, 1.1)
ax5.set_ylabel('P(Go)')
ax5.set_title('Go to avoid')
plt.tight_layout()

ax8 = ax7.twiny()
ax7.plot(data4_pred_final_mean, 'ro-', linewidth=1, markersize=3)
ax7.fill_between(np.arange(0,20),data4_pred_final_mean - data4_pred_final_std/np.sqrt(66), data4_pred_final_mean + data4_pred_final_std/np.sqrt(66), color='red', alpha=0.3)
ax7.set_xlabel('Trial')
ax7.set_xlim(-1, 21)
ax8.plot(data4_real_final_mean, 'ko-', linewidth=1, markersize=3)
ax8.fill_between(np.arange(0,20),data4_real_final_mean - data4_real_final_std/np.sqrt(66), data4_real_final_mean + data4_real_final_std/np.sqrt(66), color='black', alpha=0.3)
ax8.set_xticks([])
ax8.set_xlim(-1, 21)
plt.ylim(-0.1, 0.65)
ax7.set_ylabel('P(Go)')
ax7.set_title('NoGo to avoid')
plt.tight_layout()

plt.savefig('gng_ppc_all_weeks_er_final_old_model.pdf', dpi=600)
plt.close()
