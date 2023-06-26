# Overview

This is a repository for the tasks and data analysis used in the computational phenotype study.


## Study Structure

This study had participants perform seven tasks, broken up into three sessions, every week for 12 weeks.
Tasks were broken up such that each session took approximately the same amount of time.
Sessions were completed once a day for three consecutive days. <br>

Prior to doing the tasks for each session, subjects completed an intro questionnaire.
This questionnaire was the same for every session with the exception of the very first session,
which is slightly longer to collect demographic information.
Furthermore, during the first week only, subjects completed two extra surveys: DOSPERT and Barratt Impulsiveness.

Over the course of the weeks, each session contained the same tasks, but the order changed (i.e. Session 1 *always* included Go/No-Go and a Random Dot Kinetogram, but the order of which came first or second was switched every week).

## Tasks

All tasks were designed using jsPsych by Hanna Hillman and Jasmine Zhou.
Custom plug-ins were designed for the change detection task, numerosity comparison task.
Custom/modified css scripts were used for the slot-machine/bandit task, change detections task, and numerosity comparison task.

A list of the tasks, their abrreviations (used throughout repo), and their respective sessions are as follows:<br>
**Go/NoGo** - gng - Session 1<br>
**Random Dot Kinetogram** - rdm - Session 1<br>
**Slot Machine/ Bandit** - smb - Session 2<br>
**Numerosity Comparison** - nc - Session 2<br>
**Lottery Ticket** - lt - Session 2<br>
**Change Detection** - cd - Sesssion 3<br>
**Intertemporal Choice** - itc - Session 3<br>


## Data

Data was saved on our website server, as well as on a database managed by our website server.
To facilitate this there were 2 php scripts for each task/survey. These follow the templates provided on [jsPsych.com](https://www.jspsych.org/overview/data/#storing-data-permanently-in-a-mysql-database). <br>
Another PHP script, "save_data" allows for csvs of data to be saved
<br>
Stripped down versions can be found in /TASKS/PHP, instructions to modify are in /TASKS/README_taskCode.md

## Using this repository

The tasks in this repository will represent those used in the first week of sessions. There are some redundancies, such as the img folder, jspsych-6.0.5 folder, consent form, and intro questionnaire

**Session 1:** <br>
intro_bgr: this is the extended version that also includes demographic information<br>
bi: barratt impulsiveness survey (only collected once during first week)<br>
gng: go/no-go task<br>
pre_rdm: this was placed prior to rdm tasks to get subjects to expand their windows to full screen <br>
rdm: random dot kinetogram task<br>
consent: the consent forms subjects viewed as approved by our IRB <br>
jspsych-6.0.5: folder with the version of jsPsych used, as well as any custom plug ins, css, etc.<br>
img: folder with any images used in this session <br>

**Session 2:** <br>
intro_dsnl: this is the shortened version taken at the beginning of all other sessions<br>
dp: DOSPERT survey (only collected once during first week)<br>
smb: slot machine/bandit task <br>
nc: numerosity comparison task <br>
lt: lottery ticket task <br>
consent: the consent forms subjects viewed as approved by our IRB
jspsych-6.0.5: folder with the version of jsPsych used, as well as any custom plug ins, css, etc.<br>
img: folder with any images used in this session <br>

**Session 3:** <br>
intro_ic: this is the shortened version taken at the beginning of all other sessions<br>
itc: intertemporal choice survey task <br>
cd: change detection task <br>
consent: the consent forms subjects viewed as approved by our IRB
jspsych-6.0.5: folder with the version of jsPsych used, as well as any custom plug ins, css, etc.<br>
img: folder with any images used in this session <br>
