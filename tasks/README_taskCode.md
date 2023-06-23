# About the task code:
The tasks in this folder were coded by Hanna Hillman with the assistance of Jasmine Zhou, in Prof. Samuel J Gershman's Computational Cognitive Neuroscience Lab, Harvard University. <br>

For specific questions about the task code please reach out to Hanna Hillman at hhillman231@gmail.com
For specific questions aobut the parameters chosen for each task, please review the literature document.

## jsPsych
All tasks were coded using jsPsych, version 6.0.5 <br>

Citation: <br>
de Leeuw, J. R. (2015). jsPsych: A JavaScript library for creating behavioral experiments in a web browser. Behavior Research Methods, 47(1), 1-12. doi:10.3758/s13428-014-0458-y.


## Task Construction

**All tasks have the following sections**
- Quick change variables: a list of parameters that can be easily changed, such as the number of trials per block, number of blocks, etc.
- Fullscreen: puts their browser window into fullscree (RDM does not have this but has pre_rdm which asks them to do it manually, see pre_rdm task for clarification)
- Instructions: gives them task instructions
- Instruction comprehension questions: 2 multiple choice questions that are meant to ensure that they understand the instructions. If one or more of the subjects answers are incorrect, it will take them back to the beginning of the instructions. This can occur a maximum of three times (looping back to instructions) before it will proceed with practice block or experimental block (if there is no practice block, see experiment features below)
- Post-task feedback: series of questions after each task, there are two versions, one for subjects who have triggered a suspicion flag in any of their experimental blocks and one for participants who have not (in itc only one version, for the same reason it's bonus is a flat rate*) both versions ask the same questions:
    - Did you experience any technological issues during this task? (i.e. experiment glitch, accidentally closed window, mixed up keys, etc.)
    - Did you experience any major distractions or interruptions during this task? (If yes they can write out an explanation)
    - On the following scale from 1 to 4, how engaging (attention holding) did you find this task? (1=Not at all engaging, 2=Somewhat not engaging, 3=Somewhat engaging, 4= Very engaging)
    - Is there anything else you would like us to know? (Can be a specific concern or general suggestion for improvement)
- CSV Save: saves the data as a csv on the server given the appropriate PHP scripts (see PHP scripts section)
- Database save: saves the data on a SQL database given the appropriate PHP scripts (see PHP scripts section)
- save_wait variable: ensures that there is enough time to save, without this you may not get the data
- "Dashboard": Allows subjects to take a brief pause before continuing to the next task (the last tasks of each section have this commented out )

**Experiment Features**
|                       |Practice Block|Exp Blocks| Suspicion Flags                        | Bonus Pts |
|-----------------------|--------------|----------|----------------------------------------|-----------|
|Go/NoGo                | yes          | 3        | all_one <br> time_outs <br> error_rate | 32        |
|Random Dot Kinetogram  | yes          | 4        | all_one <br> time_outs <br> error_rate | 20        |
|Slot Machine/ Bandit   | yes          | 30       | all_one <br> error_rate                | 20        |
|Numerosity Comparison  | yes          | 2        | all_one <br> time_outs <br> error_rate | 20        |
|Lottery Ticket         | no           | 3        | all_one <br> one_off                   | 20        |
|Change Detection       | yes          | 5        | all_one <br> error_block               | 20        |
|Intertemporal Choice   | no           | 3        | all_one                                | 10*       |

/* 10 points is a flat rate, no bonus points are calculated based on performance due to the nature of the task being a subjective preference

## Suspicion Flags
Suspicion flags have their own column in the data. They vary based on task (see Experiment Features table above). <br>


**Types of Suspicion flags:**
- "all_one" = subject chose the same answer for every choice in that block*
- "time_outs" = subject did not make any active choice for any of the trials in that block*
- "error_rate" or "error_block" = subject has a suspiciously high error rate - the parameter for this is error_rate/error_cap and is assigned in quick change variables at the top of task code
- "one_off" = this is only for the lottery ticket task, it triggers when subject selects the lower amount of money in a lower/higher amount choice when the odds are 100% that they could get the higher amount


/* in the Go/NoGo task, even though not pressing space is technically a choice option, the "all_one" flag will only trigger if they press the space bar for every square presented, if they press nothing throughout the block the "time_outs" flag will be triggered


**flags during practice block**
- If a flag is triggered during the practice block, subjects will see a notification that 'perhaps they didn't understand the task completely' and then they are sent back to the instructions. Once they have repeated the instructions (they do not have to repeat the instruction comprehension questions) they will have another practice block. If a flag is triggered again, the process will repeat. This can happen up to three times, at which point they just proceed to the experimental block.
- This repeating practice block takes up a lot of space in the task codes, starting with "no_sus" and ending with "sus_procedure" variables.


## End of task code - init
At the end of the task code you will see a function that looks like the following:
>     jsPsych.init({
>      timeline: timeline,
>      on_finish: function() {window.location.href = "pre_rdm.html" + '?workerId=' + subjectId + '&hitId=' + hitId + '&assignmentId=' + assignmentId}
>      // on_finish: function() {     // for testing purposes
>      //   jsPsych.data.displayData();
>      //   jsPsych.data.get().localSave('csv','gng_0_CSV.csv');
>      // },
>    });

These lines:
>     jsPsych.init({
>      timeline: timeline,

Are basically what runs your entire task.
<br>
<br>
This line:
>      on_finish: function() {window.location.href = "pre_rdm.html" + '?workerId=' + subjectId + '&hitId=' + hitId + '&assignmentId=' + assignmentId}

Takes you to the next task script (in this case pre_rdm.html), and carries over the most important identifying information that can be used in the
that task (like subjectId).
<br>
<br>
If you were to comment out the above line, and activate these lines (which are usually commented out):
>      // on_finish: function() {     
>      //   jsPsych.data.displayData();
>      //   jsPsych.data.get().localSave('csv','gng_0_CSV.csv');
>      // },

The data from that task will be shown on the screen after you finish with it, as well as save to a CSV on your local machine.



## In the folder
To make the task code run, you need the following things in the folder with the task code html files:

**1. jsPsych-6.0.5**


**2. img folder**
This includes all of the images you need for the scripts, the img folder included in session1, session2, and sesssion3 folders on this repo are the identical


**3. consent form**
Only if you're running the intro quiz


**4. PHP Scripts**
For each HTML script run, there are two PHP scripts to save the data into the database, and a third that help save to csv

1. write_data_[task ID].php
  - this one is long, but you only need to change the top line to point to the other php script
2. database_config_[task ID].php
  - this one is short, just fill it in
3. save_data.php
  - for csv save, no changes need to be made, just include in the folder with tasks


These have been taken from jsPsych templates. Stripped versions of them can be found in the PHP folder
