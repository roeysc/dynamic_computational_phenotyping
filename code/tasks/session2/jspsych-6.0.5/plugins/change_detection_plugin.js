jsPsych.plugins["change-detection-task"] = (function() {
	var plugin = {};
	plugin.info = {
	    name: "change-detection-task",
	    parameters: {
        type: {
		      task_type: jsPsych.plugins.parameterType.STRING,
		      pretty_name: "Task type",
		      default: undefined,
		      description: '"Same" or "Different" depending on whether the second array is same or different from the first'
		    },
				stim_duration: {
		      type: jsPsych.plugins.parameterType.INT,
		      pretty_name: "Trial duration",
		      default: 1000,
		      description: "The length of time that the stimulus (circular array of squares) is presented for"
		    },
        beginning_fixation_dur: {
		      type: jsPsych.plugins.parameterType.INT,
		      pretty_name: "Fixation duration before the stimulus is shown at all",
		      default: 500,
		      description: "The length of presentation for the fixation cross before either array is shown"
		    },
        between_fixation_dur: {
		      type: jsPsych.plugins.parameterType.INT,
		      pretty_name: "Fixation duration between arrays of the same stimulus",
		      default: 1000,
		      description: "The length of presentation for the fixation cross between arrays of the same stimulus"
		    },
        num_squares: {
          type: jsPsych.plugins.parameterType.INT,
          pretty_name: "Number of squares",
          default: 7,
          description: "Number of squares in the array. Must be less than or equal to colors.length * color_reps"
        },
        radius: {
          type: jsPsych.plugins.parameterType.INT,
          pretty_name: "Radius of the circular array ",
          default: 200,
          description: "Radius of the circular array in which the squares are on"
        },
        square_len: {
          type: jsPsych.plugins.parameterType.INT,
          pretty_name: "Square side length ",
          default: 50,
          description: "Side length of te squares that are in the circular array"
        },
        colors: {
          type: jsPsych.plugins.parameterType.ARRAY,
          pretty_name: "Colors of Squares",
          default: ["Gold", "DarkTurquoise", "LightGreen", "LightPink", "Tomato", "DodgerBlue", "LightGray", "LightSalmon"],
          description: "List of colors that the squares can be "
        },
        color_reps: {
          type: jsPsych.plugins.parameterType.INT,
          pretty_name: "Max repitions for a color",
          default: 1,
          description: "Max number of times a certain colored square can appear in the array"
        },
        background_color: {
		      type: jsPsych.plugins.parameterType.STRING,
		      pretty_name: "Background color",
		      default: "gray",
		      description: "The color of the 'canvas' that the stimuli will be drawn on"
	    	},
        fixation_color: {
		      type: jsPsych.plugins.parameterType.STRING,
		      pretty_name: "Color of fixation Cross",
		      default: "white",
		      description: "The color of the fixation cross that is in the middle of the circular array"
	    	},
	 }
}

	//BEGINNING OF TRIAL
	plugin.trial = function(display_element, trial) {

		// Assign unassigned parameters to their default
		function assignParameterValue(argument, defaultValue){
			return typeof argument !== 'undefined' ? argument : defaultValue;
		}
		trial.stim_duration = assignParameterValue(trial.stim_duration, 1000);
    trial.beginning_fixation_dur = assignParameterValue(trial.beginning_fixation_dur, 500);
    trial.between_fixation_dur = assignParameterValue(trial.between_fixation_dur, 1000);
    trial.num_squares = assignParameterValue(trial.num_squares, 7);
    trial.radius = assignParameterValue(trial.radius, 200);
    trial.square_len = assignParameterValue(trial.square_len, 50);
    trial.colors = assignParameterValue(trial.colors, ["Gold", "DarkTurquoise", "LightGreen", "LightPink", "Tomato", "DodgerBlue", "LightGray", "LightSalmon"]);
    trial.background_color = assignParameterValue(trial.background_color, "gray");
    trial.fixation_color = assignParameterValue(trial.fixation_color, "white");

		// START MAKING STIMULI (ie canvas)
		// Edit the body, currently described by <body class="jspsych-display-element"> .... </body>
		var body = document.getElementsByClassName("jspsych-display-element")[0];
		var originalMargin = body.style.margin;
		var originalPadding = body.style.padding;
		var originalBackgroundColor = body.style.backgroundColor;
		var originalCursor = body.style.cursor;

		body.style.margin = 0;
	  body.style.padding = 0;
		body.style.backgroundColor = trial.background_color;
		body.style.cursor = "none"

		//Create a canvas element and append it to the DOM
		var canvas = document.createElement("canvas");
		display_element.appendChild(canvas);

		//Remove the margins and padding of the canvas
		canvas.style.margin = 0;
		canvas.style.padding = 0;

		//Get the context of the canvas so that it can be painted on.
		var ctx = canvas.getContext("2d");
		canvas.width = window.innerWidth;
		canvas.height = window.innerHeight;

    // HELPER  FUNCTIONS THAT GENERATE THE STIMULI
    // randomly outputs a key of an object
    function rand_key(obj) {
      var keys = Object.keys(obj)
      return keys[ keys.length * Math.random() << 0];
    }

    // shuffles all available colors, returns [array of colors, object of remaining colors]
    function selective_shuffle(colors, reps){
      options = {}
      for (i = 0; i <colors.length; i++){
        options[colors[i]]= reps
      }

      shuffled = [0]
      while (shuffled.length < trial.num_squares){
        potential = rand_key(options)
        if (potential != shuffled[shuffled.length - 1]){
          shuffled.push(potential)
          left = options[potential] - 1
          if (left == 0){
            delete options[potential];
          } else {
            options[potential] = left
          }
        }
      }

      last = []
      while (last.length == 0){
        first = shuffled[1]
        potential = rand_key(options)
        if (potential != first && potential!= shuffled[shuffled.length - 1]){
          last.push(potential)
          //update options
          left = options[potential] - 1
          if (left == 0){
            delete options[potential];
          } else {
            options[potential] = left
          }
        }

      }
      return [shuffled.slice(1,shuffled.length).concat(last), options]
    }

    // makes deep copy of array / object
    function deep_copy(o) {
       var output, v, key;
       output = Array.isArray(o) ? [] : {};
       for (key in o) {
           v = o[key];
           output[key] = (typeof v === "object") ? copy(v) : v;
       }
       return output;
    }

    // changes one elem of color array for the 'different' task
    // [return the change, index of change, color it originally was, color it changes to]
    function change_one(selective_shuff_output){
        colors_2 = deep_copy(selective_shuff_output[0])
        options = deep_copy(selective_shuff_output[1])

        var i_oi = Math.floor(Math.random() * trial.num_squares);
        if (i_oi === trial.num_squares - 1){
          i_after = 0
          i_before = i_oi - 1
        }
        else if (i_oi === 0){
          i_before = trial.num_squares - 1
          i_after = i_oi + 1
        }
        else {
          i_before = i_oi - 1
          i_after = i_oi + 1
        }

        delete options[colors_1[i_oi]]
        delete options[colors_1[i_before]]
        delete options[colors_1[i_after]]

        new_colour = rand_key(options)

        old_colour = colors_2[i_oi]
        colors_2[i_oi] = new_colour

        return [colors_2, i_oi, old_colour, new_colour]
    };

    // DETERMINE COLOR ARRAYS
    [colors_1, opts_1] = selective_shuffle(trial.colors, trial.color_reps)
    if (trial.task_type === 'different'){
      [colors_2, changed_i, old_color, new_color] = change_one([colors_1, opts_1])
    }
    else{
      colors_2 = deep_copy(colors_1)
			changed_i = null
			old_color = null
			new_color = null
    }

    // CUSTOMIZE STIMULI
    var x_midpoint = 1/2*window.innerWidth
    var y_midpoint = 1/2*window.innerHeight

    function generate_beginning_fixation(){
      ctx.clearRect(0, 0, window.innerWidth, window.innerHeight);
      ctx.strokeStyle = trial.fixation_color;
      ctx.lineWidth = 5;
      ctx.beginPath();
      ctx.moveTo(x_midpoint - 20, y_midpoint);
      ctx.lineTo(x_midpoint + 20, y_midpoint);

      ctx.moveTo(x_midpoint, y_midpoint - 20);
      ctx.lineTo(x_midpoint, y_midpoint + 20);
      ctx.stroke();
    }

    function generate_stim(colors) {
      ctx.lineWidth = 5;
      ctx.beginPath();
      ctx.moveTo(x_midpoint - 20, y_midpoint);
      ctx.lineTo(x_midpoint + 20, y_midpoint);

      ctx.moveTo(x_midpoint, y_midpoint - 20);
      ctx.lineTo(x_midpoint, y_midpoint + 20);
      ctx.stroke();

    	for (box_idx = 0; box_idx < trial.num_squares; box_idx++) {
    	    radians = 2 * Math.PI * box_idx / trial.num_squares;
    	    x_pos = x_midpoint + Math.sin(radians) * trial.radius
    	    y_pos = y_midpoint + Math.cos(radians) * trial.radius

    	    ctx.beginPath();

    	    for (i = 0; i < trial.num_squares; i++){
    	      if (box_idx == i){
    	        ctx.fillStyle = colors[i];
    	      }
    	    }
    	    ctx.rect(x_pos - trial.square_len/2 , y_pos - trial.square_len/2 , trial.square_len, trial.square_len);
    	    ctx.fill();
    		}
      jsPsych.pluginAPI.setTimeout(function() {
          generate_between_fixation();
    		}, trial.stim_duration);

    	};

    function generate_between_fixation(){
      ctx.clearRect(0, 0, window.innerWidth, window.innerHeight);
      ctx.lineWidth = 5;
      ctx.beginPath();
      ctx.moveTo(x_midpoint - 20, y_midpoint);
      ctx.lineTo(x_midpoint + 20, y_midpoint);

      ctx.moveTo(x_midpoint, y_midpoint - 20);
      ctx.lineTo(x_midpoint, y_midpoint + 20);
      ctx.stroke();

      jsPsych.pluginAPI.setTimeout(function() {
        generate_stim2(colors_2);
  		}, trial.between_fixation_dur);

    }

    function generate_stim2(colors) {
      generate_stim(colors_2)
      jsPsych.pluginAPI.setTimeout(function() {
        end_trial();
  		}, trial.stim_duration);
    };

    generate_beginning_fixation()

		jsPsych.pluginAPI.setTimeout(function() {
      generate_stim(colors_1);
		}, trial.beginning_fixation_dur);


		//Function to end the trial properly
		function end_trial() {
			jsPsych.pluginAPI.clearAllTimeouts();
			//Place all the data to be saved from this trial in one data object
			var trial_data = {
        "task_type": trial.task_type,
        "colors_1": colors_1.toString(),
        "colors_2": colors_2.toString(),
        "changed_index": changed_i,
        "old_color": old_color,
        "new_color": new_color
			}

			//Remove the canvas as the child of the display_element element
			display_element.innerHTML='';

			//Restore the settings to JsPsych defaults
			body.style.margin = originalMargin;
			body.style.padding = originalPadding;
			body.style.backgroundColor = originalBackgroundColor
			body.style.cursor = originalCursor

			jsPsych.finishTrial(trial_data); //End this trial and move on to the next trial

		} //End of end_trial function

	};
	return plugin; //Return the plugin object which contains the trial
})();
