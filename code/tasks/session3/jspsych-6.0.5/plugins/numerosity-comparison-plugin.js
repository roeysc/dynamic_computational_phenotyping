jsPsych.plugins["numerosity-task"] = (function() {
	var plugin = {};
	plugin.info = {
	    name: "numerosity-task",
	    parameters: {
				choices: {
	        type: jsPsych.plugins.parameterType.KEYCODE,
	        array: true,
	        pretty_name: 'choices',
	        default: [81, 80], // p and q for left and right, respectively
	        description: '[key to choose the left option, key to choose the right option]'
	      },
				trial_duration: {
		      type: jsPsych.plugins.parameterType.INT,
		      pretty_name: "Trial duration",
		      default: 1000,
		      description: "The length of stimulus (two assortments of circles) presentation"
		    },
				fixation_duration: {
		      type: jsPsych.plugins.parameterType.INT,
		      pretty_name: "Fixation cross duration",
		      default: 500,
		      description: "The length of fixation cross presentation between trials"
		    },
				response_duration: {
		      type: jsPsych.plugins.parameterType.INT,
		      pretty_name: "Feedback duration",
		      default: 500,
		      description: "The length of time that the aperture will be shown post response"
		    },
		    mean: {
		      type: jsPsych.plugins.parameterType.INT,
		      pretty_name: "Distribution's Mean",
		      default: 24,
		      description: "The mean number of circles. The number of circles for one of the arrays will be defined by this distribution."
		    },
				sd:{
					type: jsPsych.plugins.parameterType.INT,
		      pretty_name: "Distribution's Standard Deviation",
		      default: 5,
		      description: "Standard Deviation of circles in distribution"
				},
				radius:{
					type: jsPsych.plugins.parameterType.INT,
		      pretty_name: "radius of circles",
		      default: 10,
		      description: "radius of circles"
				},
		    background_color: {
		      type: jsPsych.plugins.parameterType.STRING,
		      pretty_name: "Background color",
		      default: "black",
		      description: "The color of the 'canvas' that the stimuli will be drawn on"
	    	},
				circle_color: {
		      type: jsPsych.plugins.parameterType.STRING,
		      pretty_name: "Color of Circles",
		      default: "white",
		      description: "Color of the circles"
		    },
				border_color: {
					type: jsPsych.plugins.parameterType.STRING,
					pretty_name: "Color of border",
					default: "white",
					description: "Color of the rectangular border that separates the two arrays"
				},
				border_width: {
					type: jsPsych.plugins.parameterType.INT,
					pretty_name: "width of border",
					default: 4,
					description: "The 'thickness' that the rectangular borders are drawn"
				},
				cursor_type: {
		      type: jsPsych.plugins.parameterType.STRING,
		      pretty_name: "Type of Cursor",
		      default: "none",
		      description: "Specifies HTML attribute for cursor"
		    },
				control_type: {
		      type: jsPsych.plugins.parameterType.STRING,
		      pretty_name: "control type",
		      default: "none",
		      description: "Specifies how the dot array should be controlled by, SA or radius"
		    }
	 }
}

	//BEGINNING OF TRIAL
	plugin.trial = function(display_element, trial) {

		// Assign unassigned parameters to their default
		function assignParameterValue(argument, defaultValue){
			return typeof argument !== 'undefined' ? argument : defaultValue;
		}

		trial.choices = assignParameterValue(trial.choices, [80, 81]);
		trial.trial_duration = assignParameterValue(trial.trial_duration, 1000);
		trial.fixation_duration = assignParameterValue(trial.fixation_duration, 500);
		trial.response_duration = assignParameterValue(trial.response_duration, 500);
		trial.mean= assignParameterValue(trial.mean, 65);
		trial.sd = assignParameterValue(trial.sd, 10);
		trial.radius = assignParameterValue(trial.radius, 10);
		trial.background_color = assignParameterValue(trial.background_color, "black");
		trial.circle_color = assignParameterValue(trial.circle_color, "white");
		trial.border_color = assignParameterValue(trial.border_color, "white");
		trial.border_width = assignParameterValue(trial.border_width, 4);
		trial.cursor_type = assignParameterValue(trial.cursor_type, "none");
		trial.control_type = assignParameterValue(trial.control_type, "none");

		// Set Default values for subject response
		// -1 is used if the trial times out and the subject has not pressed a valid key
		var response = {
			rt: -1,
			key: -1
		}

		// START MAKING STIMULI (ie canvas)
		// Edit the body, currently described by <body class="jspsych-display-element"> .... </body>
		var body = document.getElementsByClassName("jspsych-display-element")[0];
		var originalMargin = body.style.margin;
		var originalPadding = body.style.padding;
		var originalBackgroundColor = body.style.backgroundColor;
		var originalCursorType = body.style.cursor;

		body.style.margin = 0;
	  body.style.padding = 0;
		body.style.backgroundColor = trial.background_color;
		body.style.cursor = trial.cursor_type

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


		// DETERMINE SPECIFIC FEATURES OF STIMULI
		var circles_left; // draw a number from a normal distribution 1 for one array
		var circles_right; // draw a number from a normal distribution 2 for one array

		// Function to draw the number of circles to put in an array froma normal distribution with mean m and standard deviation sd
		function nrand(m,sd){
				if (sd == 0) {
						return m;
				} else {
				var x1, x2, rad, y1;
				do {
						x1 = 2 * Math.random() - 1;
						x2 = 2 * Math.random() - 1;
						rad = x1 * x1 + x2 * x2;
				} while(rad >= 1 || rad == 0);

				var c = Math.sqrt(-2 * Math.log(rad) / rad);
				var y = Math.round((x1 * c * (sd^2))+m);
				return y;
				}
		};

		function dist(x1, y1, x2, y2){
			a = x2-x1
			b = y2 -y1
			return Math.sqrt(a*a+b*b)
		}

		// ensure that the number of circles in the arrays are different
		while (circles_left === circles_right || circles_left < 0 || circles_right < 0){
			var circles_left = nrand(trial.mean, trial.sd)
			var circles_right = nrand(trial.mean, trial.sd)
		}

		var difference = circles_left - circles_right
		if (difference > 0){
			larger_magnitude = 'left'
		}
		else{
			larger_magnitude = 'right'
		}

		// determine radii of circles
		if (trial.control_type == 'SA'){
			if (Math.round(Math.random()) === 0){
				r1 = trial.radius
				r2 = r1 * Math.sqrt(circles_left/circles_right)
			} else {
				r2 = trial.radius
				r1 = r2 * Math.sqrt(circles_right/circles_left)
			}
		}
		else if (trial.control_type == 'radius'){
			r1 = trial.radius
			r2 = trial.radius
		}

		// Determining the size of the stimuli given the browser/screen/etc that the subject is using
		var x_midpoint = 1/2*window.innerWidth
		var y_midpoint = 1/2*window.innerHeight

		var rect_width = 350//1/4*window.innerWidth
		var rect_height = 300//1/2*window.innerHeight
		var rect_diff = 350

		var left_x = x_midpoint - 1/2*rect_diff - rect_width//1/12*window.innerWidth
		var left_y = y_midpoint - 1/2*rect_height//1/4*window.innerHeight

		var right_x = x_midpoint + 1/2*rect_diff  //2/3*window.innerWidth
		var right_y = y_midpoint - 1/2*rect_height//1/4*window.innerHeight


		generate_beginning_fixation()

		jsPsych.pluginAPI.setTimeout(function() {
      generate_circles();
		}, trial.fixation_duration);

		function generate_beginning_fixation(){
      ctx.clearRect(0, 0, window.innerWidth, window.innerHeight);
      ctx.strokeStyle = trial.circle_color;
      ctx.lineWidth = 5;
      ctx.beginPath();


      ctx.moveTo(x_midpoint - 20, y_midpoint);
      ctx.lineTo(x_midpoint + 20, y_midpoint);

      ctx.moveTo(x_midpoint, y_midpoint - 20);
      ctx.lineTo(x_midpoint, y_midpoint + 20);
      ctx.stroke();
    }


		function gen_left_centers() {
			var padding = r1 + trial.border_width
			var left_x_max = left_x + rect_width - padding
			var left_y_max = left_y + rect_height - padding

		  var circles = []
		      circle = {}
		      protection = 1000
		      counter = 0

		  while (circles.length < circles_left &&
		         counter < protection) {
		    circle = {
		      x: Math.random() * (left_x_max - left_x - padding) + left_x + padding,
		      y: Math.random() * (left_y_max - left_y - padding) + left_y + padding,
		    };
		    overlapping = false;

		    for (var i = 0; i < circles.length; i++) {
		      var existing = circles[i];
		      var d = dist(circle.x, circle.y, existing.x, existing.y)
		      if (d < 2*r1) { // overlapping
		        overlapping = true;
		        break;
		      }
		    }

		    if (!overlapping) {
		      circles.push(circle);
		    }

		    counter++;
		  }
			return circles
		}

		function gen_right_centers() {
			var padding = r2 + trial.border_width
			var right_x_max = right_x + rect_width - padding
			var right_y_max = right_y + rect_height - padding

		  var circles = []
		      circle = {}
		      protection = 1000
		      counter = 0

		  while (circles.length < circles_right &&
		         counter < protection) {
		    circle = {
		      x: Math.random() * (right_x_max - right_x - padding) + right_x + padding,
		      y: Math.random() * (right_y_max - right_y - padding) + right_y + padding,
		    };
		    overlapping = false;

		    for (var i = 0; i < circles.length; i++) {
		      var existing = circles[i];
		      var d = dist(circle.x, circle.y, existing.x, existing.y)
		      if (d < 2*r2) { // overlapping
		        overlapping = true;
		        break;
		      }
		    }

		    if (!overlapping) {
		      circles.push(circle);
		    }

		    counter++;
		  }
			return circles
		}
		// Draw the arrays of circles on the canvas
		function generate_circles(){
			var timerHasStarted = false; // reset the timer, don't start until after the circles are done being generated
			ctx.clearRect(0, 0, window.innerWidth, window.innerHeight); // ensure no past canvas features are remaining

			right_circs = gen_right_centers()
			for(var i=0; i < right_circs.length; i++){
					ctx.beginPath();
					ctx.arc(right_circs[i].x, right_circs[i].y, r2, 0 ,2*Math.PI);
					ctx.fillStyle = trial.circle_color;
					ctx.fill();
					ctx.closePath();
			}

			left_circs = gen_left_centers()
			for(var i=0; i < left_circs.length; i++){

					ctx.beginPath();
					ctx.arc(left_circs[i].x, left_circs[i].y, r1, 0 ,2*Math.PI);
					ctx.fillStyle = trial.circle_color;
					ctx.fill();
					ctx.closePath();
			}

			startKeyboardListener();
			jsPsych.pluginAPI.setTimeout(function() {
          end_trial();
    		}, trial.trial_duration);
		}

		// Determines what counts as a keyboard response, and how that response took
		function startKeyboardListener(){
			if (trial.choices != jsPsych.NO_KEYS) {
				//Create the keyboard listener to listen for subjects' key response
				keyboardListener = jsPsych.pluginAPI.getKeyboardResponse({
					callback_function: after_response, // Function to call once the subject presses a valid key
					valid_responses: trial.choices,
					rt_method: 'performance', // method to get rt
					persist: false, // false: keyboard listener will only trigger the first time a valid key is pressed. true: it has to be explicitly cancelled by the cancelKeyboardResponse plugin API.
					allow_held_key: false // false: Only register the key once, after this getKeyboardResponse function is called. (Check JsPsych docs for better info under 'jsPsych.pluginAPI.getKeyboardResponse').
				});
			}
		}

		function highlight_rect(key){
			ctx.beginPath();
			ctx.lineWidth = trial.border_width
			ctx.strokeStyle = trial.border_color;
			if (key == trial.choices[0]){
				ctx.rect(left_x, left_y, rect_width, rect_height);
				ctx.stroke();
			}
			else{
				ctx.rect(right_x, right_y, rect_width, rect_height);
				ctx.stroke();
			}
			ctx.closePath();
			jsPsych.pluginAPI.setTimeout(function() {
          end_trial()
    		}, trial.response_duration);
		}

		//Function to record the first response by the subject
		function after_response(info) {
			//If the response has not been recorded, record it
			if (response.key == -1) {
				response = info; //Replace the response object created above
				highlight_rect(response.key)
			}
			else{
				end_trial();

			}

		}

		//Function to end the trial properly
		function end_trial() {

			jsPsych.pluginAPI.clearAllTimeouts();
			//Kill the keyboard listener if keyboardListener has been defined
			if (typeof keyboardListener !== 'undefined') {
				jsPsych.pluginAPI.cancelKeyboardResponse(keyboardListener);
			}

			//Place all the data to be saved from this trial in one data object
			var trial_data = {
				"rt": response.rt, //The response time
				"key_press": response.key, //The key that the subject pressed
				"mean": trial.mean,
				"sd": trial.sd,
				'left_circles': circles_left,
				'right_circles': circles_right,
				'larger_magnitude': larger_magnitude,
				'left_radius': r1,
				'right_radius': r2
			}

			if (trial_data.key_press == trial.choices[0] && larger_magnitude == 'left'){
				trial_data.correct = true
			}
			else if (trial_data.key_press == trial.choices[1] && larger_magnitude == 'right'){
				trial_data.correct = true
			}
			else if (trial_data.key_press == -1 ){
				trial_data.correct = null
			}
			else{
				trial_data.correct = false
			}


			//Remove the canvas as the child of the display_element element
			display_element.innerHTML='';

			//Restore the settings to JsPsych defaults
			body.style.margin = originalMargin;
			body.style.padding = originalPadding;
			body.style.backgroundColor = originalBackgroundColor
			body.style.cursor = originalCursorType

			jsPsych.finishTrial(trial_data); //End this trial and move on to the next trial

		} //End of end_trial function

	};
	return plugin; //Return the plugin object which contains the trial
})();
