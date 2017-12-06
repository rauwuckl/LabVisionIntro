# Lab Vision Model Starter Pack

## Running the model:


	Stuff changed so no it can be run for example like:
	> ./Model -n white_diamond_l_vs_r_new_inh -s ../expSpecs/white_diamond_left_vs_right.json -b 40 -q 80 -e 200 -t 1 -p ../netParams/ALS_more_syn_new_inh.json

  -n foldername
	-s json file containing the training and testing list and net parameters
	-b start with epoch 40
	-q stop learning rate increasement at epoch 80
	-e train for 200 epochs
	-t test every epoch
	-p overwrite parameters with this file (has to come after -s to overwrite those params)


## Folders:

### ./Data
	This folder is where the example Gisi_Model looks for the input rate based image inputs.
./Data/MatlabGaborFilter
	All matlab scripts required to create visual inputs are located here. In particular the __filterImageSet.m__ file begins with a few lines describing how to filter a set of images into suitable input for the spiking neural networks.
./Data/Inputs_Gisi_BO
	The Gisi BO folder contains the set of stimuli that Gisbert used to investigate border ownership cells. This is just a set of examples but can give you an idea of how things are organised. This is also the folder used by the example script.

### ./Spike
	The spike folder contains the specific version of the Spike simulator that you are using. It is recommended (for this model at least) to use the nas_fix branch of spike. In order to get this version, the steps are:
	1) Ensure you have git (version control software) installed on your machine. This should already be there but if not, ask.
	2) The Spike software is hosted at: https://github.com/OFTNAI/Spike. In order to get a copy, go to your LabVisionIntro folder and run:
		>>> git clone https://github.com/OFTNAI/Spike.git
	3) After the software has been downloaded, you must change the current branch in which it is operating:
		>>> cd Spike
		>>> git checkout nas_fix
	4) It is now recommended that you also turn off the error checking in the Spike code. This ensures that the speed of the simulator is high. To do this, modify the file located at: "Spike/Spike/Backend/CUDA/Helpers/ErrorCheck.hpp". Now comment out (by placing "//" at the start of the line), line number 10. I.e. "#define CUDA_ERROR_CHECK" -> "//#define CUDA_ERROR_CHECK".
		N.B. if you decide to change any Spike code, it is recommended that you then turn this back on to check for errors.

### ./Results
	A folder created to copy results into. A set of example results from running 100 epochs is given.


## Files:

### ./CMakeLists.txt
	This file contains the details required by cmake (a tool designed to help package up software) to build a model with Spike
	-> This file should be edited if you want to compile/use a new model file (or if you want to change the model name). On line 21, this file has a list of the models which will be available to be build. If (for example) you want to add a new model, just add its name here on a new line.


### Gisi_Model.cpp
	The example model of a four layer visual model with a whole range of parameters. This model has been optimised for the BO task and will likely need re-optimisation for any alternative tasks.

### Parameters.md
	A description of the set of parameters which affect the example model. These are named according to the variable which you can look for in the model example and are then followed by a description of their function.
