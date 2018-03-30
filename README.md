# Lab Vision Model Starter Pack
## Instalation

Download the reposetory and spike: 
> git clone --recursive -j8 https://github.com/rauwuckl/LabVisionIntro.git

Compile it:
First make a folder called 'Build', then switch into it:
`mkdir Build; cd Build` 
Now compile with
`cmake ..; make -j8`

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


## Spike Version
The version of Spike included is located at github.com/rauwuckl/Spike. 
It is a fork of the branch `nasfix` in github.com/OFTNAI/Spike. I added some details but it should still largely be compatible. 
