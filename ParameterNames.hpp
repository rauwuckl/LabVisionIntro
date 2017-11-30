// define string names of the parameters, to avoid misspelling strings (-> compiler will complain not at runtime)
#define EXPERIMENT_NAME "experiment_name"
#define EPOCH_COUNT "epoch_count"
#define PRELOAD_WEIGHTS "preload_weights"
#define TEST_EVERY_N_EPOCHS "test_every_n_epochs"
#define RECORD_SPIKES_IN_TRAINING "record_spikes_in_training"
#define START_EPOCH "start_epoch"
#define ALL_LAYER_SAME_FLAG "all_layer_same_parameters"

#define EXPERIMENT_SPECS "experiment_specs"
	//part of experiment_specs
	#define NETWORK_PARAMS_FILENAME "network_params_filename"
	#define TRAINING_STIMULI_LIST "training_stimuli_list"
	#define TESTING_STIMULI_LIST "testing_stimuli_list"
	#define LEARNING_RATE_INC_STOP "stop_lr_inc_epoch"

#define NETWORK_PARAMS "network_params"// theses are actually internally stored in a separte ptree for no great reason whatsoever
  #define STDP_TAU_C "STDP_Tau_C"
  #define STDP_TAU_D "STDP_Tau_D"

	#define LAYERWISE_BIO_CONDUCTANCE_SCALING "layerwise_bio_conductance_scaling"
		#define E2E_FF "e2e_ff"
		#define E2I_L "e2i_l"
		#define I2E_L "i2e_l"
