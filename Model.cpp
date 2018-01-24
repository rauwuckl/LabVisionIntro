// An Example Model for running the SPIKE simulator
//
// Authors: Nasir Ahmad (16/03/2016), James Isbister (23/3/2016), Gisbert Teepe [small position changes only](15/05/17)

// To create the executable for this network, run:

//inlcude define statements for parameter names
#include "ParameterNames.hpp"

#include "Spike/Models/SpikingModel.hpp"
#include "Spike/Simulator/Simulator.hpp"

#include <vector>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <iomanip>
#include <ctime>
#include <sstream>

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

using namespace std;
namespace pt = boost::property_tree;


/** Function to load weights from a file into a SpikingModel.
*/
void load_weights(
		  SpikingModel* Model,		/**< SpikingModel Pointer to the model which should load weights */
		  std::string weightloc,	/**< String path to the file from which weights should be loaded */
		  bool binaryfile)		/**< Boolean flag indicating if the file is a binary file */
{
	std::ifstream weightfile;
	std::vector<float> WeightsToLoad; // This vector should ultimately hold the list of replacement weights

	if (binaryfile){
		weightfile.open (weightloc, ios::in | ios::binary);
		while( weightfile.good() )
		{
			float currentweight;
			weightfile.read((char*)&currentweight, sizeof(float));
			if (weightfile.eof()) {
				weightfile.close();
				break;
			}
			WeightsToLoad.push_back(currentweight);
		}
	} else {
		weightfile.open(weightloc);
		while( weightfile.good() )
		{
			string fileline;
			getline( weightfile, fileline);
			if (weightfile.eof()) {
				weightfile.close();
				break;
			}
			// Put each line into the float vector
			WeightsToLoad.push_back( std::stof(fileline) );
		}
	}
	// Check if you have the correct number of weights
	if (WeightsToLoad.size() != Model->spiking_synapses->total_number_of_synapses){
		printf("%d, %d\n", (int)WeightsToLoad.size(), Model->spiking_synapses->total_number_of_synapses);
		printf("The number of weights being loaded is not equivalent to the model.");
		exit(2);
	}
	// If the previous test is passed, apply the weights
	for (int i=0; i < WeightsToLoad.size(); i++){
		Model->spiking_synapses->synaptic_efficacies_or_weights[i] = WeightsToLoad[i];
	}
	printf("%ld Weights Loaded.\n", WeightsToLoad.size());
}


/** Function to equalize the mean rate of the stimuli being presented to the network.
 *	Not strictly necessary if the stimuli are set up well.
*/
void equalize_rates(
			ImagePoissonInputSpikingNeurons* input_neurons, /**< ImagePoissonInputSpikingNeuron pointer to initialized input population */
			float target)					/**< float value indicating desired mean FR */
{
	// Rates can be altered here without much issue
	int num_rates_per_image = input_neurons->total_number_of_rates_per_image;
	int num_images = input_neurons->total_number_of_input_stimuli;

	for (int image_index = 0; image_index < num_images; image_index++){
		float meanval = 0.0f;
		for (int rate_index = 0; rate_index < num_rates_per_image; rate_index++){
			meanval += input_neurons->gabor_input_rates[image_index*num_rates_per_image + rate_index];
		}
		meanval /= float(num_rates_per_image);
		// printf("%f\n", meanval);

		float multiplication_factor = target / meanval;
		for (int rate_index = 0; rate_index < num_rates_per_image; rate_index++){
			input_neurons->gabor_input_rates[image_index*num_rates_per_image + rate_index] *= multiplication_factor;
		}
	}
}

std::string get_time_string(){
	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);

	std::ostringstream oss;
	oss << std::put_time(&tm, "%m_%d-%H_%M");
	std::string str = oss.str();
	return str;
}


void save_tree(std::string filename, pt::ptree tree){
		pt::write_json(filename + ".json", tree);
}

void load_tree(std::string filename, pt::ptree & tree){
		pt::read_json(filename, tree);
}

template <typename T>
std::vector<T> as_vector(const pt::ptree &pt, string key)
{
    std::vector<T> r;
    for (auto& item : pt.get_child(key)){
        r.push_back(item.second.get_value<T>());
		}
    return r;
}

template <typename T, unsigned L>
void assert_array_equals_vector(const T (&arr)[L], const std::vector<T> &vec){
	assert(L == vec.size());
	for(int i=0; i<vec.size(); i++){
		assert((vec[i] - arr[i]) < 1.0e-15);
	}
}



pt::ptree parse_argv(int argc, char *argv[]){
	/*
	*
	* parse the given options
	*/

  pt::ptree simulation_params;
	pt::ptree network_params; // seperate tree for all network parameters. (i.e. that define the network and not the experiment)

	auto update_network_params = [&network_params, &simulation_params](){
		std::string network_params_filename = simulation_params.get<string>(EXPERIMENT_SPECS "." NETWORK_PARAMS_FILENAME);
		std::cout << "Loading the file >" << network_params_filename << "< as network_params\n";
		load_tree(network_params_filename, network_params);
	};


	int option_index = 0;
	int opt = 0;
	static struct option long_options[]{
		{EXPERIMENT_NAME, required_argument, 0, 'n'}, //have to be given
		{EPOCH_COUNT, required_argument, 0, 'e'}, //have to be given
		//optional
		{PRELOAD_WEIGHTS, required_argument, 0, 'w'},
		{RECORD_SPIKES_IN_TRAINING, no_argument, 0, 'r'},
		{TEST_EVERY_N_EPOCHS, required_argument, 0, 't'},

		//either EXPERIMENT_SPECS or the other three have to be given
		{EXPERIMENT_SPECS, required_argument, 0, 's'},//file with the required specs if not individually given by the following 3 params
		//experiment_specs
		{NETWORK_PARAMS_FILENAME, required_argument, 0, 'p'}, //could be in a konfiguration file
		{TRAINING_STIMULI_LIST, required_argument, 0, 'l'},
		{TESTING_STIMULI_LIST, required_argument, 0, 'c'},
		{LEARNING_RATE_INC_STOP, required_argument, 0, 'q'},
		{START_EPOCH, required_argument, 0, 'b'},
		{ONLY_TEST_STIMULI_FOLDER, required_argument, 0, 'T'},
		{0,0,0,0}
	};

	while(true){
		opt = getopt_long(argc, argv, "n:e:w:rt:s:p:l:c:q:b:T:", long_options, &option_index);
		if(-1 == opt) break;
		switch (opt){
			case 'n':
				simulation_params.put(EXPERIMENT_NAME, optarg);
				break;
			case 'T':
				simulation_params.put(ONLY_TEST_STIMULI_FOLDER, optarg);
				break;
			case 'e':
				simulation_params.put(EPOCH_COUNT, atoi(optarg));
				break;
			case 'w':
				simulation_params.put(PRELOAD_WEIGHTS, optarg);
				break;
			case 't':
				simulation_params.put(TEST_EVERY_N_EPOCHS, atoi(optarg));
				break;
			case 'r':
				simulation_params.put(RECORD_SPIKES_IN_TRAINING, 1);
				break;
			case 'b':
				simulation_params.put(START_EPOCH, atoi(optarg));
				break;
			case 's':
			{
				std::string experiment_specs_filename = optarg;
				if(simulation_params.get_optional<string>(EXPERIMENT_SPECS)){
					std::cout << "some experiment specifications have already been given. These will be overwritten by the file: "<< experiment_specs_filename <<"\n" << std::endl;
				}
				std::cout << "Loading experiment specifications (i.e. the network_params, training_stimuli_list and testing_stimuli_list) from the file >" << experiment_specs_filename << "<\n"<< std::endl;
				pt::ptree experiment_specs;
				load_tree(experiment_specs_filename, experiment_specs);
				simulation_params.add_child(EXPERIMENT_SPECS , experiment_specs);
				simulation_params.put(EXPERIMENT_SPECS ".expSpecs_filename" , optarg); // add name of the experiment specs file
				update_network_params();
				break;
			}
			case 'p':
				simulation_params.put(EXPERIMENT_SPECS "." NETWORK_PARAMS_FILENAME, optarg);
				update_network_params();
				break;
			case 'l':
			  simulation_params.put(EXPERIMENT_SPECS "." TRAINING_STIMULI_LIST, optarg);
			  break;
			case 'c':
			  simulation_params.put(EXPERIMENT_SPECS "." TESTING_STIMULI_LIST, optarg);
				break;
			case 'q':
			  simulation_params.put(EXPERIMENT_SPECS "." LEARNING_RATE_INC_STOP, atoi(optarg));
				break;
			default:
				std::cout << "Unkown option given" << std::endl;
				break;
		}
	}
	bool input_ok = true;
	if(simulation_params.not_found() == simulation_params.find(EXPERIMENT_NAME)){
		std::cout << "Did not find parameter " EXPERIMENT_NAME << "\n" << std::endl;
		input_ok = false;
	}

	if(simulation_params.get_optional<string>(ONLY_TEST_STIMULI_FOLDER)){
		if(!simulation_params.get_optional<string>(EXPERIMENT_SPECS "." TESTING_STIMULI_LIST)){
			cout << TESTING_STIMULI_LIST << " has to be provided \n"<<endl;
			input_ok = false;
		}
		if(!simulation_params.get_optional<string>(EXPERIMENT_SPECS "." NETWORK_PARAMS_FILENAME)){
			cout << NETWORK_PARAMS << " has to be provided \n"<<endl;
			input_ok = false;
		}
	}
	else{
		//actuall training
		if(simulation_params.not_found() == simulation_params.find(EPOCH_COUNT)){
			std::cout << "Did not find parameter " EPOCH_COUNT << std::endl;
			if(simulation_params.not_found() != simulation_params.find(PRELOAD_WEIGHTS)){
				std::cout << "Running only testing with the Preloaded weight file \n" << std::endl;
				simulation_params.put(EPOCH_COUNT, -1);
			}
			else input_ok = false;
		}
		if(simulation_params.not_found() == simulation_params.find(TEST_EVERY_N_EPOCHS)){
			std::cout << "Did not find parameter " TEST_EVERY_N_EPOCHS " using default value 1 \n"<< std::endl;
			simulation_params.put(TEST_EVERY_N_EPOCHS, 1);
		}

		if(!simulation_params.get_optional<int>(START_EPOCH)){
			simulation_params.put(START_EPOCH, 1);
		}


		if(
			!simulation_params.get_optional<string>(EXPERIMENT_SPECS "." NETWORK_PARAMS_FILENAME)||
			!simulation_params.get_optional<string>(EXPERIMENT_SPECS "." TRAINING_STIMULI_LIST)||
			!simulation_params.get_optional<string>(EXPERIMENT_SPECS "." TESTING_STIMULI_LIST)||
			!simulation_params.get_optional<int>(EXPERIMENT_SPECS "." LEARNING_RATE_INC_STOP)
		){
			std::cout << "One of " NETWORK_PARAMS_FILENAME " or " TRAINING_STIMULI_LIST " or " TESTING_STIMULI_LIST " or " LEARNING_RATE_INC_STOP " was not defined" << std::endl;
		  input_ok = false;
		}
	}


	if(!input_ok){exit(EXIT_FAILURE);}

	std::cout << "Tau C: " << network_params.get<string>(STDP_TAU_C) << std::endl;
	simulation_params.add_child(NETWORK_PARAMS, network_params);
	return simulation_params;
}

	/*
	 *
	 *	General Simulation Settings and Parameters
	 *
	 */

/*
 *	Main function in which the network is created and run
 */
int main (int argc, char *argv[]){


	pt::ptree simulation_params = parse_argv(argc, argv);
	std::cout << "parsing params done" << std::endl;


	std::string modelpath = "../";
	std::string modelname = "Model.cpp";
	std::string output_folder_path = "../output/";

	// Files/Paths relevent to the input set
	// std::string filepath_stimuli = "../Data/MatlabGaborFilter/Inputs_Gisi_BO/";
	std::string stimuli_folder = "training";

	if(simulation_params.get_optional<string>(EXPERIMENT_SPECS "." STIMULI_FOLDER)){
		stimuli_folder = simulation_params.get<string>(EXPERIMENT_SPECS "." STIMULI_FOLDER);
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cout << "!! using unusual stimuli folder: "<< stimuli_folder << endl;
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	}

	if(simulation_params.get_optional<string>(ONLY_TEST_STIMULI_FOLDER)){
		stimuli_folder = simulation_params.get<string>(ONLY_TEST_STIMULI_FOLDER);
	}

	std::string filepath_stimuli = "../Data/MatlabGaborFilter/" + stimuli_folder + "/";

	// Finally the experiment name is set up to be a combination of the stimulation name and number of epochs
	string experimentName = get_time_string() + string("_") + simulation_params.get<std::string>(EXPERIMENT_NAME);
	std::cout << "Starting experiment: " << experimentName << std::endl;

	/*
	*
	* Save All the loaded data
	*
	*/
	// Creating Relevant folders
	Simulator::CreateDirectoryForSimulationDataFiles(experimentName, output_folder_path);
	std::cout << "Created folder " << output_folder_path << "  " <<  experimentName << std::endl;
	// Copy the model to your output folder (for future reference)
	std::string source = modelpath + modelname;
	std::string destination = output_folder_path + experimentName+"/";

	//copy the .cpp file
	ifstream srce(source.c_str(), ios::binary ) ;
	ofstream dest((destination + modelname).c_str(), ios::binary ) ;
	dest << srce.rdbuf() ;
	dest.close();
	srce.close();

	if(!simulation_params.get_optional<string>(ONLY_TEST_STIMULI_FOLDER)){
		// if we are not only testing
		//copy train list
		source = filepath_stimuli + simulation_params.get<string>(EXPERIMENT_SPECS "." TRAINING_STIMULI_LIST);
		srce = ifstream(source.c_str(), ios::binary );
		dest = ofstream((destination + "training_list.txt").c_str(), ios::binary);
		dest << srce.rdbuf();
		dest.close();
		srce.close();
	}


	//copy test list
	source = filepath_stimuli + simulation_params.get<string>(EXPERIMENT_SPECS "." TESTING_STIMULI_LIST);
	srce = ifstream(source.c_str(), ios::binary );
	dest = ofstream((destination + "testing_list.txt").c_str(), ios::binary);
	dest << srce.rdbuf();
	dest.close();
	srce.close();


	save_tree((output_folder_path + experimentName + "/simulation_params"), simulation_params);


	/*
	* CORE Params that wont change very often
	*/
	// The timestep at which this network is run is entered here.
	// Note that the timestep can be set higher than usual for any period in which you want to generally test the network behaviour.
	const float presentation_time_per_stimulus_per_epoch_test = 2.0f; // seconds
	float presentation_time_per_stimulus_per_epoch_train = 0.2f; // 4.0f didn't really work; //used to be 2.0//0.2;//2.0f; // seconds

	float timestep = 0.0002;
	// float original_timestep = 0.00002;		// This value is the timestep used in Aki's spiking investigations
	bool human_readable_storage = false;

	// Since the model can be run under different connectivity styles, these booleans turn them on/off
	bool E2E_FB_ON = true;
	bool E2E_L_ON = true;
	bool E2E_L_STDP_ON = true;

	// In order to set up a sensible set of FF exc and inh values, a set of booleans have been set up to turn on/off the values
	bool inh_layer_on[] = {true, true, true, true};

	// Parameters for testing

	/*
	 *
	 *	Visual Model General Settings
	 *
	 */

	// Network Parameters
	const int number_of_layers = 4;			// This value is explicitly assumed in this model. Not recommended to change unless you understand what else may need changing in this file.
	// int max_number_of_connections_per_pair = 2; // (redundant with the fan in Counts) The maximum number of connections refers to multiple synaptic contacts pre->post
	int dim_excit_layer = 64;			// The dimension of the excitatory layers (grid with this width)
	int dim_inhib_layer = 32;			// The dimension of the inhibitory layers (as above)

	// G2E = Gabor to excitatory, E2E = excitatory to excitatory, E2I = excitatory to inhibitory, I2E = inhibitory to excitatory
	// FF = feed forward, L = Lateral, FB = Feedback

	// Measure of the radius of the Fan-in
	// float gaussian_synapses_standard_deviation_G2E_FF =  4.0;
	// float gaussian_synapses_standard_deviation_E2E_FF[number_of_layers-1] = {50.0, 50.0, 50.0}; // List for each layer, can be customized Seems VERY high
	float gaussian_synapses_standard_deviation_E2E_FB = 8.0;
	float gaussian_synapses_standard_deviation_E2E_L = 14.0;
	// float gaussian_synapses_standard_deviation_E2I_L = 4.0; // one inhibitory neuron depends on very few excitatory neurons havily (higher chance of drawing the same pre neuron)
	// float gaussian_synapses_standard_deviation_I2E_L = 8.0;

	// Fan-in Number
	// int fanInCount_G2E_FF = 90;
	// int fanInCount_E2E_FF = 90; // replaced by a json parameter, actually has to be double because of the max_number_of_connections_per_pair=2
	int fanInCount_E2E_FB = 5;
	int fanInCount_E2E_L = 30;
	// int fanInCount_E2I_L = 60;
	// int fanInCount_I2E_L = 90; // <<<< used to be 90 // this seeems high. inhibition from 90 neurons seems high

	// if (fanInCount_E2E_FF % max_number_of_connections_per_pair!=0){
	// 	printf("total_number_of_new_synapses has to be a multiple of max_number_of_connections_per_pair");
	// 	return 0;
	// }

	// Synaptic Parameters
	// Range of axonal transmission delay
	// timestep is defined above
	float min_delay = 5.0*timestep; // In timesteps
	float max_delay = 0.01; // In seconds (10ms)
	float max_FR_of_input_Gabor = 100.0f; // Hz
	float absolute_refractory_period = 0.002; // s
	//Synaptic Parameters
	float weight_range_bottom = 0.0;
	float weight_range_top = 1.0;
	float learning_rate_rho = 0.001f;

	// int stop_lr_increase_epoch = 75;

	// calculating different Connections
	float E2E_FF_minDelay = min_delay;
	float E2E_FF_maxDelay = max_delay;//3.0f*pow(10, -3);
	float E2I_L_minDelay = min_delay;
	float E2I_L_maxDelay = max_delay;//3.0f*pow(10, -3);
	float I2E_L_minDelay = min_delay;
	float I2E_L_maxDelay = max_delay;//3.0f*pow(10, -3);
	float E2E_FB_minDelay = min_delay;
	float E2E_FB_maxDelay = max_delay;
	float E2E_L_minDelay = min_delay;
	float E2E_L_maxDelay = max_delay;

	// Below are the decay rates of the variables for learning: Pre/Post synaptic activities C and D (See Ben Evans)
	// float decay_term_tau_C = 0.05; //aki_paper = 0.005
	// float decay_term_tau_D = 0.05; //aki_paper = 0.005

	// Biological Scaling Constant = How much you multiply the weights up or down for realism/stability
	// If this value is roughly on the order of the Leakage Conductance, it will be close to one input spike -> one output spike (n.b. depends on syn tau)
	float biological_conductance_scaling_constant_lambda_G2E_FF = 0.1 * 0.0001 * 0.00002; // << why product #TODO
	// float old_biological_conductance_scaling_constant_lambda_E2E_FF = 0.00005 * 0.00002;
	float biological_conductance_scaling_constant_lambda_E2E_FB = 0.1 * 0.0001 * 0.00002;
	float biological_conductance_scaling_constant_lambda_E2E_L	= 0.000001 * 0.00002;
	// float biological_conductance_scaling_constant_lambda_E2I_L	= 0.001 * 0.00002;
	// float biological_conductance_scaling_constant_lambda_I2E_L	= 0.005 * 0.00002;


	// float layerwise_biological_conductance_scaling_constant_lambda_E2E_FF[number_of_layers-1] = {
	// 	0.625f * biological_conductance_scaling_constant_lambda_E2E_FF,
	// 	0.5f * biological_conductance_scaling_constant_lambda_E2E_FF,
	// 	0.75f * biological_conductance_scaling_constant_lambda_E2E_FF};//different for the different layers


	// check with the loaded one
	std::vector<float> layerwise_conductance_scaling_E2E_FF = as_vector<float>(simulation_params, NETWORK_PARAMS "." LAYERWISE_BIO_CONDUCTANCE_SCALING "." E2E_FF);
	// assert_array_equals_vector(layerwise_biological_conductance_scaling_constant_lambda_E2E_FF, layerwise_conductance_scaling_E2E_FF);

	// float layerwise_biological_conductance_scaling_constant_lambda_E2I_L[number_of_layers] = {
	// 	1.1f * biological_conductance_scaling_constant_lambda_E2I_L,
	// 	1.625f * biological_conductance_scaling_constant_lambda_E2I_L,
	// 	0.875f * biological_conductance_scaling_constant_lambda_E2I_L, // Layer 2 different and better performance
	// 	1.6f * biological_conductance_scaling_constant_lambda_E2I_L};

	std::vector<float> layerwise_conductance_scaling_E2I_L = as_vector<float>(simulation_params, NETWORK_PARAMS "." LAYERWISE_BIO_CONDUCTANCE_SCALING "." E2I_L);
	// assert_array_equals_vector(layerwise_biological_conductance_scaling_constant_lambda_E2I_L, layerwise_conductance_scaling_E2I_L);

	// float layerwise_biological_conductance_scaling_constant_lambda_I2E_L[number_of_layers] = {
	// 	0.04f * biological_conductance_scaling_constant_lambda_I2E_L,
	// 	0.375f * biological_conductance_scaling_constant_lambda_I2E_L,
	// 	0.2f * biological_conductance_scaling_constant_lambda_I2E_L, // Layer 2 different and better performance
	// 	0.325f * biological_conductance_scaling_constant_lambda_I2E_L};

	std::vector<float> layerwise_conductance_scaling_I2E_L = as_vector<float>(simulation_params, NETWORK_PARAMS "." LAYERWISE_BIO_CONDUCTANCE_SCALING "." I2E_L);
	// assert_array_equals_vector(layerwise_biological_conductance_scaling_constant_lambda_I2E_L, layerwise_conductance_scaling_I2E_L);

	// Tau G = Synaptic Conductance Decay TIME CONSTANT for each synapse type (#1) (Seconds)
	// Most of these values are set to 150ms for trace-like learning. Other than Exc->Inh and Inh->Exc
	float decay_term_tau_g_G2E_FF	=	0.15;
	float decay_term_tau_g_E2E_FF	=	0.15;
	float decay_term_tau_g_E2E_FB	=	0.15;
	float decay_term_tau_g_E2E_L	=	0.15;
	float decay_term_tau_g_E2I_L	=	0.002;
	float decay_term_tau_g_I2E_L	=	0.025; //0.005;//In Ben's model, 0.005 v 0.025 and latter produced better result


	/*
	 *
	 *	Defining the Spiking Model
	 *
	 */

	// Create the SpikingModel
	SpikingModel* model = new SpikingModel();
	model->SetTimestep(timestep);

	LIFSpikingNeurons* lif_spiking_neurons = new LIFSpikingNeurons();
	ImagePoissonInputSpikingNeurons* input_neurons = new ImagePoissonInputSpikingNeurons();

	// Optionally make the input_neurons ignore the files and have random rates instead
	// I.e. if we set this flag each stimulus will be random noise, instead of what is actually in the file
	auto noise_stimuli_flag = simulation_params.get_optional<int>(EXPERIMENT_SPECS "." NOISE_STIMULI);
	if(noise_stimuli_flag && (*noise_stimuli_flag == 1)){
		input_neurons->make_stimuli_as_random_noise = true;
	}


	ConductanceSpikingSynapses* conductance_spiking_synapses = new ConductanceSpikingSynapses();

	model->spiking_neurons = lif_spiking_neurons;
	model->input_spiking_neurons = input_neurons;
	model->spiking_synapses = conductance_spiking_synapses;


	// STDP Rule Parameters
	evans_stdp_plasticity_parameters_struct STDP_PARAMS;
	STDP_PARAMS.decay_term_tau_C = simulation_params.get<float>(NETWORK_PARAMS "." STDP_TAU_C);
	STDP_PARAMS.decay_term_tau_D = simulation_params.get<float>(NETWORK_PARAMS "." STDP_TAU_D);
	STDP_PARAMS.model_parameter_alpha_D = 0.5;
	STDP_PARAMS.synaptic_neurotransmitter_concentration_alpha_C = 0.5*2.0f;//<< why product #TODO
	STDP_PARAMS.learning_rate_rho = learning_rate_rho;
	EvansSTDPPlasticity* evans_stdp = new EvansSTDPPlasticity(conductance_spiking_synapses, lif_spiking_neurons, input_neurons, &STDP_PARAMS);

	model->AddPlasticityRule(evans_stdp);


	model->init_backend();

	conductance_spiking_synapses->print_synapse_group_details = false;

	// Creating the input neurons
	TimerWithMessages adding_input_neurons_timer("Adding Input Neurons...\n");
	// Loading the required files
	input_neurons->set_up_rates(simulation_params.get<string>(EXPERIMENT_SPECS "." TESTING_STIMULI_LIST).c_str(), "FilterParameters.txt", filepath_stimuli.c_str(), max_FR_of_input_Gabor);

	image_poisson_input_spiking_neuron_parameters_struct image_poisson_input_spiking_group_params;
	image_poisson_input_spiking_group_params.rate = 30.0f;
	input_neurons->AddGroupForEachGaborType(&image_poisson_input_spiking_group_params);

	adding_input_neurons_timer.stop_timer_and_log_time_and_message("Input Neurons Added.", true);


	// Neuron Layer Creation
	TimerWithMessages adding_neurons_timer("Adding Neurons...\n");

	lif_spiking_neuron_parameters_struct EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS;
	EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS.group_shape[0] = dim_excit_layer;
	EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS.group_shape[1] = dim_excit_layer;
	EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS.resting_potential_v0 = -0.074f;
	EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS.threshold_for_action_potential_spike = -0.053f;
	EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS.somatic_capacitance_Cm = 500.0*pow(10, -12);
	EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS.somatic_leakage_conductance_g0 = 25.0*pow(10, -9); // << is this the leackage conductance that the
	EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS.absolute_refractory_period = absolute_refractory_period;


	lif_spiking_neuron_parameters_struct INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS;
	INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS.group_shape[0] = dim_inhib_layer;
	INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS.group_shape[1] = dim_inhib_layer;
	INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS.resting_potential_v0 = -0.082f;
	INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS.threshold_for_action_potential_spike = -0.053f;
	INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS.somatic_capacitance_Cm = 214.0*pow(10, -12);
	INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS.somatic_leakage_conductance_g0 = 18.0*pow(10, -9);
	INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS.absolute_refractory_period = absolute_refractory_period;

	vector<int> EXCITATORY_NEURONS;
	vector<int> INHIBITORY_NEURONS;
	for (int l=0;l<number_of_layers;l++){
		EXCITATORY_NEURONS.push_back(model->AddNeuronGroup(&EXCITATORY_LIF_SPIKING_NEURON_GROUP_PARAMS));
		INHIBITORY_NEURONS.push_back(model->AddNeuronGroup(&INHIBITORY_LIF_SPIKING_NEURON_GROUP_PARAMS));
		cout<<"Neuron Group "<<EXCITATORY_NEURONS[l]<<": Excitatory layer "<<l<<endl;
		cout<<"Neuron Group "<<INHIBITORY_NEURONS[l]<<": Inhibitory layer "<<l<<endl;
	}


	adding_neurons_timer.stop_timer_and_log_time_and_message("Neurons Added.", true);


	// Synapse Creation
	TimerWithMessages adding_synapses_timer("Adding Synapses...\n");


	conductance_spiking_synapse_parameters_struct G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
	G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.delay_range[0] = timestep;
	G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.delay_range[1] = timestep;
	G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.max_number_of_connections_per_pair = 1;
	G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.gaussian_synapses_per_postsynaptic_neuron = simulation_params.get<int>(NETWORK_PARAMS "." FAN_IN_COUNT "." G2E_FF);//fanInCount_G2E_FF;
	G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.biological_conductance_scaling_constant_lambda = biological_conductance_scaling_constant_lambda_G2E_FF;
	G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
	// In aki's model, learning on this set of synapses was off. Remove the line below to math that.
	G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.plasticity_vec.push_back(evans_stdp);
	G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.gaussian_synapses_standard_deviation = simulation_params.get<float>(NETWORK_PARAMS "." FAN_IN_STD "." G2E_FF);
	G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.reversal_potential_Vhat = 0.0; //Volts
	G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.decay_term_tau_g = decay_term_tau_g_G2E_FF;
	G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.weight_range_bottom = weight_range_bottom;
	G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.weight_range_top = weight_range_top;


	conductance_spiking_synapse_parameters_struct E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
	E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.delay_range[0] = E2E_FF_minDelay;//5.0*timestep;
	E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.delay_range[1] = E2E_FF_maxDelay;//10.0f*pow(10, -3);
	E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.max_number_of_connections_per_pair = 1;
	E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.gaussian_synapses_per_postsynaptic_neuron = simulation_params.get<int>(NETWORK_PARAMS "." FAN_IN_COUNT "." E2E_FF);//  fanInCount_E2E_FF;
	E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.biological_conductance_scaling_constant_lambda = -42; //this will be overwritten later biological_conductance_scaling_constant_lambda_E2E_FF;
	E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
	E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.plasticity_vec.push_back(evans_stdp);
	E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.reversal_potential_Vhat = 0.0;
	E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.decay_term_tau_g = decay_term_tau_g_E2E_FF;
	E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.weight_range_bottom = weight_range_bottom;
	E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.weight_range_top = weight_range_top;


	conductance_spiking_synapse_parameters_struct E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
	if(E2E_FB_ON){
		E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.delay_range[0] = E2E_FB_minDelay;
		E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.delay_range[1] = E2E_FB_maxDelay;
		E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.max_number_of_connections_per_pair = 1;
		E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.gaussian_synapses_per_postsynaptic_neuron = fanInCount_E2E_FB;
		E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.biological_conductance_scaling_constant_lambda = biological_conductance_scaling_constant_lambda_E2E_FB;
		E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
		E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.plasticity_vec.push_back(evans_stdp);
		E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.gaussian_synapses_standard_deviation = gaussian_synapses_standard_deviation_E2E_FB;
		E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.reversal_potential_Vhat = 0.0;
		E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.decay_term_tau_g = decay_term_tau_g_E2E_FB;
		E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.weight_range_bottom = weight_range_bottom;
		E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.weight_range_top = weight_range_top;
	}


	conductance_spiking_synapse_parameters_struct E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
	E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.delay_range[0] = E2I_L_minDelay; //5.0*timestep;
	E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.delay_range[1] = E2I_L_maxDelay; //10.0f*pow(10, -3);
	E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.max_number_of_connections_per_pair = 1;
	E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.gaussian_synapses_per_postsynaptic_neuron = simulation_params.get<int>(NETWORK_PARAMS "." FAN_IN_COUNT "." E2I_L);//fanInCount_E2I_L;
	E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.biological_conductance_scaling_constant_lambda = -42; // will be overwritten later biological_conductance_scaling_constant_lambda_E2I_L;
	E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
	E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.gaussian_synapses_standard_deviation = simulation_params.get<float>(NETWORK_PARAMS "." FAN_IN_STD "." E2I_L);//gaussian_synapses_standard_deviation_E2I_L;
	E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.reversal_potential_Vhat = 0.0;
	E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.decay_term_tau_g = decay_term_tau_g_E2I_L;
	E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.weight_range_bottom = weight_range_bottom;
	E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.weight_range_top = weight_range_top;

	conductance_spiking_synapse_parameters_struct I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
	I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.delay_range[0] = I2E_L_minDelay;//5.0*timestep;
	I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.delay_range[1] = I2E_L_maxDelay;//3.0f*pow(10, -3);
	I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.max_number_of_connections_per_pair = 1;
	I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.gaussian_synapses_per_postsynaptic_neuron = simulation_params.get<int>(NETWORK_PARAMS "." FAN_IN_COUNT "." I2E_L);// fanInCount_I2E_L;

	I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
	I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.gaussian_synapses_standard_deviation = simulation_params.get<float>(NETWORK_PARAMS "." FAN_IN_STD "." I2E_L);//gaussian_synapses_standard_deviation_I2E_L;
	I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.reversal_potential_Vhat = -70.0*pow(10, -3);
	I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.decay_term_tau_g = decay_term_tau_g_I2E_L;
	I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.weight_range_bottom = weight_range_bottom;
	I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.weight_range_top = weight_range_top;

	conductance_spiking_synapse_parameters_struct E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS;
	if(E2E_L_ON){
		E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.delay_range[0] = E2E_L_minDelay;
		E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.delay_range[1] = E2E_L_maxDelay;
		E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.max_number_of_connections_per_pair = 1;
		E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.gaussian_synapses_per_postsynaptic_neuron = fanInCount_E2E_L;
		E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.biological_conductance_scaling_constant_lambda = biological_conductance_scaling_constant_lambda_E2E_L;
		E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.connectivity_type = CONNECTIVITY_TYPE_GAUSSIAN_SAMPLE;
		if (E2E_L_STDP_ON)
		  E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.plasticity_vec.push_back(evans_stdp);
		E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.gaussian_synapses_standard_deviation = gaussian_synapses_standard_deviation_E2E_L;
		E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.reversal_potential_Vhat = 0.0;
		E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.decay_term_tau_g = decay_term_tau_g_E2E_L;
		E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.weight_range_bottom = weight_range_bottom;
		E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.weight_range_top = weight_range_top;
	}



	for (int l=0; l<number_of_layers; l++){
		if(l==0)
		  model->AddSynapseGroupsForNeuronGroupAndEachInputGroup(
				  EXCITATORY_NEURONS[l],
				  &G2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS);
		else{
			E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.gaussian_synapses_standard_deviation = simulation_params.get<float>(NETWORK_PARAMS "." FAN_IN_STD "." E2E_FF); //gaussian_synapses_standard_deviation_E2E_FF[l-1];
			E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.biological_conductance_scaling_constant_lambda = layerwise_conductance_scaling_E2E_FF[l-1];
			// for (int connection_number = 0; connection_number < max_number_of_connections_per_pair; connection_number++){
				model->AddSynapseGroup(EXCITATORY_NEURONS[l-1], EXCITATORY_NEURONS[l], &E2E_FF_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS);
			// }
			if(E2E_FB_ON)
				model->AddSynapseGroup(EXCITATORY_NEURONS[l], EXCITATORY_NEURONS[l-1], &E2E_FB_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS);
		}
		E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.biological_conductance_scaling_constant_lambda = layerwise_conductance_scaling_E2I_L[l];// layerwise_biological_conductance_scaling_constant_lambda_E2I_L[l];
		model->AddSynapseGroup(EXCITATORY_NEURONS[l], INHIBITORY_NEURONS[l], &E2I_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS);
		if (inh_layer_on[l]){
			I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS.biological_conductance_scaling_constant_lambda = layerwise_conductance_scaling_I2E_L[l]; // layerwise_biological_conductance_scaling_constant_lambda_I2E_L[l];
			model->AddSynapseGroup(INHIBITORY_NEURONS[l], EXCITATORY_NEURONS[l], &I2E_L_INHIBITORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS);
		}
		if(E2E_L_ON)
			model->AddSynapseGroup(EXCITATORY_NEURONS[l], EXCITATORY_NEURONS[l], &E2E_L_EXCITATORY_CONDUCTANCE_SPIKING_SYNAPSE_PARAMETERS);
	}

	adding_synapses_timer.stop_timer_and_log_time_and_message("Synapses Added.", true);



	/*
	 *
	 *	Finalize and Run the model
	 *
	 */

	// Load any weights before finalising the model
	if (simulation_params.not_found() != simulation_params.find(PRELOAD_WEIGHTS)){
		std::string weight_file_path = simulation_params.get<string>(PRELOAD_WEIGHTS) + "/Synapses_NetworkWeights.bin";
		std::cout << "Loading weights from: " << weight_file_path << std::endl;
		load_weights(model, weight_file_path, true);
	}
	model->finalise_model();


	float avg_rate_stimuli = simulation_params.get<float>(NETWORK_PARAMS "." AVG_RATE_STIMULI);

	/*
	 *	INITIAL TESTING
	 */

	Simulator_Options simulator_options_initial;
	simulator_options_initial.run_simulation_general_options->presentation_time_per_stimulus_per_epoch = presentation_time_per_stimulus_per_epoch_test;
	simulator_options_initial.run_simulation_general_options->apply_plasticity_to_relevant_synapses = false;
	simulator_options_initial.stimuli_presentation_options->reset_model_state_between_each_stimulus = true;
	simulator_options_initial.recording_electrodes_options->count_neuron_spikes_recording_electrodes_bool = true;
	simulator_options_initial.recording_electrodes_options->collect_neuron_spikes_recording_electrodes_bool = true;
	simulator_options_initial.file_storage_options->save_recorded_neuron_spikes_to_file = true;
	simulator_options_initial.recording_electrodes_options->count_input_neuron_spikes_recording_electrodes_bool = true;
	simulator_options_initial.recording_electrodes_options->collect_input_neuron_spikes_recording_electrodes_bool = true;
	simulator_options_initial.file_storage_options->save_recorded_input_neuron_spikes_to_file = true;
	simulator_options_initial.file_storage_options->write_initial_synaptic_weights_to_file_bool = true;
	simulator_options_initial.file_storage_options->output_directory = output_folder_path + experimentName + "/" + "initial/";
	simulator_options_initial.recording_electrodes_options->network_state_archive_recording_electrodes_bool = true;
	simulator_options_initial.recording_electrodes_options->network_state_archive_optional_parameters->human_readable_storage = human_readable_storage;

	// float presentation_time_per_stimulus_per_epoch = presentation_time_per_stimulus_per_epoch_test;
	// bool record_spikes = record_spikes_test;
	// bool save_recorded_spikes_and_states_to_file = save_recorded_spikes_and_states_to_file_test;

	// Load the desired input stimuli and equalize their rate
	input_neurons->set_up_rates(simulation_params.get<string>(EXPERIMENT_SPECS "." TESTING_STIMULI_LIST).c_str(), "FilterParameters.txt", filepath_stimuli.c_str(), max_FR_of_input_Gabor);
	equalize_rates(input_neurons, avg_rate_stimuli);
	input_neurons->copy_rates_to_device();

	// Run the untrained initial network
	Simulator* simulator_initial = new Simulator(model, &simulator_options_initial);
	simulator_initial->RunSimulation();
	delete simulator_initial;
	/*
	* END INITIAL TESTING
	*/


	/*
	*
	* TRAIN THE NETWORK AND TEST IN BETWEEN
	*/
	if(!simulation_params.get_optional<string>(ONLY_TEST_STIMULI_FOLDER)){
		//but only if we are not in the 'test only mode'

		// Loop through the number of epochs
		int number_of_epochs_train = simulation_params.get<int>(EPOCH_COUNT);
		int test_network_every_n_epochs = simulation_params.get<int>(TEST_EVERY_N_EPOCHS);
		int stop_lr_increase_epoch = simulation_params.get<int>(EXPERIMENT_SPECS "." LEARNING_RATE_INC_STOP);
		int start_epoch = simulation_params.get<int>(START_EPOCH);

		Simulator::CreateDirectoryForSimulationDataFiles(experimentName + "/training", output_folder_path);
		Simulator::CreateDirectoryForSimulationDataFiles(experimentName + "/testing", output_folder_path);

		for (int g = start_epoch; g <= start_epoch + number_of_epochs_train; g++){

			/*
			* TRAIN NETWORK
			*/

			// Load the desired input stimuli and equalize their rate
			input_neurons->set_up_rates(simulation_params.get<string>(EXPERIMENT_SPECS "." TRAINING_STIMULI_LIST).c_str(), "FilterParameters.txt", filepath_stimuli.c_str(), max_FR_of_input_Gabor);
			equalize_rates(input_neurons, avg_rate_stimuli);
			input_neurons->copy_rates_to_device();

			//increasing learning rate
			float lr_tmp = learning_rate_rho;
			if (g <= stop_lr_increase_epoch){
				lr_tmp = powf(1.1f, (g-1)) * learning_rate_rho;
			}
			else {
				lr_tmp = powf(1.1f, (stop_lr_increase_epoch)) * learning_rate_rho;
				std::cout << "fixed learning rate" << std::endl;
			}

			std::cout << "Learning rate at Epoch " << g << " : " << lr_tmp << std::endl;

			STDP_PARAMS.learning_rate_rho = lr_tmp;

			//UGLY #TODO
			bool save_spikes_while_training = false;
			auto optional_save_spikes_while_training = simulation_params.get_optional<int>(RECORD_SPIKES_IN_TRAINING);
			if(optional_save_spikes_while_training && *optional_save_spikes_while_training==1) save_spikes_while_training=true;

			// Now create a set of options for this training epoch
			Simulator_Options simulator_options_train;
			simulator_options_train.run_simulation_general_options->apply_plasticity_to_relevant_synapses = true;
			simulator_options_train.run_simulation_general_options->number_of_epochs = 1;
			simulator_options_train.run_simulation_general_options->presentation_time_per_stimulus_per_epoch = presentation_time_per_stimulus_per_epoch_train;
			simulator_options_train.run_simulation_general_options->stimulus_presentation_order_seed = g;

			// Set the output folder for this epoch
			simulator_options_train.file_storage_options->output_directory = output_folder_path + experimentName + "/training/epoch" + string(to_string(g)) + "/";
			simulator_options_train.recording_electrodes_options->network_state_archive_optional_parameters->human_readable_storage = human_readable_storage;

			simulator_options_train.recording_electrodes_options->count_neuron_spikes_recording_electrodes_bool = save_spikes_while_training;
			simulator_options_train.recording_electrodes_options->collect_neuron_spikes_recording_electrodes_bool = save_spikes_while_training;
			simulator_options_train.file_storage_options->save_recorded_neuron_spikes_to_file = save_spikes_while_training;

			simulator_options_train.stimuli_presentation_options->presentation_format = PRESENTATION_FORMAT_OBJECT_BY_OBJECT_RESET_BETWEEN_OBJECTS;
			simulator_options_train.stimuli_presentation_options->object_order = OBJECT_ORDER_RANDOM;
			simulator_options_train.stimuli_presentation_options->transform_order = TRANSFORM_ORDER_RANDOM;

			// Now finally, create the Simulator and run this epoch of training
			Simulator* simulator_train = new Simulator(model, &simulator_options_train);
			simulator_train->RunSimulation();
			delete simulator_train;



			if (g % test_network_every_n_epochs == 0 || g == number_of_epochs_train) {

				// Load the desired input stimuli and equalize their rate
				input_neurons->set_up_rates(simulation_params.get<string>(EXPERIMENT_SPECS "." TESTING_STIMULI_LIST).c_str(), "FilterParameters.txt", filepath_stimuli.c_str(), max_FR_of_input_Gabor);
				equalize_rates(input_neurons, avg_rate_stimuli);
				input_neurons->copy_rates_to_device();

				// Given that the weights shall still be those loaded when training, go ahead and carry out test
				// Set up the options for the simulator as before
				Simulator_Options simulator_options_test;
				simulator_options_test.run_simulation_general_options->number_of_epochs = 1;
				simulator_options_test.run_simulation_general_options->apply_plasticity_to_relevant_synapses = false;
				simulator_options_test.stimuli_presentation_options->reset_model_state_between_each_stimulus = true;
				simulator_options_test.run_simulation_general_options->presentation_time_per_stimulus_per_epoch = presentation_time_per_stimulus_per_epoch_test;
				simulator_options_test.recording_electrodes_options->count_neuron_spikes_recording_electrodes_bool = true;
				simulator_options_test.recording_electrodes_options->collect_neuron_spikes_recording_electrodes_bool = true;
				simulator_options_test.recording_electrodes_options->network_state_archive_recording_electrodes_bool = true; // write the network architecture to file, AND THE SO FAR LEARNED WEIGHTS
				simulator_options_test.file_storage_options->save_recorded_neuron_spikes_to_file = true;
				//simulator_options_test.file_storage_options->write_initial_synaptic_weights_to_file_bool = true;


				simulator_options_test.file_storage_options->output_directory = output_folder_path + experimentName + "/testing/epoch" + string(to_string(g)) + "/" ;
				// The network state should already be recorded in the training folder
				//simulator_options_test.recording_electrodes_options->network_state_archive_recording_electrodes_bool = true;
				// simulator_options_train.recording_electrodes_options->network_state_archive_optional_parameters->human_readable_storage = true;

				// Finally create the Simulator and run this epoch of testing
				Simulator* simulator_test = new Simulator(model, &simulator_options_test);
				simulator_test->RunSimulation();
				delete simulator_test;
			}
		}
		// funal Code
		cout << "before system call" << endl;
		system(("/Users/clemens/.virtualenvs/myscipy/bin/python /Users/clemens/Documents/Code/AnalysisToolbox/spikeAnalysisToolsV2/quick_overview.py -n " + experimentName + " &").c_str());
		cout << "after system call" << endl;
	}

	return 0;
}
