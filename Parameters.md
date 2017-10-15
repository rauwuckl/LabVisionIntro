# Parameter File

The example model contains a number of parameters. An attempt has been made to name these in a sensible way. Nonetheless they may need some description, I have outlined each model and the descriptions below. These are split according to whether they are Neuron, Synapse, or Plasticity related.


## Neurons
	resting_potential_v0: The voltage to which the LIF neuron decays by default
	threshold_for_action_potential: The voltage at which (if the membrane potential hits) the neuron is said to reach a threshold and spike
	somatic_capacitance_Cm: The capacitance of the membrane of the cell (with the leakage conductance can be used to calculate the membrane time constant)
	somatice_leakage_conductance_g0: the conductance of the membrane (the inverse of the membrane resistnace, i.e. R_membrane = (1 / somatic_leakage_conductance))	
	absolute_refractory_period: The length of time for which a neuron is not allowed to spike after it carries out an action potential



## Synapses

Key (XXX):
	G2E = Gabor to Excitatory Connections
	E2E = Excitatory to Excitatory connections
	E2I = Excitatory to Inhibitory connections
	I2E = Inhibitory to Excitatory connections
The above connection types can also be of 3 different types (YY):
	FF = FeedForward
	FB = FeedBack
	L = Lateral


The model has a number of parameters for different types of connections:
	E2E_FB_ON: boolean flag indicating if feedback connections in the network are on
	E2E_L_ON: boolean flag for if the network has lateral connections or not
	E2E_L_STDP_ON: boolean flag to turn on/off plasticity on the lateral connections (if they exist)
	
	inh_layer_on: list of booleans indicating if the inhibitory neuron population is created for each layer of the network (useful for optimising)
	max_number_of_connections_per_pair: By making multiple pre->post connections, 

	The connectivity in this network is such that there is a (gaussian) radius of connectivity probability from the layer before. This radius/fan-in style connectivity has a number of parameters:
		gaussian_synapses_standard_deviation_XXX_YY: The standard deviation of the gaussian probability profile that is applied to the previous layer to figure out which connections are most likely.
		fanInCount_XXX_YY: This is the number of connections that are sampled from the gaussian probability distribution (based on the radius above).

	The synaptic transmission delay is how much of a delay there is between when a neuron fires and when that spike reaches the synapse. The delays selected from a uniform distribution between the minimum and maximum transmission delay values.
		min_delay: This is the smallest value that the synaptic transmission delay can be (seconds)
		max_delay: this is the maximum value that the synaptic transmission delay can be (seconds)

	The synaptic weight is multiplied by the biological scaling factor (see below) before it eventually represents how much the conductance of a synapse increases by when a pre-synaptic spike reaches the synapse. Again, weight values are randomly (uniformly) sampled between the minimum and maximum value
		weight_range_bottom = the minimum weight value
		weight_range_top = the maximum weight value
		
		biological_conductance_scaling_constant_lambda_XXX_YY: (BSC) this is a multiplier for the weight so that it can be brought into a reasonable range. As mentioned in the intro lecture, if the time constants of the synapses and neurons are roughly equal, this value should be a fraction of the somatic leakage conductance. If it were equal to the somatic leakage conductance, a single pre-synaptic spike would cause a single post-synaptic spike. Another note: this model uses a very large synaptic time constant (150ms) compared to the exc neuron time constant (20ms), therefore this BSC should be under a tenth of the leakage conductance.
	
	The synaptic time constant determines the decay of the synaptic conductance curve. The larger this is, the longer it takes for the synaptic conductance to decay to zero.
		decay_term_tau_g_XXX_YY: The synaptic decay constant time in seconds


## Plasticity
	learning_rate_rho: this values is multiplied by any weight changes calculated by the STDP rule before applying that change to the synapse. This essentially speeds up/slows down learning.
	decay_term_tau_C: the decay term for the LTP side of the STDP curve. See the Ben Evans papers to understand the learning rule
	decay_term_tay_D: the decay term for the LTD side of the STDP curve. As above
	alpha_C: the height of the LTP side of the STDP curve
	alpha_D: the height of the LTD side of the STDP curve

	
