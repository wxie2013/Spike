#include "IzhikevichSpikingNeurons.h"
#include <stdlib.h>
#include <stdio.h>
#include "CUDAErrorCheckHelpers.h"


// IzhikevichSpikingNeurons Constructor
IzhikevichSpikingNeurons::IzhikevichSpikingNeurons() {

}


// IzhikevichSpikingNeurons Destructor
IzhikevichSpikingNeurons::~IzhikevichSpikingNeurons() {

}


int IzhikevichSpikingNeurons::AddGroupNew(neuron_struct *params, int group_shape[2]){

	int new_group_id = SpikingNeurons::AddGroupNew(params, group_shape);

	param_a = (float*)realloc(param_a, (total_number_of_neurons*sizeof(float)));
	param_b = (float*)realloc(param_b, (total_number_of_neurons*sizeof(float)));

	return new_group_id;
}


void IzhikevichSpikingNeurons::initialise_device_pointersNew() {
 	
 	SpikingNeurons::initialise_device_pointersNew();


 	CudaSafeCall(cudaMalloc((void **)&d_param_a, sizeof(float)*total_number_of_neurons));
 	CudaSafeCall(cudaMalloc((void **)&d_param_b, sizeof(float)*total_number_of_neurons));
	
	reset_neuron_variables_and_spikesNew();
}

void IzhikevichSpikingNeurons::reset_neuron_variables_and_spikesNew() {

	SpikingNeurons::reset_neuron_variables_and_spikesNew();	

	CudaSafeCall(cudaMemcpy(d_param_a, param_a, sizeof(float)*total_number_of_neurons, cudaMemcpyHostToDevice));
	CudaSafeCall(cudaMemcpy(d_param_b, param_b, sizeof(float)*total_number_of_neurons, cudaMemcpyHostToDevice));

}


__global__ void izhikevich_state_update(float *d_states_v,
								float *d_states_u,
								float *d_param_a,
								float *d_param_b,
								float* currentinj,
								float timestep,
								size_t total_number_of_neurons);


void IzhikevichSpikingNeurons::izhikevich_state_update_wrapper(float timestep) {

	izhikevich_state_update<<<number_of_neuron_blocks_per_grid, threads_per_block>>>(d_states_v,
																	d_states_u,
																	d_param_a,
																	d_param_b,
																	d_current_injections,
																	timestep,
																	total_number_of_neurons);

	CudaCheckError();
}


// State Update
__global__ void izhikevich_state_update(float *d_states_v,
								float *d_states_u,
								float *d_param_a,
								float *d_param_b,
								float* currentinj,
								float timestep,
								size_t total_number_of_neurons){

	// We require the equation timestep in ms:
	float eqtimestep = timestep*1000.0f;
	// Get thread IDs
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < total_number_of_neurons) {
		// Update the neuron states according to the Izhikevich equations
		float v_update = 0.04f*d_states_v[idx]*d_states_v[idx] 
							+ 5.0f*d_states_v[idx]
							+ 140 
							- d_states_u[idx]
							+ currentinj[idx];

		d_states_v[idx] += eqtimestep*v_update;
		d_states_u[idx] += eqtimestep*(d_param_a[idx] * (d_param_b[idx] * d_states_v[idx] - 
							d_states_u[idx]));
	}
	__syncthreads();
}


