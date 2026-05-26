#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <cuda_runtime.h>

__global__ void compute_zscore (
				double *d_in,
				double *d_total_mu,
				double *d_total_var,
				double *d_total_zscore,
				int samples_per_channel,
				int window_size)
{
	// 1. Allocating the thread ID
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	// 2. Setting up the constant variables
	const double alpha = 0.1;
	const double epsilon = 0.000001;

	// 3. Defining the starting and ending points for each window
	int begin_index = (tid - window_size + 1 > 0) ? (tid - window_size + 1) : 0;
    int end_index = tid;
    int curr_len = end_index - begin_index + 1;

    // 4. Initializing the per-channel variables
    double mu = d_in[begin_index];
    double var = 0.0;
    double z_score = 0.0;

    // 5. Performing the Z-score, mean and variance
    for (int i = 1; i < curr_len; ++i){
    	double frame = d_in[begin_index + i];

    	mu = alpha * frame + (1.0 - alpha) * mu;
    	var = alpha * pow(frame - mu, 2) + (1.0 - alpha) * var;
    	z_score = (frame - mu) / sqrt(var + epsilon);
    }

    // 6. Assigning the current results to respective thread ID
    d_total_mu[tid] = mu;
    d_total_var[tid] = var;
    d_total_zscore[tid] = z_score;
}

int main (int argc, char *argv[]) {
	
	// 1. Setting up channel metadata
	const int total_channels = 250;
	const int samples_per_channel = 5000000;
	const int window_size = 1024;
	const int num_GPUs = 8;
	size_t total_samples = (size_t)total_channels * samples_per_channel;
	omp_set_num_threads(num_GPUs);


	// 2. Host memory allocation
	double *h_in = (double*)malloc(total_samples * sizeof(double));
	double *h_total_mu = (double *)malloc(total_samples * sizeof(double));
    double *h_total_var = (double *)malloc(total_samples * sizeof(double));
    double *h_total_zscore = (double *)malloc(total_samples * sizeof(double));

    // 3. Reading neural input binary file
    const char* filename = "neural_data.bin";
	FILE* file = fopen(filename, "rb");
	size_t samples_read = fread(h_in, sizeof(double), total_samples, file);
	fclose(file);

	FILE *out_file = fopen("output_file.txt","w"); // output file

	// 4. GPU Threads incrementing
	for (int exp = 0; exp <= 10; exp++){
		int threads_per_block = 1 << exp;
		int blocks_per_grid = (samples_per_channel + threads_per_block - 1) / threads_per_block;

		// 5. Initialize dummy values for output variables
		for(size_t i = 0; i < total_samples; i++) {
			h_total_mu[i] = 0.0;
			h_total_var[i] = 0.0;
			h_total_zscore[i] = 0.0;
		}
	
		double start_timer = omp_get_wtime();

		// 6. Spawn OpenMP CPU threads
		#pragma omp parallel
		{
			// 7. Get current thread ID
			int tid = omp_get_thread_num(); // Get thread id
			cudaSetDevice(tid);

			// 8. Divide the total channels equally to number of available GPUs
			int channels_per_gpu = total_channels / num_GPUs; // 500 / 8 = 62
			int channels_remain = total_channels % num_GPUs; // 500 % 8 = 4

			// 9. Edge cases for if there is remainig channels, assign to first fews
			int begin_channel_index = tid * channels_per_gpu + (tid < channels_remain ? tid : channels_remain);
			int channels_remain_assign = channels_per_gpu + (tid < channels_remain ? 1 : 0); // Assign remaining channels to first few GPUs
			int end_channel_index = begin_channel_index + channels_remain_assign;

			// 10. Allocate VRAM per GPU device memory
			double *d_in, *d_total_mu, *d_total_var, *d_total_zscore;
			cudaMalloc((void**)&d_in, samples_per_channel * sizeof(double));
			cudaMalloc((void**)&d_total_mu, samples_per_channel * sizeof(double));
			cudaMalloc((void**)&d_total_var, samples_per_channel * sizeof(double));
			cudaMalloc((void**)&d_total_zscore, samples_per_channel * sizeof(double));

			// 11. Process each channel one by one
			for (int i = begin_channel_index; i < end_channel_index; ++i){
				size_t offset = (size_t)i * samples_per_channel;

				// 12. Copy input from CPU RAM to GPU VRAM
				cudaMemcpy(d_in, &h_in[offset], samples_per_channel * sizeof(double), cudaMemcpyHostToDevice);

				// 13. Trigger the kernel
				compute_zscore <<<blocks_per_grid, threads_per_block>>>(
					d_in,
					d_total_mu,
					d_total_var,
					d_total_zscore,
					samples_per_channel,
					window_size
				);

				cudaDeviceSynchronize();

				// 14. Copy output from GPU VRAM to CPU RAM
				cudaMemcpy(&h_total_mu[offset], d_total_mu, samples_per_channel * sizeof(double), cudaMemcpyDeviceToHost);
				cudaMemcpy(&h_total_var[offset], d_total_var, samples_per_channel * sizeof(double), cudaMemcpyDeviceToHost);
				cudaMemcpy(&h_total_zscore[offset], d_total_zscore, samples_per_channel * sizeof(double), cudaMemcpyDeviceToHost);
			
			}
			printf("Running from GPU %d\n", tid);

			// 15. Free resources (VRAM)
			cudaFree(d_in);
			cudaFree(d_total_mu);
			cudaFree(d_total_var);
			cudaFree(d_total_zscore);
		}
		double end_timer = omp_get_wtime();
		
		/* Debugging mathematical correctness
		for(int c = 0; c < total_channels; ++c){
			for(int s = 0; s < samples_per_channel; ++s){
				size_t index = (size_t)c * samples_per_channel + s; // 1D array offset

				fprintf(out_file, "%.6lf ", h_total_zscore[index]);
			}
			fprintf(out_file, "\n"); // New line for every channel
		}
		*/

		fprintf(out_file, "%d %d %d %.6lf\n", exp, threads_per_block, num_GPUs, (end_timer - start_timer));
	}

	fclose(out_file);
	
	// 16. Free resources (Host RAM)
	free(h_in);
	free(h_total_mu);
	free(h_total_var);
	free(h_total_zscore);
	return 0;
}