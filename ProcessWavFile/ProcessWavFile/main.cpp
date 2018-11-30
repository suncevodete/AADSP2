
#include <stdlib.h>
#include <string.h>
#include "WAVheader.h"

#define BLOCK_SIZE 16
#define MAX_NUM_CHANNEL 4
#define INPUT_GAIN  0.8912509381337456
#define MODE 0 //0 = 3_2_0; 1 = 2_0_0

double sampleBuffer[MAX_NUM_CHANNEL+1][BLOCK_SIZE];
double helpBuffer[MAX_NUM_CHANNEL+1][BLOCK_SIZE];
double gain[] = {0.6309573444801932, 0.31622776601683794, 0.6382634861905486, 0.33496543915782767}; //-4, -10. -3.9, -9.5 db

double LPF_11kHz_coeff[6] = {0.29185073257004568000, 0.58370146514009136000, 0.29185073257004568000,
							1.00000000000000000000, -0.00417302338598968350, 0.17157595366680914000};
double HPF_3kHz_coeff[6] = {0.73855434570188305000, -1.47710869140376610000, 0.73855434570188305000,
							1.00000000000000000000, -1.40750534395471780000, 0.54664949370997029000};
double HPF_5kHz_coeff[6] = {0.60074547832695790000, -1.20149095665391580000, 0.60074547832695790000,
							1.00000000000000000000, -1.03517120979351820000, 0.36781068948958456000};

double x_history0[2];
double y_history0[2];
double x_history1[2];
double y_history1[2];
double x_history2[2];
double y_history2[2];
double x_history3[2];
double y_history3[2];
double x_history4[2];
double y_history4[2];
double x_history5[2];
double y_history5[2];

double second_order_IIR(double input, double* coefficients, double* x_history, double* y_history) 
{
    double output = 0;

    output += coefficients[0] * input;        /* A0 * x(n)     */
    output += coefficients[1] * x_history[0]; /* A1 * x(n-1) */
    output += coefficients[2] * x_history[1]; /* A2 * x(n-2)   */
    output -= coefficients[4] * y_history[0]; /* B1 * y(n-1) */
    output -= coefficients[5] * y_history[1]; /* B2 * y(n-2)   */

    y_history[1] = y_history[0];    /* y(n-2) = y(n-1) */
    y_history[0] = output; /* y(n-1) = y(n)   */
    x_history [1] = x_history [0];  /* x(n-2) = x(n-1) */
    x_history [0] = input;          /* x(n-1) = x(n)   */

    return output;
}

void processing(double in[][BLOCK_SIZE], double out[][BLOCK_SIZE])
{
	double left[BLOCK_SIZE]; //ovde cuvam levi i desni kanal da bih sampleBuffer mogla da menjam
	double right[BLOCK_SIZE];

	/* mnozim sa pocetnim gainom */
	for (int m = 0; m < BLOCK_SIZE; m++) 
	{
		left[m] = in[0][m] * INPUT_GAIN;
		right[m] = in[1][m] * INPUT_GAIN;
	}
	//pomocni baferi koji ce predstavljati izlazne kanale
	double tmp_0[BLOCK_SIZE];
	double tmp_1[BLOCK_SIZE];
	double tmp_2[BLOCK_SIZE];
	double tmp_3[BLOCK_SIZE];
	double tmp_4[BLOCK_SIZE];
	double tmp_5[BLOCK_SIZE];
	
	for (int i = 0; i < BLOCK_SIZE; i++)
	{
		tmp_0[i] = second_order_IIR(left[i], LPF_11kHz_coeff, x_history0, y_history0); //ovaj ide direktno na izlaz
		tmp_1[i] = second_order_IIR(tmp_0[i], HPF_3kHz_coeff, x_history1, y_history1);
		tmp_1[i] = tmp_1[i] * gain[0]; //zavrsen out[0]
		tmp_2[i] = second_order_IIR(left[i], HPF_5kHz_coeff, x_history2, y_history2);
		tmp_2[i] = tmp_2[i] * gain[1]; //ostaje da ga saberes sa tmp_2
		//tmp_1[i] = tmp_1[i] + tmp_2[i];
		tmp_3[i] = second_order_IIR(right[i], LPF_11kHz_coeff, x_history3, y_history3); //zavrsen out[3]
		tmp_4[i] = second_order_IIR(right[i], HPF_5kHz_coeff, x_history4, y_history4); 
		tmp_4[i] = tmp_4[i] * gain[3];
		tmp_5[i] = second_order_IIR(tmp_3[i], HPF_3kHz_coeff, x_history5, y_history5);
		tmp_5[i] = tmp_5[i] * gain[2];
	}

	if (!MODE) {
		for(int j = 0; j < BLOCK_SIZE; j++)
		{
		out[0][j] = tmp_0[j];
		out[1][j] = tmp_1[j] + tmp_2[j];
		out[2][j] = tmp_2[j];
		out[3][j] = tmp_3[j];
		out[4][j] = tmp_4[j] + tmp_5[j];
		}
	} else {
		for(int j = 0; j < BLOCK_SIZE; j++)
		{ //samo 2 kanala, one gore racunas jer su ti ulaz u tmp1 i tmp5
		out[0][j] = tmp_0[j];
		out[1][j] = tmp_1[j] + tmp_2[j];
		out[2][j] = tmp_2[j];
		out[3][j] = tmp_3[j];
		out[4][j] = tmp_4[j] + tmp_5[j];
		}
	}

	
	
}

int main(int argc, char* argv[])
{
	FILE *wav_in=NULL;
	FILE *wav_out=NULL;
	char WavInputName[256];
	char WavOutputName[256];
	WAV_HEADER inputWAVhdr,outputWAVhdr;	

	// Init channel buffers
	for(int i=0; i<MAX_NUM_CHANNEL; i++)
	{
		memset(&sampleBuffer[i], 0, BLOCK_SIZE);
		memset(&helpBuffer[i], 0, BLOCK_SIZE);
	}
	// Open input and output wav files
	//-------------------------------------------------
	strcpy(WavInputName,argv[1]);
	wav_in = OpenWavFileForRead (WavInputName,"rb");
	strcpy(WavOutputName,argv[2]);
	wav_out = OpenWavFileForRead (WavOutputName,"wb");
	//-------------------------------------------------

	// Read input wav header
	//-------------------------------------------------
	ReadWavHeader(wav_in,inputWAVhdr);
	//-------------------------------------------------
	
	// Set up output WAV header
	//-------------------------------------------------	
	outputWAVhdr = inputWAVhdr;
	outputWAVhdr.fmt.NumChannels = MAX_NUM_CHANNEL + 1; // change number of channels

	int oneChannelSubChunk2Size = inputWAVhdr.data.SubChunk2Size/inputWAVhdr.fmt.NumChannels;
	int oneChannelByteRate = inputWAVhdr.fmt.ByteRate/inputWAVhdr.fmt.NumChannels;
	int oneChannelBlockAlign = inputWAVhdr.fmt.BlockAlign/inputWAVhdr.fmt.NumChannels;
	
	outputWAVhdr.data.SubChunk2Size = oneChannelSubChunk2Size*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.ByteRate = oneChannelByteRate*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.BlockAlign = oneChannelBlockAlign*outputWAVhdr.fmt.NumChannels;


	// Write output WAV header to file
	//-------------------------------------------------
	WriteWavHeader(wav_out,outputWAVhdr);

	// Processing loop
	//-------------------------------------------------	
	{
		int sample;
		int BytesPerSample = inputWAVhdr.fmt.BitsPerSample/8;
		const double SAMPLE_SCALE = -(double)(1 << 31);		//2^31
		int iNumSamples = inputWAVhdr.data.SubChunk2Size/(inputWAVhdr.fmt.NumChannels*inputWAVhdr.fmt.BitsPerSample/8);
		
		// exact file length should be handled correctly...
		for(int i=0; i<iNumSamples/BLOCK_SIZE; i++)
		{	
			for(int j=0; j<BLOCK_SIZE; j++)
			{
				for(int k=0; k<inputWAVhdr.fmt.NumChannels; k++)
				{	
					sample = 0; //debug
					fread(&sample, BytesPerSample, 1, wav_in);
					sample = sample << (32 - inputWAVhdr.fmt.BitsPerSample); // force signextend
					sampleBuffer[k][j] = sample / SAMPLE_SCALE;				// scale sample to 1.0/-1.0 range		
				}
			}
			
			
			processing(sampleBuffer, sampleBuffer);

			for(int j=0; j<BLOCK_SIZE; j++)
			{
				for(int k=0; k<outputWAVhdr.fmt.NumChannels; k++)
				{	
					sample = sampleBuffer[k][j] * SAMPLE_SCALE ;	// crude, non-rounding 			
					sample = sample >> (32 - inputWAVhdr.fmt.BitsPerSample);
					fwrite(&sample, outputWAVhdr.fmt.BitsPerSample/8, 1, wav_out);		
				}
			}		
		}
	}
	
	// Close files
	//-------------------------------------------------	
	fclose(wav_in);
	fclose(wav_out);
	//-------------------------------------------------	

	return 0;
}