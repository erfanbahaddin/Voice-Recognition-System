#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "portaudio.h"
#include "libmfcc.h"



#define SAMPLE_RATE  (44100)
#define FRAMES_PER_BUFFER (512)
#define NUM_SECONDS     (3)
#define NUM_CHANNELS    (2)
#define DITHER_FLAG     (0)
#define VERY_BIG (1e30)


#if 1
#define PA_SAMPLE_TYPE  paFloat32
typedef float SAMPLE;
#define SAMPLE_SILENCE  (0.0f)
#define PRINTF_S_FORMAT "%.8f"
#elif 1
#define PA_SAMPLE_TYPE  paInt16
typedef short SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#elif 0
#define PA_SAMPLE_TYPE  paInt8
typedef char SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#else
#define PA_SAMPLE_TYPE  paUInt8
typedef unsigned char SAMPLE;
#define SAMPLE_SILENCE  (128)
#define PRINTF_S_FORMAT "%d"
#endif



typedef struct {
    int          frameIndex;  /* Index into sample array. */
    int          maxFrameIndex;
    SAMPLE      *recordedSamples;
}paTestData;

/*========================================Record Callback Function========================================*/

static int recordCallback( const void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo* timeInfo, PaStreamCallbackFlags statusFlags, void *userData ) {
    paTestData *data = (paTestData*)userData;
    const SAMPLE *rptr = (const SAMPLE*)inputBuffer;
    SAMPLE *wptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    long framesToCalc;
    long i;
    int finished;
    unsigned long framesLeft = data->maxFrameIndex - data->frameIndex;
    (void) outputBuffer; /* Prevent unused variable warnings. */
    (void) timeInfo;
    (void) statusFlags;
    (void) userData;
    if( framesLeft < framesPerBuffer ) {
        framesToCalc = framesLeft;
        finished = paComplete;
    }
    else {
        framesToCalc = framesPerBuffer;
        finished = paContinue;
    }
    if( inputBuffer == NULL ) {
        for( i=0; i<framesToCalc; i++ ) {
            *wptr++ = SAMPLE_SILENCE;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = SAMPLE_SILENCE;  /* right */
        }
    }
    else {
        for( i=0; i<framesToCalc; i++ ) {
            *wptr++ = *rptr++;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  /* right */
        }
    }
    data->frameIndex += framesToCalc;
    return finished;
}

/*========================================Play Callback Function========================================*/

static int playCallback( const void *inputBuffer, void *outputBuffer, unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo* timeInfo, PaStreamCallbackFlags statusFlags, void *userData ) {
    paTestData *data = (paTestData*)userData;
    SAMPLE *rptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
    SAMPLE *wptr = (SAMPLE*)outputBuffer;
    unsigned int i;
    int finished;
    unsigned int framesLeft = data->maxFrameIndex - data->frameIndex;
    (void) inputBuffer; /* Prevent unused variable warnings. */
    (void) timeInfo;
    (void) statusFlags;
    (void) userData;
    if( framesLeft < framesPerBuffer ) {
        /* final buffer... */
        for( i=0; i<framesLeft; i++ ) {
            *wptr++ = *rptr++;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  /* right */
        }
        for( ; i<framesPerBuffer; i++ ) {
            *wptr++ = 0;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = 0;  /* right */
        }
        data->frameIndex += framesLeft;
        finished = paComplete;
    }
    else {
        for( i=0; i<framesPerBuffer; i++ ) {
            *wptr++ = *rptr++;  /* left */
            if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  /* right */
        }
        data->frameIndex += framesPerBuffer;
        finished = paContinue;
    }
    return finished;
}

/*=========================Pre Emphasis Function=========================*/

void preEmphasis(float input[], float output[], int n) {
    float lastInput = 0;
    int i;
    for (i = 0; i < n; i++) {
	output[i] = input[i] - (0.95 * lastInput);
	lastInput = input[i];
    }
}

/*=========================Hamming Function=========================*/

void hamming(float input[], float output[], int n) {
    int i;
    double Pi = 3.14159;
    for (i = 0; i < n; i++) {
	output[i] = 0.54 - (0.46 * cos(2 * Pi * i / n));
    }
}

/*=========================Discrete Fourier Transform Function=========================*/

void computeDFT(float input[], float output[], int n) {
    double Pi = 3.14159;
    for (int k = 0; k < n; k++) {
	float sum = 0;
	for (int t = 0; t < n; t++) {
	    double angle = 2 * Pi * t * k / n;
	    sum +=  input[t] * cos(angle);
	}
        output[k] = sum;
    }
}

/*========================================MFCC Functions========================================*/

double GetCoefficient(double* spectralData, unsigned int samplingRate, unsigned int NumFilters, unsigned int binSize, unsigned int m) {
    double result = 0.0f;
    double outerSum = 0.0f;
    double innerSum = 0.0f;
    unsigned int k, l;
    if(m >= NumFilters) {
	return 0.0f;
    }
    result = NormalizationFactor(NumFilters, m);
    for(l = 1; l <= NumFilters; l++) {
	innerSum = 0.0f;
	for(k = 0; k < binSize - 1; k++) {
	    innerSum += fabs(spectralData[k] * GetFilterParameter(samplingRate, binSize, k, l));
	}
	if(innerSum > 0.0f) {
	    innerSum = log(innerSum); // The log of 0 is undefined, so don't use it
	}
	innerSum = innerSum * cos(((m * PI) / NumFilters) * (l - 0.5f));
	outerSum += innerSum;
    }
    result *= outerSum;
    return result;
}



double NormalizationFactor(int NumFilters, int m) {
    double normalizationFactor = 0.0f;
    if(m == 0) {
	normalizationFactor = sqrt(1.0f / NumFilters);
    }
    else {
	normalizationFactor = sqrt(2.0f / NumFilters);
    }
    return normalizationFactor;
}



double GetFilterParameter(unsigned int samplingRate, unsigned int binSize, unsigned int frequencyBand, unsigned int filterBand) {
    double filterParameter = 0.0f;
    double boundary = (frequencyBand * samplingRate) / binSize;		// k * Fs / N
    double prevCenterFrequency = GetCenterFrequency(filterBand - 1);		// fc(l - 1) etc.
    double thisCenterFrequency = GetCenterFrequency(filterBand);
    double nextCenterFrequency = GetCenterFrequency(filterBand + 1);
    if(boundary >= 0 && boundary < prevCenterFrequency) {
        filterParameter = 0.0f;
    }
    else if(boundary >= prevCenterFrequency && boundary < thisCenterFrequency) {
	filterParameter = (boundary - prevCenterFrequency) / (thisCenterFrequency - prevCenterFrequency);
	filterParameter *= GetMagnitudeFactor(filterBand);
    }
    else if(boundary >= thisCenterFrequency && boundary < nextCenterFrequency) {
	filterParameter = (boundary - nextCenterFrequency) / (thisCenterFrequency - nextCenterFrequency);
	filterParameter *= GetMagnitudeFactor(filterBand);
    }
    else if(boundary >= nextCenterFrequency && boundary < samplingRate) {
	filterParameter = 0.0f;
    }
    return filterParameter;
}



double GetMagnitudeFactor(unsigned int filterBand) {
    double magnitudeFactor = 0.0f;
    if(filterBand >= 1 && filterBand <= 14) {
	magnitudeFactor = 0.015;
    }
    else if(filterBand >= 15 && filterBand <= 48) {
	magnitudeFactor = 2.0f / (GetCenterFrequency(filterBand + 1) - GetCenterFrequency(filterBand -1));
    }
    return magnitudeFactor;
}



double GetCenterFrequency(unsigned int filterBand) {
    double centerFrequency = 0.0f;
    double exponent;
    if(filterBand == 0) {
	centerFrequency = 0;
    }
    else if(filterBand >= 1 && filterBand <= 14) {
	centerFrequency = (200.0f * filterBand) / 3.0f;
    }
    else {
	exponent = filterBand - 14.0f;
	centerFrequency = pow(1.0711703, exponent);
	centerFrequency *= 1073.4;
    }
    return centerFrequency;
}

/*========================================DTW Function========================================*/

float DTW(char *file1_name, char *file2_name , int line1 , int line2 , int columns , char *debug_name){

double **globdist;
double **Dist;

double top, mid, bot, cheapest, total;
unsigned short int **move;
unsigned short int **warp;
unsigned short int **temp;

unsigned int I, X, Y, n, i, j, k;
unsigned int xsize = line1;
unsigned int ysize = line2;
unsigned int params = columns;

unsigned int debug; /* debug flag */

float **x, **y; /*now 2 dimensional*/

FILE *glob, *debug_file, *output_file;
FILE * file_name1;
FILE * file_name2;

     /* open debug file */

     if ((debug_file = fopen(debug_name,"wb")) == NULL)
       {fprintf(stderr,"Cannot open debug file\n");
       exit(1);
       }

     debug = 1;


 /* open x-parameter file */

if ((file_name1=fopen(file1_name,"rb"))==NULL)
{fprintf(stderr,"File %s cannot be opened\n",file1_name);
exit(1);
}

/* open y-parameter file */

if ((file_name2=fopen(file2_name,"rb"))==NULL)
{fprintf(stderr,"File %s cannot be opened\n",file2_name);
exit(1);
}

if (debug==1) fprintf(debug_file,"xsize %d ysize %d params %d\n",xsize,ysize,params);



/* allocate memory for x and y matrices */

if ((x = malloc(xsize * sizeof(float *))) == NULL)
     fprintf(stderr,"Memory allocation error (x)\n");

for (i=0; i < xsize; i++)
     if ((x[i] = malloc(params * sizeof(float))) == NULL)
     fprintf(stderr,"Memory allocation error (x)\n");

if ((y = malloc(ysize * sizeof(float *))) == NULL)
     fprintf(stderr,"Memory allocation error (y)\n");

for (i=0; i < ysize; i++)
     if ((y[i] = malloc(params * sizeof(float))) == NULL)
     fprintf(stderr,"Memory allocation error (y)\n");

     /* allocate memory for Dist */

if ((Dist = malloc(xsize * sizeof(double *))) == NULL)
     fprintf(stderr,"Memory allocation error (Dist)\n");

for (i=0; i < xsize; i++)
if ((Dist[i] = malloc(ysize * sizeof(double))) == NULL)
     fprintf(stderr,"Memory allocation error (Dist)\n");

     /* allocate memory for globdist */

if ((globdist = malloc(xsize * sizeof(double *))) == NULL)
     fprintf(stderr,"Memory allocation error (globdist)\n");

for (i=0; i < xsize; i++)
if ((globdist[i] = malloc(ysize * sizeof(double))) == NULL)
     fprintf(stderr,"Memory allocation error (globdist)\n");

     /* allocate memory for move */

if ((move = malloc(xsize * sizeof(short *))) == NULL)
     fprintf(stderr,"Memory allocation error (move)\n");

for (i=0; i < xsize; i++)
if ((move[i] = malloc(ysize * sizeof(short))) == NULL)
     fprintf(stderr,"Memory allocation error (move)\n");

     /* allocate memory for temp */

if ((temp = malloc(xsize * 2 * sizeof(short *))) == NULL)
     fprintf(stderr,"Memory allocation error (temp)\n");

for (i=0; i < xsize*2; i++)
if ((temp[i] = malloc(2 * sizeof(short))) == NULL)
     fprintf(stderr,"Memory allocation error (temp)\n");

     /* allocate memory for warp */

if ((warp = malloc(xsize * 2 * sizeof(short *))) == NULL)
     fprintf(stderr,"Memory allocation error (warp)\n");

for (i=0; i < xsize*2; i++)
if ((warp[i] = malloc(2 * sizeof(short))) == NULL)
     fprintf(stderr,"Memory allocation error (warp)\n");

//fprintf(stdout,"Reading input files\n");

/*read x parameter in x[]*/

for (i=0; i < xsize; i++)
{
  for (k=0; k < params; k++)
    {
if (feof(file_name1))
  {fprintf(stderr,"Premature EOF in %s\n","file1");
  exit(1);
  }

  fscanf(file_name1,"%f ",&x[i][k]);

if (debug == 1)
  fprintf(debug_file,"float_x[%d %d] = %f\n",i,k,x[i][k]);
    }
}

/*read y parameter in y[]*/

for (i=0; i < ysize; i++)
{
  for (k=0; k < params; k++)
    {

if (feof(file_name2))
  {fprintf(stderr,"Premature EOF in %s\n","file2");
  exit(1);
  }

fscanf(file_name2,"%f ",&y[i][k]);
 if (debug == 1)
   fprintf(debug_file,"float_y[%d %d] = %f\n",i,k,y[i][k]);
    }
}

//fprintf(stdout,"Computing distance matrix ...\n");

/*Compute distance matrix*/

for(i=0;i<xsize;i++) {
  for(j=0;j<ysize;j++) {
    total = 0;
    for (k=0;k<params;k++) {
      total = total + ((x[i][k] - y[j][k]) * (x[i][k] - y[j][k]));
    }

    Dist[i][j] = total;
if (debug == 1)
      fprintf(debug_file,"Dist: %d %d %.0f %.0f\n",i,j,total,Dist[i][j]);
  }
}


//fprintf(stdout,"Warping in progress ...\n");

/*% for first frame, only possible match is at (0,0)*/

globdist[0][0] = Dist[0][0];
for (j=1; j<xsize; j++)
  globdist[j][0] = VERY_BIG;

globdist[0][1] = VERY_BIG;
globdist[1][1] = globdist[0][0] + Dist[1][1];
move[1][1] = 2;

for(j=2;j<xsize;j++)
  globdist[j][1] = VERY_BIG;

for(i=2;i<ysize;i++) {
  globdist[0][i] = VERY_BIG;
  globdist[1][i] = globdist[0][i-1] + Dist[1][i];
  if (debug = 1)
    fprintf(debug_file,"globdist[2][%d] = %.2e\n",i,globdist[2][i]);

  for(j=2;j<xsize;j++) {
    top = globdist[j-1][i-2] + Dist[j][i-1] + Dist[j][i];
    mid = globdist[j-1][i-1] + Dist[j][i];
    bot = globdist[j-2][i-1] + Dist[j-1][i] + Dist[j][i];
    if( (top < mid) && (top < bot))
    {
    cheapest = top;
    I = 1;
    }
  else if (mid < bot)
    {
    cheapest = mid;
    I = 2;
    }
  else {cheapest = bot;
    I = 3;
    }

/*if all costs are equal, pick middle path*/
       if( ( top == mid) && (mid == bot))
   I = 2;

  globdist[j][i] = cheapest;
  move[j][i] = I;
  if (debug == 1) {
    fprintf(debug_file,"globdist[%d][%d] = %.2e\n",j,i,globdist[j][i]);
    fprintf(debug_file,"move j:%d:i:%d=%d\n",j,i,move[j][i]);
        }
      }
}

if (debug == 1) {
  for (j=0; j<xsize; j++) {
    for (i=0; i<ysize; i++) {
      fprintf(debug_file,"[%d %d] globdist: %.2e    move: %d    \n",j,i,globdist[j][i],move[j][i]);
    }
  }
}



X = ysize-1; Y = xsize-1; n=0;
warp[n][0] = X; warp[n][1] = Y;


while (X > 0 && Y > 0) {
n=n+1;


if (n>ysize *2) {fprintf (stderr,"Warning: warp matrix too large!");
exit(1);
}

 if (debug == 1)
   fprintf(debug_file,"Move %d %d %d\n", Y, X, move[Y][X]);

if (move[Y] [X] == 1 )
  {
  warp[n][0] = X-1; warp[n][1] = Y;
  n=n+1;
  X=X-2; Y = Y-1;
  }
else if (move[Y] [X] == 2)
  {
  X=X-1; Y = Y-1;
  }
else if (move[Y] [X] == 3 )
  {
  warp[n] [0] = X;
  warp[n] [1] = Y-1;
  n=n+1;
  X=X-1; Y = Y-2;
      }
else {fprintf(stderr,"Error: move not defined for X = %d Y = %d\n",X,Y);
}
warp[n] [0] =X;
warp[n] [1] =Y;

}


/*flip warp*/
for (i=0;i<=n;i++) {
  temp[i][0] = warp[n-i][0];
  temp[i][1] = warp[n-i][1];

}

for (i=0;i<=n;i++) {
  warp[i][0] = temp[i][0];
  warp[i][1] = temp[i][1];

}

// writing step (khodemun mizanimesh)

//fprintf(stdout,"Warping complete. Writing output file.\n");

/* open output file */
if ((output_file=fopen("details.txt","wb"))==NULL)
{fprintf(stderr,"File %s cannot be opened\n","details.txt");
 exit(1);
}

/*print warped trajectory to stdout*/
for (i=0;i<=n;i++)
     fprintf(output_file,"%d %d\n",warp[i][0]+1,warp[i][1]+1);

     fclose(output_file);

/* print global distance to globfile*/

if ((glob=fopen("glob","w"))==NULL)
     fprintf(stderr,"Cannot open file glob\n");

float result=globdist[xsize-1][ysize-1];
fprintf(glob,"%f\n",result);
fclose(glob);

fprintf(stdout,"DTW result : %.3f\n",globdist[xsize-1][ysize-1]);


return result;
}

/*##################################################     Main     ##################################################*/

int Record(char filename[]) {
	    PaStreamParameters  inputParameters, outputParameters;
	    PaStream*           stream;
	    PaError             err = paNoError;
	    paTestData          data;
	    int                 i;
	    int                 totalFrames;
	    int                 numSamples;
	    int                 numBytes;
	    SAMPLE              max, val;
	    double              average;
	    printf("patest_record.c\n"); fflush(stdout);
	    data.maxFrameIndex = totalFrames = NUM_SECONDS * SAMPLE_RATE; /* Record for a few seconds. */
	    data.frameIndex = 0;
	    numSamples = totalFrames * NUM_CHANNELS;
	    numBytes = numSamples * sizeof(SAMPLE);
	    data.recordedSamples = (SAMPLE *)malloc( numBytes ); /* From now on, recordedSamples is initialised. */
	    if( data.recordedSamples == NULL ) {
		printf("Could not allocate record array.\n");
		goto done;
	    }
	    for( i=0; i<numSamples; i++ ) data.recordedSamples[i] = 0;
	    err = Pa_Initialize();
	    if( err != paNoError ) goto done;
	    inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
	    if (inputParameters.device == paNoDevice) {
		fprintf(stderr,"Error: No default input device.\n");
		goto done;
	    }
	    inputParameters.channelCount = 2;                    /* stereo input */
	    inputParameters.sampleFormat = PA_SAMPLE_TYPE;
	    inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency;
	    inputParameters.hostApiSpecificStreamInfo = NULL;

/*========================================Record Audio========================================*/

	    err = Pa_OpenStream(&stream, &inputParameters, NULL, SAMPLE_RATE, FRAMES_PER_BUFFER, paClipOff, recordCallback, &data);
	    if( err != paNoError ) goto done;
	    err = Pa_StartStream( stream );
	    if( err != paNoError ) goto done;
	    printf("\n=== Now recording!! Please speak into the microphone ===\n"); fflush(stdout);
	    while( ( err = Pa_IsStreamActive( stream ) ) == 1 ) {
		Pa_Sleep(1000);
		printf("index = %d\n", data.frameIndex ); fflush(stdout);
	    }
	    if( err < 0 ) goto done;
	    err = Pa_CloseStream( stream );
	    if( err != paNoError ) goto done;

/*========================================Play Recorded Audio========================================*/

	    data.frameIndex = 0;
	    outputParameters.device = Pa_GetDefaultOutputDevice(); /* default output device */
	    if (outputParameters.device == paNoDevice) {
		fprintf(stderr,"Error: No default output device.\n");
		goto done;
	    }
	    outputParameters.channelCount = 2;                     /* stereo output */
	    outputParameters.sampleFormat =  PA_SAMPLE_TYPE;
	    outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device )->defaultLowOutputLatency;
	    outputParameters.hostApiSpecificStreamInfo = NULL;
	    printf("\n=== Now playing back ===\n"); fflush(stdout);
	    err = Pa_OpenStream(
		      &stream,
		      NULL, /* no input */
		      &outputParameters,
		      SAMPLE_RATE,
		      FRAMES_PER_BUFFER,
		      paClipOff,      /* we won't output out of range samples */
		      playCallback,
		      &data );
	    if( err != paNoError ) goto done;
	    if( stream ) {
		err = Pa_StartStream( stream );
		if( err != paNoError ) goto done;
		fflush(stdout);
		while( ( err = Pa_IsStreamActive( stream ) ) == 1 ) Pa_Sleep(100);
		if( err < 0 ) goto done;
		err = Pa_CloseStream( stream );
		if( err != paNoError ) goto done;
		printf("Done!\n"); fflush(stdout);
	    }

/*===================================Pre Emphasising===================================*/

	    float *preEmphasisedData = (float *)malloc(sizeof(float) * numSamples);
	    preEmphasis(data.recordedSamples, preEmphasisedData, numSamples);
	    FILE *PreEmph = fopen("Pre Emphasised Data.txt", "w");
	    for (i = 0; i < numSamples; i++) {
		fprintf(PreEmph, "%f\n", preEmphasisedData[i]);
	    }
	    fclose(PreEmph);
	    printf("Pre Emphasising Done!\n");

/*=========================Maximum Amplitude & Average=========================*/

	    max = 0;
	    average = 0.0;
	    for( i=0; i<numSamples; i++ ) {
		val = preEmphasisedData[i];
		if( val < 0 ) val = -val;
		if( val > max ) {
		    max = val;
		}
		average += val;
	    }

	    average = average / (double)numSamples;

	    printf("     Sample Max Amplitude = "PRINTF_S_FORMAT"\n", max );
	    printf("     Sample Average = %lf\n", average );

/*===================================Hamming===================================*/

	    float *usefulData = (float *)malloc(sizeof(float) * numSamples);
	    int len;
	    for (i = 0, len = 0; i < numSamples; i++) {
		if (preEmphasisedData[i] > 0.02 || preEmphasisedData[i] < -0.02) {
		    usefulData[len++] = preEmphasisedData[i];
		}
	    }
	    float *hammingData = (float *)malloc(sizeof(float) * len);
	    hamming(usefulData, hammingData, len);
	    free(preEmphasisedData);
	    free(usefulData);
	    printf("Hamming Done!\n");

/*===============Discrete Fourier Transform Implementation===============*/

	    float *DFTrecordedSamples = (float *)malloc(sizeof(float) * len);
	    computeDFT(hammingData, DFTrecordedSamples, len);
	    free(hammingData);
	    FILE *DFT = fopen("DFTsample.txt", "w");
	    for (i = 0; i < len; i++){
		fprintf(DFT, "%f\n", DFTrecordedSamples[i]);
	    }
	    fclose(DFT);
	    free(DFTrecordedSamples);
	    printf("Discrete Fourier Transform Done!\n");

/*=========================MFCC Library Implementation=========================*/

	    double *spectrum = (double *)malloc(sizeof(double) * 200000);
	    int j = 0;
	    unsigned int coeff;
	    double mfcc_result;
	    FILE *mfcc = fopen("DFTsample.txt", "rb");
	    while(fscanf(mfcc, "%lf", &spectrum[j]) != EOF) {
		j++;
	    }
	    fclose(mfcc);
	    FILE *final = fopen(filename,"w");
	    float MFCC[13];
	    for(coeff = 0; coeff < 13; coeff++) {
		mfcc_result = GetCoefficient(spectrum, 44100, 48, len, coeff);
		MFCC[coeff] = mfcc_result;
		printf("%f\n", mfcc_result);
		fprintf(final, "%f\n", mfcc_result);
	    }
	    fclose(final);
	    free(spectrum);
	    printf("MFCC Done!\n");


	    done:
	        Pa_Terminate();
	        if( data.recordedSamples ) {
		    free( data.recordedSamples );
	        }
	        if( err != paNoError ) {
		    fprintf( stderr, "An error occured while using the portaudio stream\n" );
		    fprintf( stderr, "Error number: %d\n", err );
		    fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
		    err = 1;
	        }
	        return err;
}



int main() {
    while (1) {	
	system("clear");
	printf("Choose what do you want to do : \n1.Record sample\n2.Record & compare with the previous sample\n3.Exit & Show result\n");
	int get;
	scanf("%i", &get);
	if (get == 3) {
	    system("clear");
	    break;
	}
	else if (get == 1) {
	    Record("MFCC Final 1.txt");
	}
	else if (get == 2) {
	    Record("MFCC Final 2.txt");
	}
    }
/*==============================DTW==============================*/

    float DTWresult = DTW("MFCC Final 1.txt","MFCC Final 2.txt",13,13,1,"DTWdebug");
    if( DTWresult <= 15 ){
	printf("They're possiblly match!\n");
    }
    else{
	printf("not match!\n");
    }
    return 0;
}

