
/* *************************************************************************  
 * Robert Genega
 * csynth.cpp
 *
 * 10/31/2017
 *
 *
 *
 *
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

/* **************************** CONSTANTS ***************************** */

	const double hz = 44100;
	const double ihz = 1/hz;
	const double nyquist = 22050;
	const double max_amplitude = 32760;

	const int max_samples = 44100*60*60; // 1 hour max length ?



	const double tinnitus = nyquist/1.2;

	const double C4 = 261.626;
	const double Db4 = 277.18;
	const double D4 = 293.66;
	const double Eb4=311.13;
	const double E4=329.63;
	const double F4=349.23;
	const double Gb4=369.99;
	const double G4=392.00;
	const double Ab4=415.30;
	const double A4=440.0;
	const double Bb4=466.16;
	const double B4 =493.88;

	//equal temperament constants
	double TET12 = 1.059463;
	double NTET12 = 1/TET12;



namespace little_endian_io
{
  template <typename Word>
  std::ostream& write_word( std::ostream& outs, Word value, unsigned size = sizeof( Word ) )
  {
    for (; size; --size, value >>= 8)
      outs.put( static_cast <char> (value & 0xFF) );
    return outs;
  }
}
using namespace little_endian_io;




//functions

	//n = number of samples being written
	//freq = frequency of sine wave
	//amplitude = max amplitude of sine wave
	//a = array of samples being written to
	void drawSine(int n, double freq, double amplitude, int *a)
	{
	  double w = 2*3.14159*freq; // angular velocity
	  double t;					 // timestep

	  for (int i = 0; i < n; i++)
	  {
	  	t = i*ihz; 
	    double val = sin( w*t);

	    *a += (int)(amplitude * val);

	    a++;

	  }
	}

	//n = number of samples being written
	//freq = initial frequency of sine wave
	//df = frequency change per sample
	//amplitude = max amplitude of sine wave
	//a = array of samples being written to
	void drawSine(int n, double freq, double freqf, double amplitude, int *a)
	{

	  double df = (freqf - freq)/n;
	  double k = df*hz;

	  for (int i = 0; i < n; i++)
	  {
	  	double t = i*ihz;

	    double val = sin(2*3.14159 * (freq*t + k*t*t/2) );

	    *a += (int)(amplitude * val);

	    a++;
	  }
	}
	 

	void drawSaw(int n, double freq, double amplitude, int *a)
	{
	  double period = hz/freq;

	  double step = 2*amplitude/period;

	  double y = -1*amplitude;

	  double time = 0;


	  for (int i = 0; i < n; i++)
	  {
	  	time++;

	    *a +=(int)( y );

	    y += step;

	    if( time/period >= 1 - ihz)
	    {
	      y = -1*amplitude;
	      time = period - time;
	    }

	    a++;
	  }
	}   

	void drawSquare(int n, double freq, double amplitude, int *a)
	{
	  double period = hz/freq;

	  double time = 0;

	  for (int i = 0; i < n; i++)
	  {
	  	time++;

	    if (time/period < .5)
	    {
	      *a +=(int)( amplitude );
	    }
	    else
	    {
	      *a +=(int)( -1*amplitude );
	    }

	    if( time/period >= 1 - ihz)
	    {
	      time = period - time;
	    }
	    a++;
	  }
	}

	void fadeSignal(int n, double *a)
	{
		double step = 1.0/(n-1);
		double count = 0;
		for (int i = 0; i < n; i++)
		{
			*a += (double)(count);
			count += step;
			a++;
		}
	}

	void reverse(int n, int *a)
	{
      int *left = a;
      int *right = a + n - 1;
      int xtemp;
      int ytemp;

      for (int i = 0; i < n/2; i++)
      {
      	xtemp = *left;
      	ytemp = *right;
      	*left = ytemp;
      	*right = xtemp;
      	left++;
      	right--;
      }
	}

	void sineSignal(int n, double freq, double *a)
	{
	  for (int i = 0; i < n; i++)
	  {
	    double val = sin( (3.14159*2*i*freq)*ihz);

	    *a += (double)(val);

	     a++;
	  }
	} 

	void chirpSignal(int n, double freq, double df, double *a)
	{
	  double k = df*hz;

	  for (int i = 0; i < n; i++)
	  {
	  	double t = i*ihz;

	    *a += sin(2*3.14159 * (freq*t + k*t*t/2) );

	    a++;
	  }
	}

    //alpha < .5
    //alpha > 0
	void lpfilter(int n , double alpha, int *a)
	{

	  int y = alpha * *a;
	  a++;
	  for (int i = 0; i < n; i++)
	  {
	    *a = y + alpha * (*a - y);
	    y = *a;
	    a++;
	  }
	}

	void convolve(int n, int *a, int *b)
	{
	  for (int i = 0; i < n; i++)
	  {
	    *a = *a * *b;
	    a++;
	    b++;
	  }
	}

	void convolve(int n, int *a, double *b)
	{
	  for (int i = 0; i < n; i++)
	  {
	    *a = (int)(*a * *b);
	    a++;
	    b++;
	  }
	}

	void convolve(int n, double *a, double *b)
	{
	  for (int i = 0; i < n; i++)
	  {
	    *a = (*a * *b);
	    a++;
	    b++;
	  }
	}

	void fadein(int n, int *a)
	{
		if (n == 0)
		{
			return;
		}
	  double dv = 1.0/n;
	  double filt = 0;
	  for (int i = 0; i< n; i++)
	  {
	  	*a = (*a) * filt;
	  	filt+=dv;
	  	a++;
	  }
	}

	void fadeout(int n, int *a)
	{
		if (n == 0)
		{
			return;
		}
	  double dv = 1.0/n;
	  double filt = 0;
	  a+=n-1;
	  for (int i = 0; i< n; i++)
	  {
	  	*a = (*a) * filt;
	  	filt+=dv;
	  	a--;
	  }
	}

    //average each amplitude with itself and its neighbors
    //n = number of samples being smoothed
    //d = number of smoothing iterations
	void smooth(int n, int *a, int d)
	{
		int x;	

		for (int i = 0; i < d; i++)
		{
		    int *curr = a;
		    int prev = *a;
			for (int z = 0; z < n-1; z++)
			{
				x = *curr;
				int *next = curr+1;
				*curr = (prev + x + *next)/3;
				prev = x;
				curr++;
				
			}
			*curr = (*curr + *curr + prev)/3 ;
		}
	}

    //sum a and b into a
	void sum(int n, int *a, int *b)
	{
		for (int i = 0; i < n; i++)
		{
			*a += *b;
			a++;
			b++;
		}
	}

    //copy n samples from a, and paste d consecutive times
	void copyPaste(int n, int *a, int d)
	{
		int *b = a + n;

		for (int i = 0; i < d; i++)
		{			
			for (int z = 0; z < n; z++)
			{
			  *b = *a;
			  a++;
			  b++;
		    }
		    a-=n;
		}

	}
 
    //n = number of samples being delayed
	//m = delay (in samples)
	//a = input wave
	void delay(int n, int m, int *a)
	{
	  int *b = new int[n+m];

	  int *ai = a;
	  int *bi = b;

	  b+=m;
	  for (int i = 0; i < n; i++)
	  {
	    *b += *a;
	    a++;
	    b++;
	  }
	  sum(n, ai, bi);
	  a = ai;
	}

	//creates random white noise with maximum amplitude of thresh
	void drawNoise(int n, int thresh, int *a)
	{
	  srand(time(NULL));

	  for (int i = 0; i < n; i++)
	  {
	  	int z = rand() % (thresh+1);

	  	z *= 2;

	  	z-= thresh;

	  	*a += z;

	  	

	  	a++;
	  }
	}

	void chorus(int n, int *a)
	{

	}

	void flange(int n, int *a)
	{

	}

	void hpfilter(int n, int *a)
	{

	}

	void reverb(int n, int *a)
	{

	}

	void adsr(int n, double *a)
	{

	}

	//converts s seconds into samples
	int sectosamp(double s)
	{
		return (int)(s * hz);
	}

	//converts n samples to seconds
	double samptosec(int n)
	{
		return n*ihz;
	}

    //n 			= number of samples being granulized
    //grainsize 	= number of samples per grain
    //grainSpace 	= number of samples inbetween in each grain
    //fadeSize 		= number of samples for the fade in/fade out of each grain (part of the total grain size)
    //a 			= input wave
    //b 			= output wave (make sure is long enough)
	void granulize(int n, int grainSize, int grainSpace, int fadeSize, int *a, int *b)
	{
	  //double *filtera = new double[fadeSize];

	  //sineSignal(fadeSize, hz/(2*fadeSize), filtera );

	  for (int i = 0; i < n; i++)
	  {
	  	if (i%(grainSpace + grainSize) == 0)
	  	{
	      fadein(fadeSize, a);
	  	  //convolve(fadeSize, a, filtera);

	      *b = *a;
	      b++;
	      a++;
	  	}
	  	else if (i %(grainSpace + grainSize) == grainSize-fadeSize)
	  	{
	      fadeout(fadeSize, a);

	      *b = *a;
	      b++;
	      a++;
	  	}
	  	else if(i % (grainSpace + grainSize) > 0 &&  i % (grainSpace + grainSize) < grainSize )
	  	{
	  		*b = *a;
	  		b++;
	  		a++;
	  	}
	  	else
	  	{
	  		*b = 0;
	  		b++;
	  	}
	  }
	}


	void wavetable(int n, int *a, int *b)
	{
	}

	//n 		= number of samples being printed
	//offset 	= index of first sample being printed
	//a 		= input wave
	//pixel 	= pixel array being printed to
	//height 	= height of pixel array
	//width 	= width of pixel array
	void printWave(int n, int offset, int *a, unsigned char *pixel, int height, int width)
	{
	  a += offset;

      int xres = n/width;

      int yres = 2*max_amplitude/height;

      int xsum;

      int xcount = 0;

      //set background to grey
      for (int i = 0; i < height*width*3; i++)
      {
       
        pixel[i] = 50;

      }


      for (int x = 0; x < width; x++)
      {
      	xsum = 0;

      	for (int j = 0; j < xres; j++)
      	{
          xsum = max_amplitude + a[xcount];
          xcount++;
        }	
        
        for (int y = 0; y < height; y++)
        {
        	//left border
        	if (x == 0)
        	{
        		*pixel = 255;
        		pixel++;
        		*pixel = 0;
        		pixel++;
        		*pixel = 255;
        		pixel++;
        	}
        	//right border
        	else if (x == width - 1)
        	{
        		*pixel = 255;
        		pixel++;
        		*pixel = 0;
        		pixel++;
        		*pixel = 0;
        		pixel++;
        	}
        	//bottom border
        	else if (y == 0)
        	{
        		*pixel = 0;
        		pixel++;
        		*pixel = 255;
        		pixel++;
        		*pixel = 0;
        		pixel++;
        	}
        	//top border
        	else if (y == height - 1)
        	{
        		*pixel = 0;
        		pixel++;
        		*pixel = 255;
        		pixel++;
        		*pixel = 255;
        		pixel++;
        	}

        	// white line for amplitude = 0
        	else if (y == height/2)
        	{
        		*pixel = 255;
        		pixel++;
        		*pixel = 255;
        		pixel++;
        		*pixel = 255;
        		pixel++;
        	}

        	//color area under the curve
        	else if (xsum*1.0/yres > height/2.0 && y > height/2 && y < xsum*1.0/yres)
        	{
        		*pixel = 20;
        		pixel++;
        		*pixel = 80;
        		pixel++;
        		*pixel = 200;
        		pixel++;
        		
        	}

        	//color area under the curve
        	else if (xsum*1.0/yres < height/2.0 && y < height/2 && y > xsum*1.0/yres)
        	{
        		*pixel = 20;
        		pixel++;
        		*pixel = 80;
        		pixel++;
        		*pixel = 200;
        		pixel++;
        		
        	}
        	else
        	{
        		pixel+=3;
        	}
        }
      }
	}



int main()
{
  

  char const* fname = "chirp.wav";
  char const* fname2 = "chirp.ppm";

  char const* fname3 = "fractal.ppm";

  double seconds   = 5;      // time
  double BPM = 100;            // beats per minute
  int N = hz * seconds;  // total number of samples
  int printSize = 10000;       //number of samples being printed (make printsize = 1000 for 1 sample per pixel)
  int printOffset = N-10001; 		  //index of first sample printed

  int HEIGHT = 300; //max height of image for now
  int WIDTH = 1000; //max width of image for now 

  double wsf = hz/printSize; // frequency = 1 / 1000 samples



  /* ******************************* DO NOT EDIT ******************************************  */
  
	  ofstream f( fname, ios::binary );
	  // Write the file headers
	  f << "RIFF----WAVEfmt ";     // (chunk size to be filled in later)
	  write_word( f,     16, 4 );  // no extension data
	  write_word( f,      1, 2 );  // PCM - integer samples
	  write_word( f,      2, 2 );  // two channels (stereo file)
	  write_word( f,  44100, 4 );  // samples per second (Hz)
	  write_word( f, 176400, 4 );  // (Sample Rate * BitsPerSample * Channels) / 8
	  write_word( f,      4, 2 );  // data block size (size of two integer samples, one for each channel, in bytes)
	  write_word( f,     16, 2 );  // number of bits per sample (use a multiple of 8)
	  // Write the data chunk header
	  size_t data_chunk_pos = f.tellp();
	  f << "data----";  // (chunk size to be filled in later)
  /* ******************************* DO NOT EDIT ******************************************  */
  


  /* ******************************* AUX INIT SECTION ******************************************  */
	  	  

	  int SPB = (int)(hz*60/BPM); //samples ber beat

	  int numbeats = (int)(N/ SPB);
	 


	  //SCALES
	  //c chromatic scale
	  double cchromatic[12] = { };
	  cchromatic[0] = C4;
	  for (int i = 1; i <12; i++)
	  {
	    cchromatic[i] = cchromatic[i-1]*TET12;
	  }
	 
	  //c major scale
	  double cmajor[7] = { };
	  cmajor[0] = C4;
	  cmajor[1] = cmajor[0]*TET12*TET12;
	  cmajor[2] = cmajor[1]*TET12*TET12;
	  cmajor[3] = cmajor[2]*TET12;
	  cmajor[4] = cmajor[3]*TET12*TET12;
	  cmajor[5] = cmajor[4]*TET12*TET12;
	  cmajor[6] = cmajor[5]*TET12*TET12;


	  //PITCH AND RHYTHM
	  double *lead;
	  lead = new double[numbeats];
	  double *harmony;
	  harmony = new double[numbeats*12];


	  //PROGRAM LEAD
	  double *element = lead;
	 
	  for (int i = 0; i < numbeats; i++)
	  {
	    element++;
	  }

	  for (int i = 0; i < numbeats; i++)
	  {
	    element++;
	  }


	  int *left;  // left channel
	  int *right;  // right channel
	  left = new int[N];
	  right = new int[N];

	  int *leftp = left;  //left channel pointer
	  int *rightp = right;  //right channel pointer
	  int *leftTemp, *rightTemp; //placholders

      //clean up buffer
	  for (int i = 0; i < N; i++)
	  {
	  	if (left[i] != 0)
	  	{
	  		left[i] = 0;
	  	}
	  	if (right[i] != 0)
	  	{
            right[i] = 0;
	  	}
	  }

	  int *filter1, *filter2, *filter3;  //available filters
	  double *filtera, *filterb, *filterc;

	  double f1; //working frequency
	  double a1; //working amplitude
	  int n1; //working sample size
  /* ******************************* END AUX INIT SECTION ******************************************  */




  /* ******************************* MAIN ******************************************  */
	  
	  f1 = C4*8;
	  a1 = max_amplitude/2;
	  n1 = N;

	  drawSine(n1, f1, 0, a1, leftp);

	  //drawSaw(n1, f1, a1, leftp);
	  //drawSquare(n1, f1, a1, leftp);
	  //drawSine(n1, f1, a1, leftp);



	  //sandbox section
	  /*
	  f1 =C4/3;
	  a1 = max_amplitude/2;
	  n1 = N/10;

	  

	  filter1 = new int[n1];
	  filter2 = new int[n1];

	  //clean up buffer
	  for (int i = 0; i < n1; i++)
	  {
	  	if (filter1[i] != 0)
	  	{
	  		filter1[i] = 0;
	  	}
	  	if (filter2[i] != 0)
	  	{
            filter2[i] = 0;
	  	}
	  }


	  //drawSaw(n1, f1, a1, filter1);
	  //drawSaw(n1, f1, a1, filter2);

	  //drawSine(n1, f1, a1, filter1);
	  //drawSine(n1, f1, a1, filter2);

	  //drawSine(n1*100, f1, 0.05, a1, filter1);
	  //drawSine(n1*100, f1, 0.05, a1, filter2);

	  //drawSquare(n1*100, f1, a1, filter1);
	  //drawSquare(n1*100, f1, a1, filter2);

	  //drawNoise(n1/2, a1, filter1);
	  //drawNoise(n1/2, a1, filter2);

	  //fadeout(n1/2, filter1);
	  //fadeout(n1/2, filter2);

	  //smooth(n1, filter1, 200);
	  //smooth(n1, filter2, 200);

	  //sum(n1, leftp, filter1);
	  //sum(n1, rightp, filter2);



      //KICK DRUM
      //
	  drawSine(N/10, f1, -.005, a1, leftp);
	  drawSine(N/10, f1, -.005, a1, rightp);

	  leftp += (int)(N/20);
	  rightp += (int)(N/20);

	  fadeout(n1/2, leftp);
	  fadeout(n1/2, rightp);
	  //

	  //SNARE DRUM
      //
	  drawSine(N/40, 200, a1/1.5, leftp);
	  drawSine(N/40, 200, a1/1.5, rightp);

	  fadeout(N/40, leftp);
	  fadeout(N/40, rightp);

	  drawNoise(N/15, a1/8, leftp);
	  drawNoise(N/15, a1/8, rightp);

	  fadein(50, leftp);
	  fadein(50, rightp);

	  fadeout(N/15, leftp);
	  fadeout(N/15, rightp);


	  leftp += (int)(N/20);
	  rightp += (int)(N/20);

	  fadeout(N/20, leftp);
	  fadeout(N/20, rightp);
	  */




	  //sine wave with linear frequency change
	  /*
	  f1 = cmajor[0]*2;
	  a1 = max_amplitude;
	  n1 = N;
	  double df = f1/(n1-1);

	  drawSine(n1, f1, -df, a1, leftp);
	  drawSine(n1, f1, -df, a1, rightp);
	  leftp+=n1;
	  rightp+=n1;  
	  */




	  //checks for out of bounds amplitudes (errors)
      for (int i = 0; i < N; i++)
      {
      	if ( right[i] > max_amplitude || right[i] < -1*max_amplitude )
        {
          printf("%d: %d\n", i , right[i]);
          right[i] = 0;
        }
      }
      for (int i = 0; i < N; i++)
      {
      	if ( left[i] > max_amplitude || left[i] < -1*max_amplitude )
        {
          printf("%d: %d\n", i , left[i]);
          left[i] = 0;
        }
      }
  /* ******************************* END MAIN ******************************************  */




  /* ******************************* IMAGE WRITE ******************************************  */
	  //
      unsigned char pixels[HEIGHT*WIDTH*3];


      printWave(printSize, printOffset, left, pixels, HEIGHT, WIDTH);

   	  //converts u char * to u char [] [] [] ... i know... whatever...
		  unsigned char pixelsf[HEIGHT][WIDTH][3];
		  int iter = 0;
		  for (int w = 0; w < WIDTH; w++)
		  {
		  	for (int h = HEIGHT-1; h >= 0; h--)
		  	{
		  		for (int c = 0; c < 3; c++)
		  		{
		  			pixelsf[h][w][c] = pixels[iter];
		  			iter++;
		  		}
		  	}
		  }

	  FILE *foo = fopen(fname2, "wb");
	  fprintf(foo, "P6\n%d %d\n%d\n", WIDTH, HEIGHT, 255);
      fwrite(pixelsf, 1, HEIGHT*WIDTH*3, foo);
	  fclose(foo);
	  //
  /* ******************************* END IMAGE WRITE ******************************************  */




  /* ******************************* IMAGE WRITE 2 ******************************************  */
	  /*
	  int HEIGHT2 = 691;
      int WIDTH2 = 1000;
      unsigned char pixels2[HEIGHT2][WIDTH2][3];


      for (int x = 0; x<WIDTH2; x++)
	  {
	    for (int y = 0; y<HEIGHT2; y++)
	    {

		  int i = 0;
		  int max = 100;

		  double a = 0;
		  double b = 0;
							

		  while (((a*a + b*b) < 4) && i < max)
		  {
			//double temp = a*a - b*b + (x/400.0)-2.5;
			//b = 2*a*b + (y/400.0) - 1;

		    double temp = a*a - b*b + (x/400.0)-2.5;
		    b = 2*a*b + (y/400.0) - 1;

			a = temp;

			i++;
	      }
				
		  pixels2[y][x][0] = i%255;
		  pixels2[y][x][1] = (i*20)%255;
	      pixels2[y][x][2] = (2*i)%255;
			
	    }

	  }



	  FILE *fr = fopen(fname3, "wb");
	  fprintf(fr, "P6\n%d %d\n%d\n", WIDTH2, HEIGHT2, 255);
      fwrite(pixels2, 1, HEIGHT2*WIDTH2*3, fr);
	  fclose(fr);
	  */
  /* ******************************* END IMAGE WRITE 2 ******************************************  */




  /* ******************************* DO NOT EDIT ******************************************  */

	  //FINISHED CALCULATING VALUE
	  for (int n = 0; n < N; n++)
	  {
	    write_word( f, left[n], 2 );
	    write_word( f, right[n], 2);
	  }

	  // (We'll need the final file size to fix the chunk sizes above)
	  size_t file_length = f.tellp();

	  // Fix the data chunk header to contain the data size
	  f.seekp( data_chunk_pos + 4 );
	  write_word( f, file_length - data_chunk_pos + 8 );

	  // Fix the file header to contain the proper RIFF chunk size, which is (file size - 8) bytes
	  f.seekp( 0 + 4 );
	  write_word( f, file_length - 8, 4 );
	  return 0; 
  /* ******************************* DO NOT EDIT ******************************************  */
}
