#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <tcl/tcl.h>

#define WIDDY 882
#define THRESHY 0.1

unsigned char *raw = NULL, *newraw = NULL;
double * pitches=NULL, *rmses=NULL, * speeds=NULL;
long int * inmarkers=NULL, *outmarkers=NULL;

unsigned long int samples, bytes, desiredLength, frames, newSamples=0;
int input_markerIndex = 0, output_markerIndex = 0;

int sampleAt(int index, int chooseLeft) {
    int val;
    if( index<0 || index>=samples ) return 0;
    if (chooseLeft)
        val = raw[index * 4] + (raw[index * 4 + 1] << 8);
    else
        val = raw[index * 4 + 2] + (raw[index * 4 + 3] << 8);
    if (val > 32767)
        val -= 65536;
    return val;
}

void sampleSet(int index, int chooseLeft, int val) {
    int offs;
    if( index<0 || index>=samples ) return;

    offs = chooseLeft ? 0 : 2;

    raw[index*4 + offs] = val&0xff;
    raw[index*4 + offs + 1] = (val&0xff00)>>8;
}

void samplePut(int index, int chooseLeft, int what, int maxSamp) {
    int val;
    int pos = chooseLeft ? 0 : 2;

    if( index<0 || index>=maxSamp ) return;

    if (chooseLeft)
        val = newraw[index * 4] + (newraw[index * 4 + 1] << 8);
    else
        val = newraw[index * 4 + 2] + (newraw[index * 4 + 3] << 8);
    if (val > 32767)
        val -= 65536;

    val += what;

    if( val>32767 || val<-32768 ) 
	fprintf(stderr,"samplePut at %d out of bounds (%d)\n", index, val );

    newraw[index * 4 + pos] = (val & 0xff);
    newraw[index * 4 + pos + 1] = ((val >> 8) & 0xff);
}

double dotsq(int index, int offset) {
    double sum = 0.0, x, y;
    for (int i = index; i < index + WIDDY; i++) {
        x = (double)sampleAt(i, 1);
        y = (double)sampleAt(i + offset, 1);
        sum += (x - y) * (x - y);
    }
    return sum;
}

/*  makes sense to set 0 phase at time t=0
 *
 */
void fund(int index, double period, double * eye, double * queue) {
    double sum = 0.0, x, y;
    double inph=0.0, quph=0.0, omega, freq;
    int arraylen;

    arraylen = ((int)(WIDDY/period)); 
    arraylen = ((int)(arraylen*period));

    for (int i = index; i < index + arraylen; i++) {

	omega = 6.283185307179586*i/period;
		
        x = (double)sampleAt(i, 1);
	inph += x*cos(omega);
	quph += x*sin(omega);
    }
    *eye = 2*inph/arraylen;
    *queue = 2*quph/arraylen;
}

void defund( int index, double period, double amt ) {
	double eye, que, omega, x;

	fund( index, period, &eye, &que );
	
    	for (int i = index; i < index + WIDDY; i++) {

		omega = 6.283185307179586*i/period;

        	x = (double)sampleAt(i, 1);
	        x -= amt*eye*cos(omega);
	        x -= amt*que*sin(omega);
	
		sampleSet( i, 1, (int) x );		
		sampleSet( i, 0, (int) x );		
	}
}


double gauss( void ) {
        double rsq, v1, v2;

        do {
                v1 = (rand()&65535)/32768.0-1.0;
                v2 = (rand()&65535)/32768.0-1.0;
                rsq = v1*v1+v2*v2;
        } while( rsq==0.0 || rsq>=1.0 );
 
        double fac=sqrt(-2.0*log(rsq)/rsq);
        return v1*fac;
}

void noise( int index, double amt ) {
	double eye, que, omega, x;

    	for (int i = index; i < index + WIDDY; i++) {

        	x = (double)sampleAt(i, 1);
	        x *= 1.0+amt*gauss();
	
		sampleSet( i, 1, (int) x );		
		sampleSet( i, 0, (int) x );		
	}
}


double pitchy(int index, int *argmin, double *lm, double *lam) {
    int i, tflag = 0;
    double sum = 0.0, min = 1.0, dee, deep, odeep = 1, omin, nmin;
    double a, b, c;

    *argmin = 0;
    for (i = 1; i < WIDDY; i++) {
        dee = dotsq(index, i);
        sum += dee;
        deep = dee * i / sum;

        if (deep < THRESHY) {
            if (tflag == 0)
                tflag = 1;
            if (deep < min) {
                min = deep;
                omin = odeep;
                *argmin = i;
            }
        }
        if (i == (*argmin) + 1)
            nmin = deep;
        if (deep >= THRESHY && tflag > 0)
            break;
        odeep = deep;
    }

    if ((*argmin) != 0) {
        i = (*argmin);
        a = (omin + nmin) / 2 - min;
        b = (min - omin) + a * (1 - 2 * i);
        c = min - b * i - a * i * i;

        *lam = -b / (2 * a);
        *lm = a * (*lam) * (*lam) + b * (*lam) + c;
    } else {
        *lm = 1.0;
        *lam = 0.0;
    }

    return min;
}



void allocate( void ) {

    if (newraw!=NULL) {free(newraw); newraw=NULL; }  // no longer allocate newraw here.
 
    // Recalculate number of frames for new audio
    frames = (samples + WIDDY - 1) / WIDDY;
    
    // Free old arrays if they exist
    if (rmses!=NULL) free(rmses);
    if (pitches!=NULL) free(pitches);
    if (speeds!=NULL) free(speeds);
        
    // Allocate new arrays with correct size
    rmses = (double*)malloc(frames * sizeof(double));
    pitches = (double*)malloc(frames * sizeof(double));
    speeds = (double*)malloc(frames * sizeof(double));

    // Initialize new arrays to default values
    for (int i = 0; i < frames; i++) {
        rmses[i] = 1.0;
        pitches[i] = 1.0;
        speeds[i] = 1.0;
    }

    //Free and reallocate marker arrays
    if (inmarkers!=NULL) free(inmarkers);
    if (outmarkers!=NULL) free(outmarkers);

    desiredLength = 10 * samples / WIDDY + 1;
    inmarkers = (long int*)malloc(desiredLength * sizeof(long int));
    outmarkers = (long int*)malloc(desiredLength * sizeof(long int));
    input_markerIndex = 0;
    output_markerIndex = 0;

    return;
}

long read_wav( char * filename ) {
    unsigned char data[5] = "data";
    FILE * fp;

    fp = fopen(filename,"r");
    do {
        fread(data, 4, 1, fp);
    } while (strcmp((char *)data, "data"));
    fread(data, 4, 1, fp);
    bytes = data[0] + (data[1] << 8) + (data[2] << 16) + (data[3] << 24);

    raw = malloc(bytes);
    fread(raw, 1, bytes, fp);

    samples = bytes / 4;

    allocate();

    return samples;
}

void write_wav( char * filename, char * buf, long samples ) {

	FILE * fp;
	long size;
        unsigned char data[5] = "RIFF";

	fp = fopen( filename, "wb" );
	fwrite( data, 4, 1, fp );

	size = 36 + samples*4;
	data[0] = size&255; data[1] = (size>>8)&255; data[2] = (size>>16)&255; data[3] = (size>>24)&255;
	fwrite( data, 4, 1, fp );

	sprintf( (char *) data, "WAVE" ); fwrite( data, 4, 1, fp );
	sprintf( (char *) data, "fmt " ); fwrite( data, 4, 1, fp );
	data[0] = 16; data[1] = 0; data[2] = 0; data[3] = 0; 
	fwrite( data, 4, 1, fp );
	
	data[0] = 1; data[1] = 0; data[2] = 2; data[3] = 0; 
	fwrite( data, 4, 1, fp );

	size = 44100;
	data[0] = size&255; data[1] = (size>>8)&255; data[2] = (size>>16)&255; data[3] = (size>>24)&255;
	fwrite( data, 4, 1, fp );
		
	size = 176400;
	data[0] = size&255; data[1] = (size>>8)&255; data[2] = (size>>16)&255; data[3] = (size>>24)&255;
	fwrite( data, 4, 1, fp );
		
	data[0] = 4; data[1] = 0; data[2] = 16; data[3] = 0; 
	fwrite( data, 4, 1, fp );
	
	sprintf( (char *) data, "data" ); fwrite( data, 4, 1, fp );
	size = samples*4;
	data[0] = size&255; data[1] = (size>>8)&255; data[2] = (size>>16)&255; data[3] = (size>>24)&255;
	fwrite( data, 4, 1, fp );

	fprintf(stderr, "About to write\n");
	fwrite( buf, samples*4, 1, fp );
	fprintf(stderr, "After to write\n");
	
	fclose(fp);
}

/*  bangdecide:  decide how many times (0 or more) a given input tag/window
 *               should be copied into the output.
 *               
 *               This is based on speed[frame], where speed==1.0 means 
 *               make one copy.  A speed of 2.3 will cause successive frames
 *               to copy something like 3, 2, 2, 2, 3, 2, 2, 3... times 
 *               Basically bang-bang control to get 2.3
 */
int bangdecide( int tag ) {
	static double amount=-1.0, pee=0.0;
	static long count=0, sum=0, base=0;	// that's a lot of statics.

	int out;

	long index = inmarkers[tag];	// get position in wav file
	int framn = index/WIDDY;	// get frame number

	if( amount != speeds[framn] ) { // change in speed!
		amount = speeds[framn];
		base = (int) amount;
		pee = amount-base;
		count = sum = 0;
	}

	out = (sum < pee*count) ? base+1 : base;
	sum += out;
	count++;
	return out;
}

/*  trimangle: triangle synthesis that duplicates or drops frames to control the
 *             output speed  
 *
 *             ignore outmarkers, use only inmarkers, 
 *
 *
 */

long sizeOfNew( void ) {
	int i;
	long size=0;
	double factor=0;

	for( i=0; i<frames; i++ ) {
		factor += WIDDY * speeds[i];
	}

	size =  (long) (factor);
	size = ((size/WIDDY)+2)*WIDDY;

    //fprintf( stderr, "SON:  newraw = %d, allocated to %ld\n", input_markerIndex, size );

	return size;
}

long trimangle( void ) {
    long samplesPut=0, outMark=0, nrs;
    int k=0, closest;
    // Triangle window synthesis using only input markers

    nrs = sizeOfNew();    // find out in real time how big output buffer needs to be

    if( newraw == NULL || (newSamples)<nrs ) {
	newraw = (unsigned char *)realloc( newraw, 4*(nrs) );   // nrs is temporary size of 
	for( int i=0; i<nrs*4; i++ ) newraw[i]=0;               //     new buffer, in samples
    }
	
   // fprintf( stderr, "TRIMANGLE:  newraw = %p, allocated to %ld\n", newraw, nrs );

    for (int j = 1; j < input_markerIndex - 1; j++ ) {

	int flag = bangdecide( j );  // copy over the window flag times

    	for( int jj=0; jj<flag; jj++ ) {
        	outMark += inmarkers[j]-inmarkers[j-1];

		closest = j;
        	int start = inmarkers[closest - 1];
        	int center = inmarkers[closest];
        	int end = inmarkers[closest + 1];
        	int length = end - start;

        	for (int n = 0; n < length; n++) {
            		double position = (double)n / (length - 1);
            		double multiplier = 1.0 - fabs(2.0 * position - 1.0);

            		int inputIndex = start + n;
            		int outputIndex = outMark - (length / 2) + n;

            		if (inputIndex >= 0 && inputIndex < samples && outputIndex >= 0 && outputIndex < nrs) {
                		int sample = sampleAt(inputIndex, 1);
                		int frame = inputIndex / WIDDY;
                		double rmsScale = (frame < frames) ? rmses[frame] : 1.0;
                		int value = (int)(sample * multiplier * rmsScale);
                		samplePut(outputIndex, 1, value, nrs);
                		samplePut(outputIndex, 0, value, nrs);
				if( samplesPut<outputIndex ) samplesPut=outputIndex;
            		}
        	}
    	}
    }
    return samplesPut+1;
}

long triangle( void ) {
    long samplesPut=0, nrs;
    int k=0;
    // Triangle window synthesis

    nrs = ((samples/WIDDY) + 2)* WIDDY;    // create needed output buffer on the fly

    if( newraw == NULL || (newSamples)<nrs ) {
	newraw = realloc( newraw, 4*(nrs) );
	for( int i=0; i<nrs*4; i++ ) newraw[i]=0;
    }
   // fprintf( stderr, "TRIANGLE:  newraw = %p, allocated to %ld\n", newraw, nrs );

    for (int j = 1; j < output_markerIndex - 1; j++) {
        long int outMark = outmarkers[j];

        while( k<input_markerIndex-1 && inmarkers[k+1] <= outMark ) k++;
        int minDist = outMark-inmarkers[k];
        int closest = k;
        if( minDist > (inmarkers[k+1]-outMark) ) {
            minDist = inmarkers[k+1]-outMark;
            closest = k+1;
        }

//	fprintf( stderr, "Tag j=%d (%ld), closest=%d (%ld), numSamples=%ld \n", j, outMark, closest, inmarkers[closest], samples );

        int start = inmarkers[closest - 1];
        int center = inmarkers[closest];
        int end = inmarkers[closest + 1];
        int length = end - start;

        for (int n = 0; n < length; n++) {
            double position = (double)n / (length - 1);
            double multiplier = 1.0 - fabs(2.0 * position - 1.0);

            int inputIndex = start + n;
            int outputIndex = outMark - (length / 2) + n;

            if (inputIndex >= 0 && inputIndex < samples && outputIndex >= 0 && outputIndex < nrs) {
                int sample = sampleAt(inputIndex, 1);
                int frame = inputIndex / WIDDY;                    
                double rmsScale = (frame < frames) ? rmses[frame] : 1.0;  
                int value = (int)(sample * multiplier * rmsScale);
                samplePut(outputIndex, 1, value, nrs);
                samplePut(outputIndex, 0, value, nrs);
		if( samplesPut<outputIndex ) samplesPut=outputIndex;
            }
        }
    }
	fprintf(stderr, "triangle done, samplesPut=%ld\n", samplesPut);
    return samplesPut+1;
}

// RMS Computation
double RMSAt(int index, int chooseLeft){
  double sum = 0.0;
  int count = 0;
  for (int i = index; i< index + WIDDY && i<samples; i++){
    int s = sampleAt(i, chooseLeft);
    sum += (double)(s*s);
    count++;
  }
  if (count == 0) return 0.0;
  return sqrt(2*sum/(double)count);
}


int SetInputMarkers(void){

  int i, j, val, posa, delta, newdelta, input_marker = 0, output_marker = 0, def_per = 100, percent = -25, per, mod_mark;
  double score, lm, lam, lfreq, lpfreq;

  for (i = 0; i < samples - WIDDY; i += WIDDY) {
     // fprintf(stderr, "before pitch\n" );
      score = pitchy(i, &val, &lm, &lam);
      lfreq = lpfreq = 0.0;
     // fprintf(stderr, "pitch %d\n", val );

      if (lam != 0.0) {
          lfreq = 44100.0 / lam;
          lpfreq = 44100.0 / val;
      }

      while (input_marker < i + WIDDY && input_markerIndex < desiredLength) {
          if (val > (WIDDY/10))
              def_per = val;
          input_marker += def_per;
          if (input_marker >= samples) break;
          inmarkers[input_markerIndex++] = input_marker;
      }
    //  fprintf(stderr, "frame %d of %ld: %d\n", i, samples, input_marker);
    }
      fprintf(stderr, "frames done\n" );
return input_markerIndex;
}


int setOutputMarker(void){
  int i, j, val, pos, posa, delta, newdelta, input_marker = 0, output_marker = 0, def_per = 100, percent = -25, per, mod_mark;
  double score, lm, lam, lfreq, lpfreq;
  int minnie;

  outmarkers[0] = inmarkers[0];
  output_markerIndex = 1;

  for( ;; ) {
      posa = outmarkers[output_markerIndex-1];  // posa is most recently set outmarker

	 //fprintf( stderr, "out[%d] = %d %ld\n ", output_markerIndex-1, posa, desiredLength );
     //fprintf( stderr, "new delta: %d \n", newdelta);

      // find inmarker closes to posa, store in mod_mark;
      minnie = 100000;
      for( int k=mod_mark=0; k<input_markerIndex; k++ ) {
	
		pos = labs(posa - inmarkers[k] );
		if( pos < minnie ) {
			minnie = pos;
			mod_mark = k;
		}
      }
	// fprintf( stderr, "best match in[%d] = %ld samples = %ld \n", mod_mark, inmarkers[mod_mark], samples );

      if (mod_mark >= input_markerIndex-1) {
              delta = inmarkers[input_markerIndex-1]
                    - inmarkers[input_markerIndex-2];
      } else {
              delta = inmarkers[mod_mark+1]
                    - inmarkers[mod_mark];
      }
	// fprintf( stderr, "inmarkers[%d]=%ld ", mod_mark+1, inmarkers[mod_mark+1] );
      // get multiplier
      newdelta = (int) ( round( ((double)delta) / pitches[(inmarkers[mod_mark]/WIDDY)] ));

	// fprintf( stderr, "%d--%d\n", delta, newdelta );

      if( posa+newdelta >= samples ) break;
      outmarkers[output_markerIndex++] = posa+newdelta;
  }
	 fprintf( stderr, "\n\n" );
  return output_markerIndex;
}



/**  Tcl boilerplate code below:
 **
 **  new main function that launches the interpreter;
 **  Cmd functions to wrap sample, pitchy and readwav commands
 **  app_init to register these new commands
 **
 **  The functions below may look complex, but they're just
 **  housekeeping, that invoke calls to sampleAt(), read_wav()
 **  and pitchy(), feed the arguments from the user, and package
 **  up the return values back to the user.
 **/

int Tcl_AppInit(Tcl_Interp *intp);

int main( int argc, char ** argv ) {
	Tcl_Main( argc, argv, Tcl_AppInit);
	return 0;
}

/*  Wrapper to call sampleAt()
 *
 */
int MySample_Cmd(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {

    int chooseLeft=1, samp;
    long int index;
    // Check for correct number of arguments, get arguments from user
    switch(objc) {
	case 3:  Tcl_GetIntFromObj(interp, objv[2], &chooseLeft);
        case 2:  Tcl_GetLongFromObj(interp, objv[1], &index);
		 break;
        default:
        Tcl_WrongNumArgs(interp, 1, objv, "position {left}");
        return TCL_ERROR;
    }

    samp = sampleAt( index, chooseLeft );   // call the actual function

    Tcl_SetObjResult(interp, Tcl_NewIntObj(samp) );  // set return value

    return TCL_OK;
}

/*  Wrapper to call read_wav()
 *
 */
int MyReadWav_Cmd(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {

    int chooseLeft=1, samp;
    long int samplef;
    char * filename;

    // Check for correct number of arguments, get filename from user
    switch(objc) {
	case 2:  filename = Tcl_GetString(objv[1]);
	         break;
        default:
        Tcl_WrongNumArgs(interp, 1, objv, "filename");
        return TCL_ERROR;
    }

    samplef = read_wav( filename );   // call read_wav()

    Tcl_SetObjResult(interp, Tcl_NewLongObj(samplef) ); // return # of samples read

    return TCL_OK;
}

/*  Wrapper for pitchy()
 *
 */
int MyPitchy_Cmd(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {

    int chooseLeft=1, samp;
    int argmin;
    double lm, lam;

    long int index;
    char * filename;

    // Check for correct number of arguments, get argument from user
    switch(objc) {
	case 2:  Tcl_GetLongFromObj(interp, objv[1], &index);
	         break;
        default:
        Tcl_WrongNumArgs(interp, 1, objv, "position");
        return TCL_ERROR;
    }

    pitchy(index, &argmin, &lm, &lam);  // call pitchy()

    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(lam) );  // return period

    return TCL_OK;
}


//RMS Tcl command
int MyRMS_Cmd(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]){
  long int index;
  int chooseLeft = 1;
  double value;

  if (objc < 2 || objc > 3){
    Tcl_WrongNumArgs(interp, 1, objv, "position {left}");
    return TCL_ERROR;
  }

  if (Tcl_GetLongFromObj(interp, objv[1], &index) != TCL_OK) return TCL_ERROR;
  if (objc == 3 && Tcl_GetIntFromObj(interp, objv[2], &chooseLeft) != TCL_OK) return TCL_ERROR;

  value = RMSAt(index, chooseLeft);
  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(value));
  return TCL_OK;
}

int MyFund_Cmd(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]){
  long int index;
  int chooseLeft = 1, argmin;
  double value, pitch, lm, lam;
  double eye, que;

  if (objc != 2 ){
    Tcl_WrongNumArgs(interp, 1, objv, "position {left}");
    return TCL_ERROR;
  }

  if (Tcl_GetLongFromObj(interp, objv[1], &index) != TCL_OK) return TCL_ERROR;
  pitchy(index, &argmin, &lm, &lam);  // call pitchy()
  if( lam > 0.0 ) {
        fund( index, lam, &eye, &que );
  } else {
	value = 0.0;
	eye = que = 0.0;
  }
  
  Tcl_Obj *rptr = Tcl_NewListObj(0, NULL);
  Tcl_ListObjAppendElement(interp, rptr, Tcl_NewDoubleObj(eye));
  Tcl_ListObjAppendElement(interp, rptr, Tcl_NewDoubleObj(que));
  Tcl_SetObjResult(interp, rptr);
  return TCL_OK;
}

int MyDefund_Cmd(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]){
  long int index;
  int chooseLeft = 1, argmin;
  double value, pitch, lm, lam;
  double eye, que;

  if (objc != 3 ){
    Tcl_WrongNumArgs(interp, 1, objv, "position amount");
    return TCL_ERROR;
  }

  if (Tcl_GetLongFromObj(interp, objv[1], &index) != TCL_OK) return TCL_ERROR;
  if (Tcl_GetDoubleFromObj(interp, objv[2], &value) != TCL_OK) return TCL_ERROR;

  pitchy(index, &argmin, &lm, &lam);  // call pitchy()
  if( lam > 0.0 ) {
        defund( index, lam, value );
  } 
  return TCL_OK;
}

int MyNoise_Cmd(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]){
  long int index;
  int chooseLeft = 1, argmin;
  double value, pitch, lm, lam;
  double eye, que;

  if (objc != 3 ){
    Tcl_WrongNumArgs(interp, 1, objv, "position amount");
    return TCL_ERROR;
  }

  if (Tcl_GetLongFromObj(interp, objv[1], &index) != TCL_OK) return TCL_ERROR;
  if (Tcl_GetDoubleFromObj(interp, objv[2], &value) != TCL_OK) return TCL_ERROR;

  noise( index, value );
  return TCL_OK;
}



int MySetRMS_Cmd(ClientData clientData, Tcl_Interp* interp, int objc, Tcl_Obj* const objv[]) {
    
    long frame;
    double value;
    
    if (objc != 3) {
        Tcl_WrongNumArgs(interp, 1, objv, "frame value");
        return TCL_ERROR;
    }

    if (Tcl_GetLongFromObj(interp, objv[1], &frame) != TCL_OK) return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[2], &value) != TCL_OK) return TCL_ERROR;

    if (frame >= 0 && frame < frames) {
        rmses[frame] = value;
        Tcl_SetResult(interp, "RMS value set", NULL);
    }
    else {
        Tcl_SetResult(interp, "Frame index out of range", NULL);
        return TCL_ERROR;
    }

    return TCL_OK;
}


int MyTag_Cmd(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {

	int markers;

    markers = SetInputMarkers();
    Tcl_SetObjResult(interp, Tcl_NewLongObj((long)markers) ); // return # of samples read
    return TCL_OK;
}

int MyOutTag_Cmd(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {

    int markers;

    markers = setOutputMarker();
    Tcl_SetObjResult(interp, Tcl_NewLongObj((long)markers) ); // return # of samples read
    return TCL_OK;
}

int MySetLong_Cmd(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {

    long index;
    double num;

    switch(objc) {
	case 3:  Tcl_GetLongFromObj(interp, objv[1], &index);
	         Tcl_GetDoubleFromObj(interp, objv[2], &num);
	         break;
        default:
        Tcl_WrongNumArgs(interp, 1, objv, "position speed (1.0 is normal)");
        return TCL_ERROR;
    }
    speeds[(index/WIDDY)] = num;
    return TCL_OK;
}

int MySetPitch_Cmd(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {

    long index;
    double multy;

    switch(objc) {
	case 3:  Tcl_GetLongFromObj(interp, objv[1], &index);
	         Tcl_GetDoubleFromObj(interp, objv[2], &multy);
	         break;
        default:
        Tcl_WrongNumArgs(interp, 1, objv, "position value");
        return TCL_ERROR;
    }
    pitches[(index/WIDDY)] = multy;
    return TCL_OK;
}

int MyTri_Cmd(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {

    newSamples = triangle();
    Tcl_SetObjResult(interp, Tcl_NewLongObj((long)newSamples) ); // return # of samples read
    return TCL_OK;
}

int MyTrimangle_Cmd(ClientData clientData, Tcl_Interp* interp, int objc, Tcl_Obj* const objv[]) {

    newSamples = trimangle();
    Tcl_SetObjResult(interp, Tcl_NewLongObj((long)newSamples)); // return # of samples read
    return TCL_OK;
}

int MyCommit_Cmd(ClientData clientData, Tcl_Interp* interp, int objc, Tcl_Obj* const objv[]) {

    if (raw!=NULL) {free(raw); raw=NULL;}

    raw = (unsigned char*) malloc(newSamples * 4);
    memcpy(raw, newraw, newSamples * 4);
    samples = newSamples;

    allocate();

    //Tcl_SetResult(interp, "Commit executed succesfully", NULL);
    Tcl_SetObjResult(interp, Tcl_NewLongObj(samples) ); // return # of samples read

    return TCL_OK;
}

int MyWrite_Cmd(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {

    int chooseLeft=1, samp;
    long int samplef;
    char * filename;

    // Check for correct number of arguments, get filename from user
    switch(objc) {
	case 2:  filename = Tcl_GetString(objv[1]);
	         break;
        default:
        Tcl_WrongNumArgs(interp, 1, objv, "filename");
        return TCL_ERROR;
    }

    fprintf( stderr, "wRITEWAVE_CMd:  newraw = %p, ready to write %ld\n", newraw, newSamples );
    write_wav( filename, (char *) newraw, newSamples );   // call read_wav()

    return TCL_OK;
}

/*  This just registers the commands.
 *
 *  To create new commands, just add a line here, and define a wrapper
 *  function like the ones above.
 */
int Tcl_AppInit(Tcl_Interp *intp) {
	if (Tcl_Init(intp) == TCL_ERROR ) return TCL_ERROR;

	Tcl_CreateObjCommand( intp, "readwav", MyReadWav_Cmd,
			(ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );
	Tcl_CreateObjCommand( intp, "sample", MySample_Cmd,
			(ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );
	Tcl_CreateObjCommand( intp, "pitchy", MyPitchy_Cmd,
			(ClientData)NULL, (Tcl_CmdDeleteProc *)NULL );
  Tcl_CreateObjCommand(intp, "RMS", MyRMS_Cmd,
          (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateObjCommand(intp, "fund", MyFund_Cmd,
          (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateObjCommand(intp, "defund", MyDefund_Cmd,
          (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateObjCommand(intp, "noise", MyNoise_Cmd,
          (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateObjCommand(intp, "setrms", MySetRMS_Cmd,
           (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(intp, "tag", MyTag_Cmd,
            (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateObjCommand(intp, "outtag", MyOutTag_Cmd,
            (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateObjCommand(intp, "setpitch", MySetPitch_Cmd,
            (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateObjCommand(intp, "setlong", MySetLong_Cmd,
            (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateObjCommand(intp, "triangle", MyTri_Cmd,
            (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateObjCommand(intp, "trimangle", MyTrimangle_Cmd,
            (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(intp, "commit", MyCommit_Cmd,
            (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
  Tcl_CreateObjCommand(intp, "writewav", MyWrite_Cmd,
            (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

	return TCL_OK;
}
