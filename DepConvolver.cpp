#include <assert.h>
#include <signal.h>
#include <string.h>
#ifdef WIN32
#else
#include <unistd.h>
#endif

#include <iostream>
#include <fstream>
#include <sys/ioctl.h>
#include <stdio.h>
#include <math.h>

#include <ThreadedDSP.hpp>
#include <buffer.hpp>
#include <rtp.hpp>
#include <log.hpp>
#include <histogram.hpp>

using namespace DAES67;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Configuration - size and structure details
//


#define	SHM_NAME	"DepConvolver"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Data structures and Static Variables
//

char sScreen[256*256];
ThreadedDSP DSP;
float fIn [DSP.MaxChannels] = {};
float fOut[DSP.MaxChannels] = {};
bool bKill;
float Filt[DSP.MaxLength]={}; 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Event handlers
//
static bool g_running = true;
static void signal_handler(int sig)
{
	(void) sig;
	g_running = false;
	bKill = true;
	signal(SIGINT, signal_handler);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UI Display - Text Meters
//
#define paint(x,y,c)    { if ((x)<nWidth && (y)<nHeight) s[(nHeight-(y)-1)*(nWidth+1)+(x)]=c; }
#define row(x,y,xx,c)   { for (int __n=(x); __n<(x)+(xx); __n++) paint(__n,y,c); }
#define col(x,y,yy,c)   { for (int __n=(y); __n<(y)+(yy); __n++) paint(x,__n,c); }
void meters(char *s, int nWidth, int nHeight, int nIn, int nOut, float *pfIn, float *pfOut)
{
	memset(s,' ',nHeight*(nWidth+1)); for (int n=0;n<nHeight;n++) { s[n*(nWidth+1)+nWidth]='\n'; };
	s[nHeight*(nWidth+1)]=0;

	int   nBar   = nHeight-3;
	float fScale = 1.0F/5.0F; 

	col(3,3,nBar,'|');
	col(3+nIn+1,3,nBar,'|');
	col(3+nIn+3,3,nBar,'|');
	col(3+nIn+3+nOut+1,3,nBar,'|');
	row(3,2,nIn+2,'-');
	row(3+nIn+3,2,nOut+2,'-');

	for (int n=0; n<9; n++) 
	{ 
		int nH = nHeight-(int)(n*10*fScale+0.5)-1; 
		if (nH<3) break;
		if (n>0) paint(0,nH,'-'); paint(1,nH,'0'+n); paint(2,nH,'0');
		if (n>0) paint(8+nIn+nOut,nH,'-'); paint(9+nIn+nOut,nH,'0'+n); paint(10+nIn+nOut,nH,'0');
	}

	static char S[] = { ' ', '_', '.', 'x', 'X'};

	for (int n=0; n<nIn; n++) 
	{ 
		paint(4+n,0,(n+1)%10+'0');
		paint(4+n,1,(n+1)/10+'0');
		float h = nBar+fScale*20.0F*log10f(pfIn[n]+1E-12F);
		if (h>0)
		{
			col  (4+n,3,(int)h,'X');
			paint(4+n,3+(int)h,S[(int)((h-(int)h)*5)]);
		}
	};
	
	for (int n=0; n<nOut; n++) 
	{ 
		paint(4+nIn+3+n,0,(n+1)%10+'0'); 
		paint(4+nIn+3+n,1,(n+1)/10+'0');
		float h = nBar+fScale*20.0F*log10f(pfOut[n]+1E-12F);
		if (h>0)
		{
			col  (4+nIn+3+n,3,(int)h,'X');
			paint(4+nIn+3+n,3+(int)h,S[(int)((h-(int)h)*5)]);
		}
	};
}

static int nWait;



void dep(void)
{

	DAES67::Buffer daes67;

	if (!daes67.connect("DanteEP")) { fprintf(stderr,"Cannot open Dante\n"); exit(1); };

	struct sched_param param;
	param.sched_priority = sched_get_priority_max(SCHED_FIFO)-60;
    int ret = pthread_setschedparam(pthread_self(), SCHED_FIFO, &param);

	int 	 nTx 	 = daes67.get()->audio.tx;
	int 	 nRx 	 = daes67.get()->audio.rx;
	int32_t  Temp[DSP.MaxBlock] = {};

	daes67.clear();
	while (g_running && DSP.Running())
	{
		uint64_t nPeriod;
		int nPeriodsPerBlock = DSP.BlockSize() / daes67.get()->clock.period;
		for (int n=0; n<nPeriodsPerBlock; n++) nPeriod = daes67.wait(10);

		if (nPeriod==0) { printf("RESETTING IN DEP LOOP\n"); daes67.clear(); continue; };

		int chunk;
		int 		 nBlock             = DSP.BlockSize();
		int          nSamplesPerPeriod  = daes67.get()->clock.period;
		int          nSamplesPerChannel = daes67.get()->audio.samples;
		int          nLatency           = DSP.Latency();
		unsigned int nRxHead = (unsigned int) (((nPeriod-nPeriodsPerBlock)*nSamplesPerPeriod           ) % nSamplesPerChannel);		// Use data one Block before the
		unsigned int nTxHead = (unsigned int) (((nPeriod-nPeriodsPerBlock)*nSamplesPerPeriod + nLatency ) % nSamplesPerChannel);		// actual Dante heads

		if (daes67.get()->clock.periods - nPeriod > 0)
			Histogram(DSP.Chrono_N(0)).add((daes67.get()->clock.periods - nPeriod)*nSamplesPerPeriod); 		// Log the call time against sample clock

		chunk = nSamplesPerChannel - nRxHead;
		if (chunk>=nBlock)  for (int n = 0; n < nRx; n++) DSP.Input (n, (int32_t *)daes67.RX(n)+nRxHead,1);
		else
		{
			for (int n = 0; n < nRx; n++)
			{
				memcpy(Temp,       daes67.RX(n)+nRxHead,    chunk    *sizeof(int32_t));
				memcpy(Temp+chunk, daes67.RX(n),       (nBlock-chunk)*sizeof(int32_t));
				DSP.Input(n, Temp, 1);
			} 
		}
		
		//Debug("PROCESS Loop with data %d",*daes67.RX(0));
		DSP.Process();
		DSP.Finish();

		chunk = nSamplesPerChannel - nTxHead;

		if (chunk>=nBlock)  for (int n = 0; n < nTx; n++) DSP.Output(n, (int32_t *)daes67.TX(n)+nTxHead,1);
		else
		{
			for (int n = 0; n < nTx; n++)
			{
				DSP.Output(n, Temp, 1);
				memcpy(daes67.TX(n)+nTxHead,Temp,         chunk     *sizeof(int32_t));
				memcpy(daes67.TX(n),        Temp+chunk,(nBlock-chunk)*sizeof(int32_t));
			} 
		}

		
/*		chunk = std::min(nSamplesPerChannel-(int)nRxHead,nBlock);
		memcpy(Temp,      (int32_t *)Buffers.getDanteRxChannel(0)+nRxHead,    chunk    *sizeof(int32_t));
		if (chunk<nBlock) memcpy(Temp+chunk,(int32_t *)Buffers.getDanteRxChannel(0),        (nBlock-chunk)*sizeof(int32_t));
		chunk = std::min(nSamplesPerChannel-(int)nTxHead,nBlock);
		memcpy((int32_t *)Buffers.getDanteTxChannel(0)+nTxHead,Temp,         chunk     *sizeof(int32_t));
		if (chunk<nBlock) memcpy((int32_t *)Buffers.getDanteTxChannel(0),        Temp+chunk,(nBlock-chunk)*sizeof(int32_t));
*/		

		int time = (daes67.get()->clock.periods - nPeriod)*nSamplesPerPeriod - nLatency + nBlock;
		Histogram(DSP.Chrono_N(1)).add(time);  

	}
	Debug("DEP Thread Finished");
}


int main(int argc, char * argv[])
{
	std::string g;

	// PARSE COMMAND LINE OPTIONS
	bool bClear=false, bReset=false, bDSP=false, bSilent=false, bLogY=false;
	int  nDSP[]={ 16, 16, 64, 32, 144, 1, 1}, nTestLength=0;
	const char *sFile="";

	setlogmask(LOG_UPTO(LOG_DEBUG));

	while (argc>1)
	{
		argv++; argc--;		// Skip executable
		int nVal;
		if (*argv[0]=='-') switch (tolower(argv[0][1]))
		{
			case 'c' : bClear=true;  break;
			case 'r' : bReset=true;  break;
			case 's' : bSilent=true; break;			
			case 'l' : bLogY=true;   break;			
			case 'd' : bDSP=true;    for (int n=0;n<7;n++) { argv++; argc--; if (!argc) break; sscanf(argv[0],"%d",nDSP+n); }; break; 
			case 'x' : argv++; argc--; if (!argc) break; sscanf(argv[0],"%d",&nVal); nTestLength=nVal; break;  	
			case 'k' : bKill=true; break;
			default: printf("Argument Error : Usage  \nDepConvolver [ -s(ilent) -c(lear) -(r)eset -l(ogY) -d(sp) rx tx block blocks latency filters threads -x testlength -(k)ill ] filter_file\n"); exit(0);
		}
		else { sFile=argv[0]; break; };
	}

	//bDSP = true;
	//nTestLength = 32;

	signal(SIGINT, signal_handler);
	struct winsize w;
	char *s = (char *)calloc(1024*256,sizeof(char));

	if ( ( bDSP && !DSP.Create(nDSP[2],nDSP[3],nDSP[0],nDSP[1],nDSP[5],0,SHM_NAME,nDSP[6],nDSP[4])) ||
	     (!bDSP && !DSP.Attach(SHM_NAME) ) )	
	{	std::cerr << "Unable to " << (bDSP?"create":"attach") << " shared memory for Convolver at " << SHM_NAME << std::endl; exit(0); };

	if (bKill) { DSP.Stop(); exit(0); }; 

	if (DSP.Owner())
	{
		Histogram(DSP.Chrono_N(0)).configure("COUNT OF LATE CALL EVENTS (samples)",0,800);
		Histogram(DSP.Chrono_N(1)).configure("COUNT OF DSP OUTPUT TIMES (samples)",-400,400);
	    std::thread *pthread = new std::thread(dep);
	}

	Filt[0]=1.0F;
	if (bReset) for (int n=0; n<DSP.Inputs(); n++) for (int m=0; m<DSP.Outputs(); m++) DSP.LoadFilter(n,m); 
	if (nTestLength) for (int n=0; n<DSP.MaxFilters; n++) DSP.LoadFilter(((n*(std::max(nDSP[0],nDSP[1])-1))/nDSP[1])%nDSP[0],(n*(std::max(nDSP[0],nDSP[1])-1))%nDSP[1],std::min((int)(sizeof(Filt)/sizeof(float)),nTestLength),Filt);

	if (bClear || bReset) DSP.Chrono_Reset();
	
	std::ifstream file;  
	if (sFile[0])
	{
		std::ifstream file(sFile);  
		if (!file.is_open()) { std::cerr << "WARNING Could not load filters" << std::endl; };
		while (file.is_open() && !file.eof())
		{
			std::string str;
			std::getline(file,str);
			if (str.find("FILTER")!=std::string::npos)
			{
				int nIn=0, nOut=0, nLength=0;
				sscanf(str.c_str(),"FILTER %*s Length= %d In= %d Out= %d",&nLength,&nIn,&nOut);
				if (nLength && nIn && nOut)
				{
					for (int n=0; n<nLength; n++) file >> Filt[n];
					DSP.LoadFilter(nIn-1,nOut-1,nLength,Filt);
					printf("FILTER Loaded IN=%-2d  OUT=%-2d  Length=%6d\n",nIn,nOut,nLength);
				} 
			}
		}
	}
	

	while(g_running && DSP.Running())
	{
		usleep(200000);
		ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
		if (!bDSP && bSilent) { bKill=true; break; }; 
		if (w.ws_row<20 || bSilent) continue;
		for (int n=0; n<DSP.Inputs(); n++) fIn[n] = DSP.PeakIn(n); 
		for (int n=0; n<DSP.Outputs(); n++) fOut[n] = DSP.PeakOut(n);

		meters(s, w.ws_col, w.ws_row-2 - 4*(w.ws_row/5), DSP.Inputs(), DSP.Outputs(), fIn, fOut);
		int ret = system("clear");
		printf("%s [%s] Block %12ld   Taps %8d\n",DSP.Owner()?"Running":"Attached",DSP.SHM_Name(),DSP.Count(),DSP.Taps());
		printf("%s",s);


//		Calls.text(w.ws_row/3, s, true);
//		printf("\n\n%s\n\n",s);


		Histogram(DSP.Chrono_CallTime()).text(w.ws_row/5-2, s, bLogY);
		printf("%s",s);
		Histogram(DSP.Chrono_Load()).text(w.ws_row/5-2, s, bLogY);
		printf("%s",s);
		Histogram(DSP.Chrono_N(0)).text(w.ws_row/5-2, s, bLogY);
		printf("%s",s);
		Histogram(DSP.Chrono_N(1)).text(w.ws_row/5-2, s, bLogY);
		printf("%s",s);
	}
	usleep(50000);
	if (!DSP.Owner() && DSP.Running()) printf("DSP STATE %12ld Rx=%d Tx=%d N=%d M=%d L=%d F=%d T=%d   Filters=%d UsedTaps=%d\n",DSP.Count(), DSP.Inputs(),DSP.Outputs(),DSP.BlockSize(),DSP.Blocks(),DSP.Latency(),DSP.Filters(),DSP.Threads(),DSP.Filters(),DSP.Taps());
	if (!bKill)	std::cerr << "Unexpected Fail : Possibly Insufficient CPU for sustaining DSP" << std::endl;
	Debug("Main loop exited and finishing normally");
	return 0;
}


#ifdef __cplusplus
extern "C"
{
#endif

#include <stdbool.h>

    const char *       ATS_CORE_GIT_HASH_FULL = "";
    const char *       ATS_CORE_GIT_HASH_SHORT = "";
    const bool         ATS_CORE_GIT_STATE_CHANGES = false;
    const char *       ATS_CORE_GIT_COMPONENT_TAG_HASH = "";
    const char *       ATS_CORE_GIT_COMPONENT_TAG_HASH_SHORT = "";
    const char *       ATS_CORE_GIT_COMPONENT_TAG = "";
    const char *       ATS_CORE_GIT_COMPONENT_TAG_VERSION = "";
    const bool         ATS_CORE_GIT_COMMITS_SINCE_COMPONENT_TAG = false;
    const char *       ATS_CORE_GIT_COMPONENT_TAG_VERSION_TYPE = "";
    const char *       ATS_CORE_GIT_COMPONENT_TAG_VERSION_FULL = "";
    const char *	   ATS_CORE_GIT_VERSION_MAJOR = "";
    const char *       ATS_CORE_GIT_VERSION_MINOR = "";
    const char *       ATS_CORE_GIT_VERSION_PATCH = "";
    const char *       ATS_CORE_GIT_VERSION_SUFFIX = "";

#ifdef __cplusplus
}
#endif