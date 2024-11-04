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

#include "ThreadedDSP.hpp"
#include <daes67.hpp>
#include <log.hpp>
#include <RtMidi.h>
#include <math.h>

using namespace DAES67;

#define	SHM_NAME	"DepConvolver"

ThreadedDSP DSP;
float Filt[DSP.MaxLength]={}; 
float fVol = 1.0F;

int   nFocus = 0;				// Focus button

// FOCUS fprintf("%4.3f, ",(abs(mod(a+180,360)-180)/180).^1.5)
// BLIND fprintf("%4.3f, ",(abs(mod(a,360)-180)/180).^1.5)
//                  -45,  -45,  -40,  -40,  -35,  -35,  -25,  -25,  -20,  -20,  -15,  -15,  -10,  -10,   -5,   -5,    5,    5,   10,   10,   15,   15,   20,   20,   25,   25,   35,   35,   40,   40,   45,   45,   60,   60,   60,   60,  120,  120,  120,  120,  240,  240,  240,  240,  300,  300,  300,  300,  -30,    0,   30,   90,  150,  210,  270 };
float fFocus[] = { 0.12, 0.12, 0.10, 0.10, 0.09, 0.09, 0.05, 0.05, 0.04, 0.04, 0.02, 0.02, 0.01, 0.01, 0.00, 0.00, 0.00, 0.00, 0.01, 0.01, 0.02, 0.02, 0.04, 0.04, 0.05, 0.05, 0.09, 0.09, 0.10, 0.10, 0.12, 0.12, 0.19, 0.19, 0.19, 0.19, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.19, 0.19, 0.19, 0.19, 0.07, 0.00, 0.07, 0.35, 0.76, 0.76, 0.35 };
float fBlind[] = { 0.65, 0.65, 0.69, 0.69, 0.72, 0.72, 0.80, 0.80, 0.84, 0.84, 0.88, 0.88, 0.92, 0.92, 0.96, 0.96, 0.96, 0.96, 0.92, 0.92, 0.88, 0.88, 0.84, 0.84, 0.80, 0.80, 0.72, 0.72, 0.69, 0.69, 0.65, 0.65, 0.54, 0.54, 0.54, 0.54, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.19, 0.54, 0.54, 0.54, 0.54, 0.76, 1.00, 0.76, 0.35, 0.07, 0.07, 0.35 };

void SetVol(float f) 
{ 
	fVol = f; 
	float *fdB;
	if (nFocus >= 0) fdB = fFocus;
	else             fdB = fBlind;

	for (int n=0; n<DSP.Outputs(); n++) 
	{
		
		if (n<sizeof(fFocus)/sizeof(float))
		{
		   if (nFocus<0) DSP.SetGainOut(n,f * powf(10.0, -abs(nFocus)*(fdB[n]-0.2)/20.0F));		// Boost the rear channels by up to 20dB
		   else          DSP.SetGainOut(n,f * powf(10.0, -abs(nFocus)*(fdB[n]-0.06)/20.0F));    // Boost the cente channels by up to 6
		}
		else                                DSP.SetGainOut(n,f); 
	}
};
 


static bool g_running = true;
static void signal_handler(int sig)
{
	(void) sig;
	g_running = false;
	signal(SIGINT, signal_handler);
}

int main(int argc, char * argv[])
{
	const char *sFile[DSP.MaxGroups]= {};
	bool bReset=false;
	int  nGroups = 0;


	setlogmask(LOG_UPTO(LOG_DEBUG));

	while (argc>1)
	{
		argv++; argc--;		// Skip executable
		int nVal;
		if (*argv[0]=='-') switch (tolower(argv[0][1]))
		{
			case 'r' : bReset=true;  break;
			default: printf("Argument Error : Usage  \nda_controller [ -r(eset) ] filter_base\n"); exit(0);
		}
		else { sFile[nGroups++]=argv[0]; continue; };
	}

	signal(SIGINT, signal_handler);

	if (!DSP.Attach(SHM_NAME) )
	{	std::cerr << "Unable to attach shared memory for Convolver at " << SHM_NAME << std::endl; exit(0); };

	if (bReset) 
	{
		SetVol(0.0F);
		DSP.ClearFilters();
	}

	int  nTaps = 0;
	int  nSets[nGroups] = {};
	for (int g=0; g<nGroups; g++)
	{
		if (sFile[g]==0) continue;
		SetVol(0.0F);
		char filename[512];
		printf("\nLOADING GROUP %d BASE %s\n",g,sFile[g]);
		for (int s=0; s<DSP.MaxSets && g_running; s++) 
		{
			snprintf(filename,512,"%s_%04d.txt",sFile[g],s);
			std::ifstream file(filename);  
			if (!file.is_open()) continue;
			printf("   SET %4d  FILTERS ",s);
			int nFilters = 0;
			std::this_thread::sleep_for(std::chrono::milliseconds(20));
			DSP.SetFilterSet(s, g);
			std::this_thread::sleep_for(std::chrono::milliseconds(20));
			while (file.is_open() && !file.eof() && g_running)
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
						if (nFilters%10==0) { printf("."); fflush(stdout); };
						nFilters++;
						nTaps += nLength;
					} 
				}
			}
			printf("   %d\n",nFilters);
			nSets[g]++;
			file.close();
		}
		DSP.SetFilterSet(0,g);
	}
	int nRots     = (nSets[0]/4);
	int nRotScale = 1;
	int nRotMax   = nRots*nRotScale;
	int nRotCent  = (nRotMax+1)/2 - 1;

	printf("\nLOADED TOTAL OF %d SETS AND %d TAPS  - INFERRING ROTATIONS %d  CENTRE %d  SCALE %d\n",nSets[0],nTaps,nRots,nRotCent,nRotScale);
	SetVol(1.0F);

	//DSP.SetFilterSet(Set.centre);	
  	
	RtMidiIn  midiin;
  	RtMidiOut midiout;
	std::this_thread::sleep_for(std::chrono::milliseconds(1000));

	while(g_running && DSP.Running())
	{
		int nPorts = midiin.getPortCount();
		for ( unsigned int i=0; i<nPorts; i++ ) 
		{
			Debug("Detected MIDI Input device %s",midiin.getPortName(i).c_str());
			if (midiin.getPortName(i).find("RME ARC") == 0)
			{
				midiin.openPort(i);
				midiin.ignoreTypes(true, true, true);
				Debug("Opened MIDI Input device %s",midiin.getPortName(i).c_str());
				break;
			}
		}
		nPorts = midiin.getPortCount();

		for ( unsigned int i=0; i<midiout.getPortCount(); i++ ) 
		{
			Debug("Detected MIDI Output device %s",midiout.getPortName(i).c_str());
			if (midiout.getPortName(i).find("RME ARC") == 0)
			{
				midiout.openPort(i);
				Debug("Opened MIDI Output device %s",midiout.getPortName(i).c_str());
				break;
			}
		}

		if (!midiin.isPortOpen())
		{
			Debug("No ARC RME Detected - Waiting to scan ports again");
			while (g_running && (midiin.getPortCount() == nPorts))
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(100));
			}
			continue;
		}

		int nRot = nRotCent;
		int Row1 = 0;
		int Row2 = 0;

		bool bWrap = false;
		bool bReset = false;

		float fVolume = 1.0F;			// An output volume
		bool  bVolume = false;			// Toggle for the volume button (hold and turn)
		bool  bVolChanged = false;		// Set if the volume has changed

		bool  bFocus = false;			// Toggle for the focus button (hold and turn)
		bool  bFocusChanged = false;	// Set if the focus has changed

		std::vector<unsigned char> message;
		message.resize(3);
		message[0] = 0x90;
		message[2] = 0x00;
		for (int b=0x36; b<0x45; b++)
		{
			message[1] = b;
			midiout.sendMessage(&message);
		}
		message[2] = 0x7F;
		message[1] = 0x36 + Row1; midiout.sendMessage(&message);
		message[1] = 0x3A + Row2; midiout.sendMessage(&message);

		DSP.SetFilterSet(nRot/nRotScale);		

		while (g_running && DSP.Running())
		{
			std::vector<unsigned char> message;
		 	midiin.getMessage(&message);
			if (message.size()>0)
			{	
				//printf("Got message   Length %2d   %02x %02x %02x %02x %02x\n",(int)message.size(),message[0],message[1],message[2],message[3],message[4]);

				int row1 = Row1, row2 = Row2, rot = nRot;

				if (message[0]==0x90 && message[1]==0x42)						// TALK BUTTON for change or reset focus
				{
					if (message[2]==0x7F) { bFocus = true;  bFocusChanged = false; };
					if (message[2]==0x00) { bFocus = false; if (!bFocusChanged) { nFocus = 0; SetVol(fVol); }; };
				}

				if (message[0]==0x90 && message[1]==0x43 && message[2]==0x7F)	nRot = nRotCent;	// SPKRB button resets the rotation

				if (message[0]==0x90 && message[1]==0x44)											// DIM button for change or reset volume
				{
					if (message[2]==0x7F) { bVolume = true;  bVolChanged = false; };
					if (message[2]==0x00) { bVolume = false; if (!bVolChanged) SetVol(1.0); };
				}
				
				if (message[0]==0xB0 && message[1]==0x10)
				{
					if (bVolume)
					{
						if (message[2] & 0x40) fVol /= 1.0798F;		// 20dB a rotation
						else                   fVol *= 1.0798F;	    
						if (fVol>3.2)   fVol = 3.2;
						if (fVol<0.03)  fVol = 0.03;
						SetVol(fVol);
						bVolChanged = true;
					}
					else if (bFocus)
					{
						if (message[2] & 0x40) nFocus -= 1;
						else                   nFocus += 1;
						nFocus = std::max(-100,std::min(100,nFocus));
						bFocusChanged = true;
						SetVol(fVol);
					}
					else
					{
						if (message[2] & 0x40) nRot -= 1;
						else                   nRot += 1;
						if (bWrap) nRot = (nRot+nRotMax)%(nRotMax-1);
						else       nRot = std::max(0,std::min(nRot,nRotMax-1));
					}
				}

				if (message[0]==0x90 && message[2]==0x7F)
				{
					int button = message[1] - 0x36;

					if (button>=0 && button<4)
					{
						message[1] = 0x36 + Row1;
						message[2] = 0x00;
						midiout.sendMessage(&message);
						message[1] = 0x36 + button;
						message[2] = 0x7F;
						midiout.sendMessage(&message);
						Row1 = button;
					}

					if (row1 != Row1 || row2 != Row2 )
					{
						
						if (bReset) nRot = nRotCent;
					}

				}
				int r = nRot/nRotScale;
				DSP.SetFilterSet(Row1*nRots + r);		
				printf("ROW1 %d  ROW2 %d  ROT %d    FOCUS %3d   NROT %3d   VOL %5.3f    DSP SET %d\n",Row1,Row2,r,nFocus,nRot,fVol,DSP.GetFilterSet());
			}	
			else std::this_thread::sleep_for(std::chrono::milliseconds(1));
			
			static int probe=0;
			if  ((probe++ % 5000 == 0) && midiin.getPortCount() != nPorts)
			{
				Debug("Change in MIDI devices - rescan");
				break;
			}
		}

		midiin.closePort();
		midiout.closePort();
	}
	return 0;
}

