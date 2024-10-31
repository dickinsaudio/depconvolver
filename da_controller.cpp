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
		if (n<sizeof(fFocus)/sizeof(float)) DSP.SetGainOut(n,f * powf(10.0, -abs(nFocus)*(fdB[n]+0.1)/20.0F));		// Boost the solo channels by up to 10dB
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
						printf("."); fflush(stdout);
						nFilters++;
						nTaps += nLength;
					} 
				}
			}
			printf("   %d\n",nFilters);
			file.close();
		}
		DSP.SetFilterSet(0,g);
	}
	printf("\nLOADED TOTAL OF %d TAPS\n",nTaps);
	SetVol(1.0F);

	struct {
		int  modes;			// Number of modes 					1 .. 4
		int  snaps;			// Number of snaps in each mode 	1 .. 8
		int  rotation;		// Rotation total steps				>= 1
		int  centre;		// Centre of rotation				0 .. rotation-1
		bool wrap;			// Wrap around						true/false
		bool reset;			// Reset to centre on a snap change true/false
	} Set = { 1, 1, 145, 71, false, false };

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



		
		int nRot   = Set.centre;

		int Row1 = 0;
		int Row2 = 0;
		int Row3 = 0;

		float fVolume = 1.0F;			// An output volume
		bool  bVolume = false;			// Toggle for the volume button (hold and turn)
		bool  bVolChanged = false;		// Set if the volume has changed

		bool  bFocus = false;			// Toggle for the focus button (hold and turn)
		bool  bFocusChanged = false;	// Set if the focus has changed

		int   nReverb = 0;				//

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
		message[1] = 0x3E + Row3; midiout.sendMessage(&message);

		while (g_running && DSP.Running())
		{
			std::vector<unsigned char> message;
		 	midiin.getMessage(&message);
			if (message.size()>0)
			{	
				//printf("Got message   Length %2d   %02x %02x %02x %02x %02x\n",(int)message.size(),message[0],message[1],message[2],message[3],message[4]);
					
				if (message[0]==0x90 && message[1]==0x42)						// TALK BUTTON for change or reset focus
				{
					if (message[2]==0x7F) { bFocus = true;  bFocusChanged = false; };
					if (message[2]==0x00) { bFocus = false; if (!bFocusChanged) { nFocus = 0; SetVol(fVol); }; };
				}

				if (message[0]==0x90 && message[1]==0x43 && message[2]==0x7F)	nRot = Set.centre;	// SPKRB button resets the rotation

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
						nFocus = std::max(-200,std::min(200,nFocus));
						bFocusChanged = true;
						SetVol(fVol);
					}
					else
					{
						if (message[2] & 0x40) nRot -= 1;
						else                   nRot += 1;
						if (Set.wrap) nRot = (nRot+Set.rotation)%Set.rotation;
						else          nRot = std::max(0,std::min(Set.rotation-1,nRot));
					}
				}

				if (message[0]==0x90 && message[2]==0x7F)
				{
					int button = message[1] - 0x36;
					int row1 = Row1, row2 = Row2, row3 = Row3;

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
					if (button>=4 && button<8)
					{
						message[1] = 0x36 + Row2 + 4;
						message[2] = 0x00;
						midiout.sendMessage(&message);
						message[1] = 0x36 + button;
						message[2] = 0x7F;
						midiout.sendMessage(&message);
						Row2 = button-4;
					}
					if (button>=8 && button<12)
					{
						message[1] = 0x36 + Row3 + 8;
						message[2] = 0x00;
						midiout.sendMessage(&message);
						message[1] = 0x36 + button;
						message[2] = 0x7F;
						midiout.sendMessage(&message);
						Row3 = button-8;
					}


					if (row1 != Row1 || row2 != Row2 || row3 != Row3)
					{
						if (Set.reset) nRot = Set.centre;
					}

				}

				DSP.SetFilterSet(nRot);				
				printf("ROW1 %d  ROW2 %d  ROW3 %d   FOCUS %4d   ROT %3d   VOL %5.3f \n",Row1,Row2,Row3,nFocus,nRot,fVol);
			}	
			else std::this_thread::sleep_for(std::chrono::milliseconds(1));
			if  (midiin.getPortCount() != nPorts)
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

