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


using namespace DAES67;

#define	SHM_NAME	"DepConvolver"

ThreadedDSP DSP;
float Filt[DSP.MaxLength]={}; 
float fVol = 1.0F;
void SetVol(float f) { fVol = f; for (int n=0; n<DSP.Outputs(); n++) DSP.SetGainOut(n,f); };
 
static bool g_running = true;
static void signal_handler(int sig)
{
	(void) sig;
	g_running = false;
	signal(SIGINT, signal_handler);
}

int main(int argc, char * argv[])
{
	const char *sFile="";
	bool bReset=false;

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
		else { sFile=argv[0]; break; };
	}

	signal(SIGINT, signal_handler);

	if (!DSP.Attach(SHM_NAME) )
	{	std::cerr << "Unable to attach shared memory for Convolver at " << SHM_NAME << std::endl; exit(0); };

	if (bReset) 
	{
		SetVol(0.0F);
		for (int s=0; s<DSP.MaxSets; s++) 
		{
			DSP.SetFilterSet(s);
			for (int n=0; n<DSP.Inputs(); n++) for (int m=0; m<DSP.Outputs(); m++) DSP.LoadFilter(n,m); 
		}
	}

	if (sFile[0])
	{
		SetVol(0.0F);
		char filename[512];
		for (int s=0; s<DSP.MaxSets && g_running; s++) 
		{
			snprintf(filename,512,"%s_%04d.txt",sFile,s);
			std::ifstream file(filename);  
			if (!file.is_open()) continue;
			printf("SET %4d  FILTERS ",s);
			DSP.SetFilterSet(s);
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
						printf(".");
					} 
				}
			}
			printf("\n");
			usleep(10000);
		}
	}
	DSP.SetFilterSet(0);
	SetVol(1.0F);

	struct {
		int  modes;			// Number of modes 					1 .. 4
		int  snaps;			// Number of snaps in each mode 	1 .. 8
		int  rotation;		// Rotation total steps				>= 1
		int  centre;		// Centre of rotation				0 .. rotation-1
		bool wrap;			// Wrap around						true/false
		bool reset;			// Reset to centre on a snap change true/false
	} Set = { 1, 1, 145, 71, false, false };

	DSP.SetFilterSet(Set.centre);	
  	
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

		int nRot  = Set.centre;
		int nMode = 0;
		int nSnap = 0;
		float fVolume = 1.0F;		// An output volume
		bool bVolume = false;		// Toggle for the volume button (hold and turn)
		bool bVolChanged = false;	// Set if the volume has changed

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
		message[1] = 0x36;
		midiout.sendMessage(&message);
		message[1] = 0x36+4;
		midiout.sendMessage(&message);

		while (g_running && DSP.Running())
		{
			std::vector<unsigned char> message;
		 	midiin.getMessage(&message);
			if (message.size()>0)
			{	
				printf("Got message   Length %2d   %02x %02x %02x %02x %02x\n",(int)message.size(),message[0],message[1],message[2],message[3],message[4]);

				if (message[0]==0x90 && message[1]==0x43 && message[2]==0x7F)	nRot = Set.centre;	// DIM button resets the rotation

				if (message[0]==0x90 && message[1]==0x44)
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
					int mode=nMode, snap=nSnap;

					if (button>=0 && button<Set.modes)
					{
						message[1] = 0x36 + nMode;
						message[2] = 0x00;
						midiout.sendMessage(&message);
						message[1] = 0x36 + button;
						message[2] = 0x7F;
						midiout.sendMessage(&message);
						mode = button;
					}
					if (button>=4 && button<4+Set.snaps)
					{
						message[1] = 0x36 + nSnap + 4;
						message[2] = 0x00;
						midiout.sendMessage(&message);
						message[1] = 0x36 + button;
						message[2] = 0x7F;
						midiout.sendMessage(&message);
						snap = button-4;
					}

					if (mode!=nMode || snap!=nSnap)
					{
						nMode = mode;
						nSnap = snap;
						if (Set.reset) nRot = Set.centre;
					}

				}

				DSP.SetFilterSet(nRot);				
				printf("MODE %d   SNAP %2d    ROT %2d\n",nMode,nSnap,nRot);
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

