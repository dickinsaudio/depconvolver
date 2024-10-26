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
		for (int s=0; s<DSP.MaxSets; s++) 
		{
			DSP.SetFilterSet(s);
			for (int n=0; n<DSP.Inputs(); n++) for (int m=0; m<DSP.Outputs(); m++) DSP.LoadFilter(n,m); 
		}
	}

	if (sFile[0])
	{
		char filename[512];
		for (int s=0; s<DSP.MaxSets; s++) 
		{
			snprintf(filename,512,"%s_%04d.txt",sFile,s);
			std::ifstream file(filename);  
			if (!file.is_open()) continue;
			DSP.SetFilterSet(s);
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
						printf("SET %4d  FILTER Loaded IN=%-2d  OUT=%-2d  Length=%6d\n",s,nIn,nOut,nLength);
					} 
				}
			}
		}
	}
	DSP.SetFilterSet(0);
	

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

		int nRot  = 0;
		int nMode = 0;
		int nSnap = 0;

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
				if (message[0]==0xB0 && message[1]==0x10)
				{
					if (message[2] & 0x40) nRot -= 1;
					else                   nRot += 1;
				}

				if (message[0]==0x90 && message[2]==0x7F)
				{
					int button = message[1] - 0x36;

					if (button==14) nRot = 0;					// DIM button resets the rotation

					if (button>=0 && button<4)
					{
						message[1] = 0x36 + nMode;
						message[2] = 0x00;
						midiout.sendMessage(&message);
						message[1] = 0x36 + button;
						message[2] = 0x7F;
						midiout.sendMessage(&message);
						nMode = button;
					}
					if (button>=4 && button<12)
					{
						message[1] = 0x36 + nSnap + 4;
						message[2] = 0x00;
						midiout.sendMessage(&message);
						message[1] = 0x36 + button;
						message[2] = 0x7F;
						midiout.sendMessage(&message);
						nSnap = button-4;
					}
				}

				DSP.SetFilterSet((nRot+16)%16);				
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

