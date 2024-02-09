#include <RtMidi.h>


void arc(void)
{
  	RtMidiIn  midiin;
  	RtMidiOut midiout;
 
	Debug("Starting ARC Thread");
	std::this_thread::sleep_for(std::chrono::milliseconds(1000));
	while (g_running)
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

		while (g_running)
		{
			std::vector<unsigned char> message;
		 	midiin.getMessage(&message);
			if (message.size()>0)
			{	
				printf("Got message   Length %2d   %02x %02x %02x %02x %02x\n",(int)message.size(),message[0],message[1],message[2],message[3],message[4]);
				midiout.sendMessage(&message);
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
	Debug("Exiting ARC Thread");
}
