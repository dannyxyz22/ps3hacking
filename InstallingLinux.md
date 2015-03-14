Check this link for InstallingFedora9

# Introduction #

The chosen distribution was Yellow Dog Linux. Most of the details on how to install were taken from IBM's [page](http://www.ibm.com/developerworks/power/library/pa-linuxps3-1/)


# Steps #

1. Download Yellow Dog Linux from TerraSoft. http://www.terrasoftsolutions.com/support/downloads/
File yellowdog-5.0.1-phoenix-20070511-PS3.iso was used for the installation.

edit: try this link for Fedora 9 http://www.bohmer.net/ps3_stuff/install-fedora-8-on-PS3.html

2. Format Ps3 so that 10 GB are reserved for Ps3 OS and 50 GB for Linux.

3. Upgrade Ps3 firmware to version 2.10 using pen drive.

4. Install linux through kboot and **installtext** command.

# Hints #

Holding the power on/off Ps3 touch button switches back to GameOS mode, in case someone is stuck in Other OS mode.

Ps3 changes menu color along the day, don't be worried.

# Problems #

When using petitboot 0.2 from official site, I had the following error:
"Corrupted data". The problem is some kind of error with browsers or with the server of petitboot. The file should have around 3MB while it has around 7MB when this is done. The solution was downloading it through linux's command wget. It's simple, just go to the command line and type:
> wget http://ozlabs.org/~jk/projects/petitboot/downloads/bin-0.2/otheros.bld
The downloaded file should have around 3MB.

After installing kboot, during the boot, the message
_INPUT IRQ STATUS -32 RECEIVED_ was looped because of not being able to use the USB keyboard. After reading some threads it was discovered that some keyboards were not supported by default procedures. Trying a different USB keyboard, with a couple less features, the message was gone and installation proceeded smooth.