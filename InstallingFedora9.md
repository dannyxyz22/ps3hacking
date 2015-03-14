# Procedures #

Download PowerPC Fedora 9 DVD iso

Download petitboot instead of kboot

Save petitboot as otheros.bld in /PS3/otheros of a pendrive

Install Other OS in Ps3 Menu (format the partition before).

Reboot. Insert DVD. Choose DVD and installation from petitboot.

Make sure you make eth0 available during FC9 installation.

After installation make sure you
```
yum install openmpi openmpi-dev gcc
```

Download the sdk iso from IBM.

Copy it to FC9. Mount it to /mnt/cell for instance

rpm -ivh  cell-install-3.1.0-0.0.noarch.rpm from mounted iso.

Go to /opt/cell/

cellsdk install

# Problems #

## Petitboot with around 7 MB - Corrupted data ##
When using petitboot 0.2 from official site, I had the following error: "Corrupted data". The problem is some kind of error with browsers or with the server of petitboot. The file should have around 3MB while it has around 7MB when this is done. The solution was downloading it through linux's command wget. It's simple, just go to the command line and type:
```
    wget http://ozlabs.org/~jk/projects/petitboot/downloads/bin-0.2/otheros.bld
```

The downloaded file should have around 3MB.
But the version that really worked was this one:
ps3-petitboot-rc-08.06.27.bld
The auto start works really nice in this version so it should automatically get the network working without a problem.

## Root password doesn't work ##
Try using the caps when typing the root password EVEN THOUGH you might not have typed it in caps lock. For some reason the root password was with caps.
A
```
passwd
```
command as root will be enough to get your password back to previous caps.