Only one Ps3 was available, so the cluster was built upon a Ps3 and a Core2 Duo notebook.

Mounting the nfs was a problem... actually it was not mounted.

Configuring the mpi was easy. Only configuring /etc/openmpi-default-hostfile . The 4 cores were working (2 from Ps3 and 2 from Core2 Duo), as seen from _top_ command.

This test didn't use the spes.

The same file was compiled on both machines and give the same names. Note that the default account for other machines is _root_ so, I had to put the files to be run in /root

One **remarkable** note was: average time for Ps3 to execute the pi example program was 35.67 seconds on one core, with the input 200000000, 29.19 seconds with 2 cores (PPE), notebook was 5.57 seconds with 2 cores, 8.66 seconds with only one.

| **Assembly** | **Seconds** |
|:-------------|:------------|
|Ps3 _one core (PPE)_ | 35.67 |
|Ps3 _two cores (PPE)_ | 29.19 |
|Core2 duo _one core_ | 8.66 |
|Core2 duo _two cores_ | 5.57 |

The four cores were delayed by some password input, so, it ended like: 15.21 seconds.

Programs were run like:
`time mpirun -np 4 time ./pi < in.txt`
in which in.txt = 200000000

So, a Core2 Duo (Intel(R) Core(TM)2 Duo CPU T7300  @ 2.00GHz), seems to be around 4 times faster than the Ps3s PowerPcs. Although differences in compilers are also noticed. Ps3 is running a 2.6.25-14.fc9.ppc64