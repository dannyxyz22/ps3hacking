# Introduction #

This post talks about some ideas given by people who work at wind tunnel in Sao Jose dos Campos - SP. They work on CFD codes that run on a cluster and had already parallelized some CFD simulations.

# Details #

These are topics talked on meeting at wind tunnel:
  * **Metis**: software used for partitioning CFD grid in smaller grids in order to distribute those grids among processors on cluster;
  * **MPICH2**: MPI version used by them. They have never tried another MPI version to make a benchmark and check each performance. MPICH2 already has a Fortran 90 compiler. Besides, they also talked about a software called [Jumpshot](http://www.mcs.anl.gov/research/projects/perfvis/software/viewers/index.htm), which is a Java-based visualization tool for doing postmortem performance analysis;
  * **Language**: They all program in Fortran 77;
  * **CFD open source**: There is an open source CFD community called [frecfd](http://www.frecfd.com);
  * **Ensight - CEI**: software for image view and animation;
  * **Tecplot**: CFD and numerical simulation software package used in post-processing simulation results;
  * **Their code**: As they program in Fortran, they spend a lot of memory to run their software, because it is not used dynamic memory allocation
  * **C++ code**: There was a project that involved many universities and [Embraer](http://www.embraer.com.br) to translate Fortran code to C++ and improve it. That project finished on 2006 and have never been used again by wind tunnel group. Maybe (or not) Embraer group uses that C++ version of CFD code;
  * **F2C**: software for converting Fortran code to C++. It works well for small and simple Fortran codes, but is not so good for complex codes;