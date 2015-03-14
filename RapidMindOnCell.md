# Introduction #
In this wiki, a comparison will be made between parallel program _x_ non-parallel program, both of them running on Cell with Rapidmind. For this, the sum of some numbers will be made as shown below:

| inputArray | outputArray |
|:-----------|:------------|
| { {0, 1, 2} , {1, 2, 3} } |	{(0+1+2) = 3, (1+2+3) = 6 } |
| { {1, 2, 3} , {2, 3, 4} } |	{(1+2+3) = 6, (2+3+4) = 9 } |
| { {2, 3, 4} , {3, 4, 5} } |	{(2+3+4) = 9, (3+4+5) = 12} |
| ... | ... |

This operation will be done many times repeatedly in order to make more calculations and increase execution time so it can be measured.

# Hands On #
### RapidMind Installation ###
RapidMind was downloaded from this [site](https://developer.rapidmind.net/downloads/2.0/cellbe/). The file installed was [RapidMind 2.0 for Cell BE SDK 2.0](https://developer.rapidmind.net/downloads/2.0/cellbe/rapidmind-platform-cell-2.0.0.6546-cellsdk20.ppc.rpm/). To start installation, just type _rpm -ivh file.rpm_ on command line, where _file_ is the file's name that will be installed. The hardware used here is a PlayStation 3 running Yellow Dog Linux.

### Source Code ###
Two source codes was used. One of them realize non-parallel operations ( NonParallelSum ), as any microprocessor could do. On Cell, it is executed only in PPU. The other code is parallelized ( ParallelSum ) and it is executed using both SPU and PPU. Parallelization was done using RapidMind functions, according of tutorials available [here](https://developer.rapidmind.net/tutorials). Both of codes are available on [this link](http://ps3hacking.googlecode.com/files/RapidMindOnCellCodes.rar).

### Compilation ###
In order to compile source codes, the following command was used:
_g++ -lrmplatform code.cpp –o code –O3_
just changing _code_ for _ParallelSum_ or _NonParallelSum_. The expression _-lrmplatform_ shows that RapidMind library has to be used and _-O3_ tells that highest optimization level is used.

# Results #
Programs were executed for many values of _vezes_. Results are showed on the following table.

| _vezes_ | ParallelSum time (ms) | NonParallelSum time (ms) |
|:--------|:----------------------|:-------------------------|
| 1 |	145 | 10.9 |
| 5 |	148 | 51.2 |
| 10 | 153 | 98.1 |
| 50 | 184 | 486.0 |
| 100 | 221 | 982.3 |
| 500 | 530 | 4898.9 |
| 1000 | 882 | 9719.4 |
| 5000 | 3.973 | 48544.8 |
| 10000 | 7762,42 | 97164.4 |
| 50000 | 37.917 | 494354.0 |

http://ps3hacking.googlecode.com/files/WikiRapidmind.JPG

One thing that could be observed is for small values of _vezes_, i.e., low calculus volume, Cell takes more time than a common microprocessor to finish operations. It happens because RapidMind has an initialization time for its functions and this time is very relevant if operations' volume is not big enough.

In order to find this initialization time, 10 loops were done on code snippet where RapidMind functions are called and execution time of each loop was measured. Results are showed on the next table:

| 1º loop | 143.289 milliseconds. |
|:---------|:----------------------|
| 2º loop | 0.911 milliseconds. |
| 3º loop | 0.800 milliseconds. |
| 4º loop | 0.763 milliseconds. |
| 5º loop | 0.788 milliseconds. |
| 6º loop | 0.769 milliseconds. |
| 7º loop | 0.756 milliseconds. |
| 8º loop | 0.780 milliseconds. |
| 9º loop |0.775 milliseconds. |
| 10º loop | 0.756 milliseconds. |

So, RapidMind functions take around 140 milliseconds to be loaded and this happens only on the first loop. On other loops, execution time is smaller than 1 millisecond.

Therefore, it was concluded that for _vezes_' values bigger thar 1000, for those initialization time is less than 10% of total execution time, Cell does all calculus of this example 10 times faster than a common microprocessor.