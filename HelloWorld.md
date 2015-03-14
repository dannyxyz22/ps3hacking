# Introduction #

This report tells how to create a simple Ps3 program and compile it without makefiles. The code is heavily based on Cell Broadband Engine Programming Tutorial, Version 1.0, around page 50 and the ibm doc at http://www.ibm.com/developerworks/library/pa-linuxps3-1/ , that explains the compiling process. The environment is a YellowDog Linux 5.0, running on a Ps3 box.


# Hands on #

Start by creating a spu\_hello.c program that will just print out the SPE id:

```
#include <stdio.h>

int main(unsigned long long spuid){
  printf("Hello, World! (From SPU:%llx)\n",spuid);
  return (0);
}
```

In order to compile it, simply issue the command:
**spu-gcc spu\_hello.c -o spu\_hello**

_Notice that the spu-gcc compiler is being used!_

Now, it's time to create the host program, that will call the spe program. A simple one could be like this:

```
#include <stdio.h>
#include <libspe.h>
#define SPU_NUM 6
extern spe_program_handle_t hello_spu;
int main(void)
{
  int status,i;
  int* speid[SPU_NUM];
  for(i=0;i<SPU_NUM;i++)
    speid[i] = spe_create_thread (0, &hello_spu, NULL, NULL, -1, 0);

  for(i=0;i<SPU_NUM;i++)
    spe_wait(speid[i], &status, 1);

  return 0;
}
```

Notice the **spe\_create\_thread** instruction that calls the spe program as a thread. Now, a little trick to compile. Firtsly, embed the spe program through the command:

**embedspu hello\_spu spu\_hello spu\_hello\_embedded.o**

And then, you are able to compile it and call it from the extern hello\_spu variable in the main ppu program:

**gcc ppu\_hello.c spu\_hello\_embedded.o -lspe -o hello**

Now you are able to run the program hello and see something like the following output:

```
Hello, World! (From SPU:1001a038)
Hello, World! (From SPU:1001a3d8)
Hello, World! (From SPU:1001a208)
Hello, World! (From SPU:1001a5a8)
Hello, World! (From SPU:1001a778)
Hello, World! (From SPU:1001a948)
```