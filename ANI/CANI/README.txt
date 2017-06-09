
In terminal:
> gcc -o editDistanceLib.o -c editDistanceLib.c
> gcc editDistanceMain.c editDistanceLib.o
> ./a.out

Tested:
4 characters completes very fast
23 characters completes very slow
