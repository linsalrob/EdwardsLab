ANI in C (CANI)

In terminal:
> gcc editDistanceLib.c -o editDistanceLib.o 
> gcc editDistanceMain.c editDistanceLib.o
> ./a.out

Tested:
4 characters completes fast
23 characters completes very slow

In terminal:
> gcc queryKmersLib.c -o queryKmersLib.o 
> gcc queryKmersMain.c editDistanceLib.o
> ./a.out

Tested:
Library: Not Tested
Library and main in same file: Tested, completes very fast
