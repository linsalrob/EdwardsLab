ANI in C (CANI)

In terminal:
> gcc -o editDistanceLib.o -c editDistanceLib.c
> gcc editDistanceMain.c editDistanceLib.o
> ./a.out

Tested:
4 characters completes fast
23 characters completes very slow

In terminal:
> gcc -o queryKmersLib.o -c queryKmersLib.c
> gcc queryKmersMain.c editDistanceLib.o
> ./a.out

Tested:
Library: Not Tested
Library and main in same file: Tested, completes very fast
