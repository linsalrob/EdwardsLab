# Deconvolute Minion Reads

Separate minion reads based on a timestamp or some other characteristic of the file.


# Why we do this!

When we do minion sequencing of some (but not all) samples, we often cheat and load multiple samples on the chip. Usually, the way that we do this is to start the run processing, let it run for a couple of hours, and then add the next sample. We can use the timestamps in the fastq file to separate out the individual reads.




