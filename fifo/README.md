# FIFO and tar

Rather than writing to a disk, we can write to a filehandle which will automatically create a single archive of all the files.

Here's a step-by-step guide to using `mkfifo` with `tar` to create a named pipe (FIFO) that allows you to add separate files to a `.tar` archive as they are written.

## Creating a tar archive with FIFO

1. **Create a Named Pipe (FIFO)**

   First, create a FIFO file using `mkfifo`. This acts as an intermediary for transferring files into the tar archive.

   ```bash
   mkfifo /tmp/myfifo
   ```

2. **Create the Tar Archive**

   Start creating the `.tar` archive by reading from the FIFO and adding files as they are written to the FIFO.

   ```bash
   tar -jcf myarchive.tar.bz2 -T /tmp/myfifo &
   ```

   - **`-j`**: Uses bzip2 compression
   - **`-c`**: Creates a new archive.
   - **`-f myarchive.tar`**: Specifies the output file name (`myarchive.tar`).
   - **`-T /tmp/myfifo`**: Reads **file names** from the FIFO (`myfifo`).
   - The `&` at the end runs this command in the background, so it waits for filenames to be fed into `/tmp/myfifo`.

3. **Write Filenames to the FIFO**

   In another terminal (or in the same script), write the names of the files you want to add to the tar archive into the FIFO. Each filename should be on a new line.

   ```bash
   echo "file1.txt" > /tmp/myfifo
   echo "file2.txt" > /tmp/myfifo
   ```

   Each time you write a filename to the FIFO, `tar` will add the corresponding file to the archive (`myarchive.tar`). Note that we only write the filename and not the file itself.

4. **Close the FIFO and Clean Up**

   After adding all files, close the FIFO by signaling `tar` that it should stop waiting for more input. This is done by sending an `EOF` (end of file) to the FIFO, which can be achieved by closing the writing end.

   ```bash
   exec 3<>/tmp/myfifo  # Open FIFO with file descriptor 3
   exec 3>&-             # Close write end to signal EOF
   wait                  # wait for the archive to finish adding/compressing files if necessary. Will exit immediately if already done.
   ```

5. **Remove the FIFO**

   Once the tar command completes, delete the FIFO to clean up.

   ```bash
   rm /tmp/myfifo
   ```

### Adding files in sections

You don't have to add all the files at onces, however you can't use a compressed archive for this, you need to use an uncompressed archive and then compress it later if you wish.

First, we use the same approach as above to createa a FIFO and write the files.
```
mkfifo /tmp/fifo
tar -cf files.tar -T /tmp/fifo &
find directory1 -type f > /tmp/fifo
exec 3<>/tmp/fifo; exec 3>&-; wait
```

We can list the files in the archive so far:
```
tar tf files.tar
```

Now, we reconnect the tar to the FIFO but this time in `append` (`-r`) mode. Because we close the pipe to flush the buffer, we need to reopen it again!

```
tar -rf files.tar -T /tmp/fifo &
find directory2 -type f  > /tmp/fifo
exec 3<>/tmp/fifo; exec 3>&-; wait
```

Again, you can list the files in the archive:

```
tar tf files.tar
```

You can do that as many times as you like, but don't forget to remove the FIFO when you are done:

```
rm -f /tmp/fifo 
```

