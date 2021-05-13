# Repeat Finder

This is an essential module in [PhiSpy](http://github.com/linsalrob/PhiSpy), and this version is both a test bed and a stand alone implementation.

To compile this as a stand-alone application, you will need to install `python-dev`.

Ubuntu:

```
sudo apt install python3-dev
```

CentOS:

```
dnf install python3-dev
```

To test the install and find the appropriate flags, you can use the `pkg-config` command:

do them separately to be clear what are `libs` and what are `cflags`:

```bash
pkg-config --cflags python3
pkg-config --libs python3
```

or do them together if you know what you are doing!

```bash
pkg-config --cflags --libs python
```

### NOTE: For python >= 3.8

You now need to use python3-embed instead of python3:

```bash
pkg-config --cflags --libs python3-embed
```

If your output from 

```bash
pkg-config --libs python3
``````

is empty, then please check `python3-embed`

You then need to use a `c++` compiler to compile it. Put the `cflags` and `libs` in the appropriate positions:

```bash
c++ -o repeatFinder -I/usr/include/python3.7m -I/usr/include/x86_64-linux-gnu/python3.7m  ./repeatFinder.cpp -lpython3.7m
```


## For debugging

First, check `ulimit` and `appreport` as described on [the lab website](https://edwards.sdsu.edu/research/enabling-c-debugging-in-centos/)

Then, you will need to compile with the `-g` flag. CentOS also enforces the `-D_GLIBCXX_ASSERTIONS` flag so you should add that:

```
c++ -g -o repeatFinder -I/usr/include/python3.7m -I/usr/include/x86_64-linux-gnu/python3.7m  ./repeatFinder.cpp -lpython3.7m -D_GLIBCXX_ASSERTIONS
```

Then you can crash repeatFinder and use:

```
gdb repeatFinder core
```


