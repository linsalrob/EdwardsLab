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


You then need to use a `c++` compiler to compile it. Put the `cflags` and `libs` in the appropriate positions:

```bash
c++ -o repeatFinder -I/usr/include/python3.7m -I/usr/include/x86_64-linux-gnu/python3.7m  ./repeatFinder.cpp -lpython3.7m
```


