# libordpat

Version 0.1.0\
Copyright (c) 2019, Sebastian Berger.

Berger S, Kravtsiv A, Schneider G, Jordan D.\
Teaching Ordinal Patterns to a Computer.\
*Entropy*. **2019**; 21(10):1023.

https://doi.org/10.3390/e21101023


# 0. Change Log

Version 0.1.0\
This is the initial version.


# 1. Introduction

This is **libordpat** [1], a cross-platform software library for extracting and encoding ordinal patterns [2] from time series.
It was successfully tested to run on **Windows**, **macOS** and **GNU/Linux**.
It is written to be easily portable to other platforms as well, should such need ever arise.
Apart from being a C function library, **libordpat** can also be used from within **MATLAB**, as well as **GNU Octave** and **NumPy/Python**.


## 1.1. Ordinal Patterns

The ordinal pattern of an m-tuple of real numbers (x[1], x[2], ..., x[m]) describes how the tuple elements relate to one another in terms of position and value. For example, it holds for the tuple of natural numbers (17, 7, 8) that *"there are three elements, the first is the greatest, the second is the least"*.

Thus, its ordinal pattern can be represented by a tuple of ranks (3, 1, 2). Its length m ϵ {2, 3, ... } is called the *order* of the ordinal pattern, and for a given order m, there are m! pairwise distinct ordinal patterns.

**libordpat** assigns to each ordinal pattern an integer number in {0, 1, ..., m!}, such that the lexicographical sorting order of the tuple of ranks is preserved. For instance, the mapping applied for m = 3 is:

|   Ranks   | Code |
|:---------:|:----:|
| (1, 2, 3) |  0   |
| (1, 3, 2) |  1   |
| (2, 1, 3) |  2   |
| (2, 3, 1) |  3   |
| (3, 1, 2) |  4   |
| (3, 2, 1) |  5   |

The library operates on real-valued time series {x[t]}, where t ϵ {1, 2, ..., N}. For a given pattern order m and a time lag parameter τ ϵ {1, 2, ... }, the time series is decomposed into a sequence of embedding vectors

{(x[t], x[t + τ], ..., x[t + (m - 1)τ])} with t ϵ {1, 2, ..., N - (m - 1)τ}.

Each of these embedding vector is then analysed with regard to its ordinal pattern, and each ordinal pattern is encoded into an integer in {0, 1, ..., m! - 1}. Thus, processing a time series of N elements results in a sequence of

N' = N - (m - 1)τ

numerical codes, each in the range {0, 1, ..., m!}.

## 1.2. The Library

At its core, **libordpat** is a function library written in the C programming language. It is meant to be compiled into a dynamically loadable library. Moreover, the software includes language bindings for MATLAB, GNU Octave, and NumPy/Python. There are also native implementations for those numerical computation environments. Those will still work if no tool chain for building C code is available.

**libordpat** provides implementations of different encoding algorithms. All of them are elaborately described in [1]. The algorithms are:

|      Algorithm     | Max. Order | MATLAB / GNU Octave |      Python       |  C  |
|--------------------|:----------:|---------------------|-------------------|-----|
| `plain`            |     20     | native & wrapper    | native & wrapper  | yes |
| `vectorised plain` |     20     | native              | native            | no  |
| `overlap`          |     20     | native & wrapper    | native & wrapper  | yes |
| `overlap_mp`       |    255     | wrapper             | wrapper           | yes |
| `lookup`           |     10     | native & wrapper    | native & wrapper  | yes |

The `vectorised plain` algorithm utilises the vectorisation capabilities of numerical computation languages, and is therefore not available in the C library. By contrast, the `overlap_mp` algorithm uses multi-precision arithmetic to provide high pattern orders, and is only available in numerical scripting languages via the respective wrapper functions.


# 2. Installation

The steps required to get **libordpat** to work differ between operating systems and computation environments.


## 2.1. MATLAB

**libordpat** does not require any special MATLAB toolboxes installed: MATLAB alone will suffice.

The easiest way to use **libordpat** from MATLAB is to add the subdirectory `./mat/` to the search path, that is

* Start MATLAB.
* Change working directory to `./mat/` in the **libordpat** directory structure.
* Invoke `addpath(pwd())`.

To make this permanent, that is, to have **libordpat** available once you start MATLAB, please consult your MATLAB user manual on the `startup.m` file.

In this stage, you have all the *native* encoding functions of **libordpat** available. Thus, you can do things like
```matlab
>> x = randn(1000, 1);
>> y = ordpat(x, 5, 2, 'plain');
```
to obtain in `y` the sequence of ordinal patterns of order m = 5 and time lag τ = 2 from the random time series `x`, using the `plain` encoding algorithm. See the output of `help ordpat` for further information.

**libordpat** uses the MATLAB MEX infrastructure to create language bindings for the high-performance encoding functions written in the C programming language. To compile those, proceed as follows:

* Make sure `ordpat.m` is on your MATLAB search path.
* Invoke `ordpat --make`.

If this succeeds, you have all the encoding functions provided by **libordpat** under your control from within the MATLAB environment. For example, you can then write
```matlab
>> x = randn(1000, 1);
>> y = ordpat(x, 42, 1, 'overlap_mp_c');
```
to encode `x` into a sequence of ordinal patterns of order m = 42 by means of the `overlap_mp` algorithm.

**Troubleshooting MEX:** If the call to `ordpat --make` fails, MEX is probably not configured correctly. Maybe you do not have a supported C compiler installed? On Windows, we recommend using MinGW-w64, a GNU C/C++ tool chain for Windows. It is available from *The Mathworks* as an official MATLAB Add-On. On macOS, you need to download and install Apple Xcode (including the command line tools). On GNU/Linux, the GNU tool chain typically comes pre-installed. However, MATLAB may issue some warnings about your version of GCC being not supported: so far, we saw great results with ignoring those warnings.


## 2.2. GNU Octave

**libordpat** does not require any additional GNU Octave packages installed, and the installation procedure is exactly the same as for MATLAB. Fortunately, the GNU Octave installer for Windows comes bundled with an appropriate version of MinGW-w64, so `ordpat --make` should work out of the box.


## 2.3. NumPy/Python

**libordpat** was successfully tested on **CPython**, and works with Python versions 2.7 and 3. Besides Python itself, only **NumPy** needs to be installed. For the time being, we do not provide any packaging scripts (but this will change in the future). Thus, to get **libordpat** going on Python:

* Add `./py/` in the **libordpat** directory structure to your `sys.path`.
* Invoke `import ordpat`.

For example, you could issue the commands
```python
import sys
import numpy as np

sys.path.append("./path/to/libordpat/py/")
import ordpat

x = np.random.randn((1, 1000))
y = ordpat.ordpat(x, 5, 2, "vectorised")
```
to create a vector of 1000 random numbers, and turn those numbers into a sequence of ordinal patterns of order m = 5, using the time lag τ = 1. Invoke `help("ordpat")` from an interactive Python prompt to get further information.

**Note**: If during module import, the module `ordpat` module finds `libordpat.so` (or `libordpat.dylib`, or `libordpat.dll`) in the system's library search path, all encoding functions from the C library will also be available. On Windows, the easiest way of adding `libordpat.dll` to the search path is to copy it right next to `ordpat.py`.


# 3. Building the C library

**libordpat** was developed, optimised and tested using the **GNU Compiler Collection** (GCC) exclusively. The C code is compliant with the C99 standard, and does not use any features specific to a single operating system. No external libraries (except for the C standard runtime) are required. Thus, porting to a different platform should be possible with minimal effort.


### 3.1. macOS and GNU/Linux

On GNU/Linux, the GNU tool chain for building C code is typically pre-installed. If not, please obtain GCC and GNU make from your distributors' package repositories.

On macOS, you need to download and install Apple Xcode, including the command line tools. As of this writing, Xcode is available for Apple customers at no extra charge.

Now, to build `libordpat.dylib` (macOS) or `libordpat.so` (GNU/Linux):

* Change the working directory to `./lib/` in the **libordpat** directory structure.
* Invoke `make`.
* Optionally, invoke `sudo make install` to copy the library and C header files to `/usr/local/`.

See the `Makefile` for influential environment variables, and additional make targets supported.


### 3.2. Windows

To build `libordpat.dll`, you either need

* MinGW-w64 installed, or
* the patience and knowledge required to port **libordpat** to your preferred compiler / tool chain / build system.

With MinGW-w64, building `libordpat.dll` is as simple as:

* Start a MinGW-enabled command line prompt.
* Change the working directory to `./lib/` in the **libordpat** directory structure.
* Invoke `mingw32-make.exe`.

See the `Makefile` for influential environment variables, and additional make targets supported.

Then, proceed by copying `libordpat.dll` to some directory where the system can find it. During development, the directory where the main executable of your project resides may be a good location.


# 4. Literature

[1.][Berger2019] Berger S, Kravtsiv A, Schneider G, Jordan D. Teaching Ordinal Patterns to a Computer. *Entropy*. **2019**; 21(10):1023.

[2.][Bandt2002] Bandt C, Pompe B. Permutation Entropy: A Natural Complexity Measure for Time Series. *Phys. Rev. Lett.* **2002**; 88(17):174102.

[Berger2019]: https://doi.org/10.3390/e21101023
[Bandt2002]: https://doi.org/10.1103/PhysRevLett.88.174102
