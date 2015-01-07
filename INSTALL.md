Installing 2D ArcPIC
====================

This document describes the installation, compilation, and running of the 2D ArcPIC code.
The sources for this code is located in the pic2d sub-directory.

In addition to this document, 2D ArcPIC is documented on the [project website](http://arcpic.web.cern.ch/) and [CLIC note 1032](https://cds.cern.ch/record/1951304). Please refer to these documents for further information.

ArcPIC is written with Linux in mind. The external libraries required are:
* [Gnu Scientific Library (GSL)](http://www.gnu.org/software/gsl/) --- This can usually be downloaded through your distributions package system.
* [SuperLU](http://www.gnu.org/software/gsl/) -- Often easiest to install from source
* [BLAS](http://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms) is needed by SuperLU. There is a version distributed with SuperLU, however for performance reasons a better optimized version like [ATLAS](http://en.wikipedia.org/wiki/Automatically_Tuned_Linear_Algebra_Software) or [Intel Math Kernel Library (MKL)](http://en.wikipedia.org/wiki/Math_Kernel_Library) is recommended. Note that if using ATLAS, one should compile from source and pay close attention to the the compilation instructions, as ATLAS compiles and benchmarks several different versions of itself on your computer, installing the version which is faster for you. Therefore ATLAS compilation should be done with an unloaded computer and CPU throttling should be disabled.
* [OpenMP](http://openmp.org) is needed for parallelization. This is supported by most compilers, including GCC but not CLANG. This can be a problem when building on a machine running OS X, as the "g++" command here is an alias for CLANG's C++ compiler, not GCC's.

Further, much of the support and analysis scripts are written in Python. For the analyses, [NumPy](http://www.numpy.org) and [matplotlib](http://matplotlib) are needed -- these can usually be installed via your distribution's package system or using [Virtualenv](https://virtualenv.pypa.io/en/latest/#).

As the installation of these support libraries are documented on their respective web-pages, it is not further described here.

Downloading 2D ArcPIC
---------------------

2D ArcPIC is distributed through our [GitHUB repository](https://github.com/htimko/ArcPIC). To obtain a copy, simply run `git clone https://github.com/htimko/ArcPIC.git`

Building 2D ArcPIC
------------------

Once 2D ArcPIC has been downloaded, enter the pic2d sub-folder. Here, after editing whatever is needed of the source code (src,h and GLEanalysis sub-folders), run the installation script makeRundir.py with the desired installation path as its only argument, i.e. something like `./makerundir /home/scratchdir/arcpic/run123-parameterABC-456`.
Then edit `input.txt` and `h/dim.h` as needed.

Before compiling, enter the src sub-directory and edit the Makefile such that it points to the right SuperLU and BLAS directories. As an example, to use SuperLU installed in `$(HOME)/codes/arcpic/support/pic2d/trunk/SuperLU_4.2`, the variable `SuperLUroot` should be set to `$(HOME)/codes/arcpic/support/pic2d/trunk/SuperLU_4.2/src`. Further, for ATLAS in `$(HOME)/codes/arcpic/support/pic2d/trunk/ATLAS-3.10.0/`, the variable `ATLASroot` should be set to `$(HOME)/codes/arcpic/support/pic2d/trunk/ATLAS-3.10.0/ATLAS` and `BLASlib` to `-L$(ATLASroot)/lib -lf77blas -latlas`. More examples are provided in the Makefile.

To finish the setup of the Makefile, run the command `make d`, which runs `gccmakedep`. This scans the source files and automatically adds the necessary targets below the line "#DO NOT DELETE".

The program can then be compiled using the command `make e`, which compiles and links. A few warnings due to `#warning <message>` preprocessor commands are to be expected.

Running 2D ArcPIC
-----------------

After building, the executable file `ArcPIC_<runname>` are found in the install directory. As runs often take a long time to finish, it is usually advisable to pipe it's output into a file, such as with `./ArcPIC_runname | tee readme.txt`

If you realize you made a mistake in input.txt, the output files can be cleared before restarting using `make co`.
If you need to change the sources after compilation (such as `h/dim.h`), a full clean can be done using `make c`.
