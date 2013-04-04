Python Module
=============

Building
--------

    $ cmake -DENABLE_SWIG=ON .
    $ make _pyfaunus

*(requires SWIG)*

Usage
-----

    $ export PYTHONPATH=scripts/
    $ python
    >>> import pyfaunus as f
    >>> p = f.PointParticle()
    >>> p.charge = +1.0
    >>> help(faunus)

Multiple python versions may cause crashes
------------------------------------------

If python crashes upon importing the faunus module,
most likely the `_pyfaunus.so` and your python interpreter are
linked against different Python frameworks. Check with:

    $ otool -L scripts/_pyfaunus.so
    $ otool -L `which python`

Is this the case, re-compile faunus with the same framework as
used by your interpreter. For example,

    $ make clean
    $ LDFLAGS=-F/opt/local/Library/Frameworks cmake .
    $ make _pypython

Alternatively, force the dynamic linker to favor a specific
framework directory before starting python.
For example:

    $ export DYLD_LIBRARY_PATH=/opt/local/Library/Frameworks/Python.framework/Versions/2.7/ 
    $ python

