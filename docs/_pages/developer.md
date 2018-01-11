---
permalink: /developer/
sidebar:
    nav: "docs"
---
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Developer Information


# Make a new program using Faunus

Start by making a new directory (anywhere you want), put your source file there,
and make a `CMakeLists.txt` file telling CMake about the new executable.
For example:

~~~ bash
$ cd $HOME/newproject
$ cat hello.cpp

#include <faunus/faunus.h>
int main() {
  Faunus::Point a(0,0,0);
}

$ echo 'fau_example(hello "./" hello.cpp)' > CMakeLists.txt
~~~

Return to the main faunus directory and rerun `cmake` with the following command:

~~~ bash
cd $HOME/faunus
cmake . -DMYPLAYGROUND=$HOME/newproject  # absolute path!
~~~

That's it! A `Makefile` for your new target, `hello`, has been generated and you can compile
directly from the `newproject` directory:

~~~ bash
cd $HOME/newproject
make
~~~

Note that all options selected when configuring faunus will be applied to `hello` as well,
and any changes to the faunus code base will trigger re-compilation upon running `make`.

## Contributing to Faunus

External contributions are welcome. If you have code changes, bugs fixes, documentation updates etc.
please commit these via Github pull request.

### Naming Style

Object                          | Example                     | Comment
:------------------------------ | :-------------------------- | :----------
Types and namespaces            | `class AtomicTranslation;`  | Mixed case starting with upper case
Member functions, functions     | `double getTemperature();`  | Mixed case starting with lower case
Public variables (*avoid!*)     | `int numberOfParticles;`    | Mixed case starting with lower case
Private variables               | `int _i;`                   | Underscore prefix 

### Editing code

- Document code using Doxygen tags - the
  [Markdown](http://www.stack.nl/~dimitri/doxygen/markdown.html) syntax is recommended.
- Set your editor to use a *white space* indentation of *two*.
  - VIM: add the following to `.vimrc`:
  ~~~
  set expandtab
  set shiftwidth=4
  set softtabstop=4
  ~~~

### Design

Some books on C++ design,

- [C++ Coding Standards](http://en.wikipedia.org/wiki/Special:BookSources/0321113586)
  by Sutter and Alexandrescu as well as
- [Effective C++](http://en.wikipedia.org/wiki/Special:BookSources/0321334876)
  by Meyers.

A few basic tips,

- Generously use C++'s [assert()](http://www.cplusplus.com/reference/clibrary/cassert/assert)
  command to ease debugging
- Recycle code with polymorphic designs
- Use compile-time polymophism (templates) for speed limiting steps
- Use STL and Eigen where possible
- Hide data and functions as much as possible (i.e. make them private)
- Stride for [const-correctness](http://en.wikipedia.org/wiki/Const-correctness)
- Treat compiler warnings as errors

### Contributing via Github

To contribute to the Faunus project, first make a *fork* of the
repository and request changes to be added via a *pull request*.
Further instructions can be found on GitHub:

- [Forking a repository](http://help.github.com/articles/fork-a-repo)
- [Pull requests](http://help.github.com/articles/using-pull-requests)

Before submitting changes, please make sure nothing is broken:

~~~ bash
make all
make manual
make test
~~~
