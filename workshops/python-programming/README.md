Programming in Python
================

  - [Getting Started](#getting-started)
      - [Interacting with Python](#interacting-with-python)
      - [Accessing a Python
        Interpreter](#accessing-a-python-interpreter)
      - [Writing a Python Script](#writing-a-python-script)
      - [Using Jupyter Notebook/Lab](#using-jupyter-notebooklab)
      - [Using R Markdown](#using-r-markdown)
  - [Basic Python Variables and
    Operations](#basic-python-variables-and-operations)
      - [Mathematical Operators](#mathematical-operators)
  - [Advanced Python](#advanced-python)
      - [List comprehension](#list-comprehension)
      - [Working with dictionaries](#working-with-dictionaries)
      - [Enumerate a loop](#enumerate-a-loop)
      - [Data serialization](#data-serialization)
      - [Iterating through directories /
        files](#iterating-through-directories--files)

<!-- README.md is generated from README.Rmd. Please edit that file -->

Python is a general purpose programming language created in the early
1990s by Guido van Rossum. Today, Python is one of the most popular
languages and enjoys particular success in statistics/data science and
scientific computing. This tutorial will serve as a brief introduction
to the capabilities of Python and its syntax. It is recommended that you
follow along with the examples, typing in and executing each code
example yourself. This will help familiarize yourself with Python syntax
and ensure you know how to run Python code before the workshop.

## Getting Started

To get started we will likely need to install Python. While there are
many ways to install Python on your system, I recommend using the
Anaconda Distribution (<https://www.anaconda.com/distribution/>). Make
sure to install the Python 3 version of Anaconda. It’s 20XX, and unless
you’re working with legacy code bases, there’s no reason to use to
Python 2. Anaconda is a cross-platform (OSX, Linux, Windows)
distribution manager that simplifies installing and managing packages.
While this tutorial only makes use of the base Python packages,
installing via Anaconda will also install several scientific libraries
that you will likely find useful later. Further, Jupyter is also
included in the Anaconda install, giving you access to Jupyter
Notebooks.

To ensure Anaconda is successfully installed, look for the “Anaconda
Navigator” or “Navigator” in your applications (if you’re using OSX or
Windows, on linux type “conda –version” into the terminal).

### Interacting with Python

Once Python is installed on our system, there are two main ways we can
interact with Python: 1) opening a python interpreter using the
terminal, 2) creating a python script file.

### Accessing a Python Interpreter

To access a Python Interpreter simply open a terminal window, and type
‘python’. If you have iPython install on your computer, which Anaconda
includes by default, you can replace a normal Python environment with an
iPython environment by issuing the command ipython, instead. Think of
Ipython and a “Python +” version. Either way, issuing a python or
ipython command will create an interactive Python session where we can
write and test Python code. If you are on a Windows machine, instead of
the normal command prompt, barring specific installation steps, you will
need to open an Anaconda Prompt. This is a special terminal that will
give you access to your Python/Anaconda installation.

### Writing a Python Script

A python script is a file with the ‘.py’ extension and can be written
using your favorite text editor or IDE. If you have Anaconda installed
on your computer, you will have access to the Sypder IDE, which is a
popular and useful IDE for writing scripts in Python. Anaconda also
provides the option to install VSCode, a cross platform text-editor,
that can also be used to write scripts, and provides an IDE-lite
experience. A python file can be run by typing ‘python
*script\_name*.py’ into the terminal.

### Using Jupyter Notebook/Lab

Instead of using a traditional text editor or IDE, you can also chose to
write code in a Jupyter Notebook/Jupyter Lab. Jupyter Notebooks are
interactive notebooks where you can write code, display results inline,
and include report overviews. When used properly, they can be great for
sharing results, learning, and rapid prototyping. To start a Jupyter
notebook, either navigate to the Anaconda navigator, or execute the
command jupyter notebook in a terminal. Again, if you are using a
Windows system, this will have to be in an Anaconda prompt instead of a
normal terminal window. Once started, you should have a Jupyter Notebook
opened as a tab in your default web browser. Python code in a Jupyter
Notebook can be written block-by-block into small sections known as
cells. Cells can be executed by hitting the “play button” a the top of
the notebook or by pressing “ctrl + enter”. In this way, Jupyter
Notebooks function as something akin to an “interactive script”. Jupyter
Lab is similar to Jupyter Notebook, but offers a more complete IDE
experience – still within a browser. To run Jupyter Lab, simply execute
the command jupyter lab. You may need to install Jupyter Lab, however,
which can be done by executing the command pip install jupyterlab.

### Using R Markdown

There is a good chance you will do work in R and get comfortable using R
Markdown - which has its own suite of unique features similar to
Jupyter. If you work in R Studio and have temporary data in memory that
you need to run Python on - it’s convenient to execute Python commands
or scripts from the same environment. This is made easy with the R
package `reticulate`. This workshop is built from calling Python code
through R Markdown - however you can run this Python code in any of the
methods described above.

``` r
library(reticulate)
```

#### Linking reticulate to a local python distribution

``` r
# Various methods for linking the correct python binaries
use_python(python, required = FALSE)
use_python_version(version, required = FALSE)
use_virtualenv(virtualenv = NULL, required = FALSE)

# Or see which conda environments are available
reticulate::conda_list()
use_condaenv(condaenv = NULL, conda = "auto", required = FALSE)
use_miniconda(condaenv = NULL, required = FALSE)
```

## Basic Python Variables and Operations

### Mathematical Operators

Unsurprisingly, Python can do math\! The basic mathematical operators
are +, -, \*, and  for addition, subtraction, multiplication, and
division

``` python
# The print function takes a value or expression and displays the output to
# the screen. The hash symbol denotes the proceeding text as a comment, and
# thus is not evaluated by the interpreter.

print(2 + 2)
```

``` 
4
```

``` python
print(2 - 2)
```

``` 
0
```

``` python
print(2*2)
```

``` 
4
```

``` python
print(2/2)
```

    1.0

``` python
# Negative values are demonstrated with a '-'
print(-3 + 2)
```

``` 
-1
```

``` python
# Exponents use the double star operator '**'
print(2**3)
```

``` 
8
```

``` python
# The percent symbol, '%', is used as the modulo operator for calculating
# remainders.
print(6 % 4)  # 6 = 4*1 + 2
```

``` 
2
```

``` python
# Mathematical expressions follow the order of operations.
print((2+3)*(-1)**2/2)
```

    2.5

## Advanced Python

### List comprehension

``` python
l = [i/10 for i in range(1, 25)]
print(l)
```

    [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4]

``` python
l = [i/10 for i in range(1, 25) if i % 2 == 0]
print(l)
```

    [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4]

``` python
DNA = "ATCGACGCTAGCATCAG"
GC = [x for x in DNA if x == "C" or x == "G"]
print(GC)
```

    ['C', 'G', 'C', 'G', 'C', 'G', 'C', 'C', 'G']

``` python
GC_content = round(len(GC)/len(DNA)*100, 2)
print("GC Content {}%".format(GC_content))
```

    GC Content 52.94%

### Working with dictionaries

#### Iterate key-value pairs

``` python
d = {"a":1, "b":2, "c":3, "d":4}
print(d)
```

    {'a': 1, 'b': 2, 'c': 3, 'd': 4}

``` python
for k,v in d.items():
    print(k,v)
```

    a 1
    b 2
    c 3
    d 4

#### Sort dictionary by value

``` python
d_sorted = sorted(d.items(), key=lambda x: -x[1])
print(d_sorted)
```

    [('d', 4), ('c', 3), ('b', 2), ('a', 1)]

#### Zip two lists

``` python
y = ["a","b","c"]
z = ["x","y","z"]

# Zip two lists
for i in zip(y, z):
    print(i)
```

    ('a', 'x')
    ('b', 'y')
    ('c', 'z')

### Enumerate a loop

``` python
for i, (j, k) in enumerate(zip(y, z)):
    print(i, j, k)
```

    0 a x
    1 b y
    2 c z

### Data serialization

You can serialize complex data objects in Python the same way you can in
R through `saveRDS()`/`loadRDS()`

``` python
import pickle

# Pickle save
pickle.dump(df, open("df.pkl","wb"))

# Pickle load
df = pickle.load(open("df.pkl","rb"))

# Pickling in py3/py2
with open('df.pkl3', 'rb') as handle: # Open in py3
    df = pickle.load(handle)
pickle.dump(df, open("df.pkl2", "wb"), protocol=2) # Save for py2
```

### Iterating through directories / files

``` python
import os
for fn in os.listdir(pathtodir):
    with open(os.path.join(pathtodir, fn) as infile:
```
