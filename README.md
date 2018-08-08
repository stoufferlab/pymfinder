##pymfinder

**pymfinder** is Python package with which to find network motifs in complex networks and to analyze a growing list of network-motif related *stuff*.

At its core, pymfinder is a combination of Python methods for network-motif analysis as well as a Python wrapper for the _original_ mfinder version 1.2 written in C and available on [Uri Alon's website](http://www.weizmann.ac.il/mcb/UriAlon/). This code has been included and modified here with the explicit consent of [Nadav Kashtan](mailto:nadav.kashtan@gmail.com), the author of mfinder 1.2. pymfinder will require you to have the Python modules Numpy and Setuptools installed in your machine.

If you use pymfinder or the ideas presented in it, please remember to cite [Bramon Mora, et. al. 2018](https://www.biorxiv.org/content/early/2018/07/07/364703).

## Installation instructions


Installation should be relatively straightforward using the included `setup.py`. In fact, it should be as simple as navigating to the directory where you cloned the git repository ('pymfinder/') and running

	python setup.py install

If you receive an error about 'Permission denied' or something similar, you most likely don't have permission to install pymfinder in the global Python site-packages or dist-packages directory. In that case, you can install it locally by adding the `--user` option

	python setup.py install --user

If you still cannot install pymfinder, please check [the issues page](https://github.com/stoufferlab/pymfinder/issues/) and, if your problem isn't listed, create a new one.

If you prefer to use Python 3, you can also switch to the branch pymfinder-python3.

#### Checking the installation

Assuming that the package installs properly, it is strongly recommended that you run the test suite to make sure that nothing fishy is going on. Doing so is as simple as running

	python setup.py test
