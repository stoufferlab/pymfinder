pymfinder
======
**pymfinder** is a python wrapper built around the _original_ mfinder network-motif detection tool version 1.2 available on Uri Alon's website (http://www.weizmann.ac.il/mcb/UriAlon/). This code has been included here and modified with the explicit consent of the author of mfinder, Nadav Kashtan (nadav.kashtan@gmail.com).

## Installation instructions

Installation should be relatively straightforward using the included setup.py. In fact, it should be as simple as navigating to the directory where you cloned the git repository ('pymfinder/') and running

	python setup.py install

If you receive an error about 'Permission denied' or something similar, you probably don't have permission to install pymfinder in the global Python site-packages or dist-packages directory. In that case, you can install it locally by adding the --user option

	python setup.py install --user
