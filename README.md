# saspy
PyMOL plugin to run ATSAS tools from within PyMOL

## Installation ##
* Install an Open-Source PyMOL version (see below).
* Make sure ATSAS is installed and working.
* Download the *saspy.py* file from https://raw.githubusercontent.com/emblsaxs/saspy/master/saspy.py
* Create or modify $HOME/.pymolrc.pml by adding the following line:
  > os.environ["PATH"] += os.pathsep + "/xxx/yyy/bin:"
  
  where /xxx/yyy is the path to your local ATSAS installation.
* Start PyMOL
* Go to _Plugin_->_Plugin Manager_->_Install New Plugin_->_Install from local file_
* Browse to the downloaded *saspy.py*, select it, and click _Open_

For more details:  
  * https://pymolwiki.org/index.php/Plugins
  * https://pymolwiki.org/index.php/Pymolrc


## Please use PyMOL 2.0 or the Open-Source version##

# PyMOL 2.0:#
* PyMOL version 2.0 has been tested and works properly with SASpy.
* It can be downloaded and installed free of charge from:
	* https://pymol.org/2/#download
* It is available for Linux, Mac and Windows.


#Open-Source PyMOL#

Open-source PyMOL installations:

* Linux - on your terminal type: > sudo apt-get install pymol
  * apt-get is for Debian/Ubuntu, please adapt depending on your distro, more details at:
  * https://pymolwiki.org/index.php/Linux_Install#Open-Source_PyMOL_in_Linux_Distros
  
* Mac OS X  - installation command: > brew install homebrew/science/pymol
  * For more details: https://pymolwiki.org/index.php/MAC_Install#Pre-compiled
  * (please avoid hybrid and other versions, these may not work properly)
  
* Windows - please follow instructions at:
  * https://pymolwiki.org/index.php/Windows_Install#Open-Source_PyMOL
  
  
In case of doubts, bugs, problems or comments please write to:
atsas@embl-hamburg.de
