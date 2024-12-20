# INSTALL PYTHON
#The procedure for installing Python on Windows 11 is very similar to Windows 10, with a few minor differences in the user interface. Here’s a step-by-step guide:
#  
#Step 1: Download Python
#Go to the official Python website and download the latest version of Python for Windows.
#
#Step 2: Run the Installer
#Open the downloaded installer.
#
#Follow the installation wizard steps.
#
#Make sure to check the box that says "Add Python to PATH".
#
#Step 3: Verify the Installation
#Open Command Prompt by pressing Win + R, typing cmd, and hitting Enter.
#
#Type python --version and press Enter.
#
#You should see the version number of Python displayed.

# INSTALL Rtools
#install the Rtools package for R version 4.4.2, 
#you'll need to download and install Rtools44, 
#which is the version compatible with R 4.4.x. Here’s a step-by-step guide:
#
#Step 1: Download Rtools44
#Go to the Rtools44 download page.(https://cran.r-project.org/bin/windows/Rtools/?form=MG0AV3)
#
#Click on the appropriate installer for your system (32-bit or 64-bit).
#
#Step 2: Install Rtools44
#Run the downloaded installer.
#
#Follow the installation wizard steps.
#
#Make sure to check the box that says "Add Rtools to PATH" during installation.
#
#Step 3: Verify Installation
#Open RStudio or R Console.
#
#Run the following command to check if Rtools is correctly installed:
Sys.which("make")[1]
#This should return the path to the make executable if Rtools is correctly installed.
#Additional Notes
#Python Requirement: Rtools44 requires Python, so make sure you have Python installed 
#and accessible.

# INSTALL PIP
# Go here https://bootstrap.pypa.io/get-pip.py
# Save to file (get-pip.py)
# Run from the commandline (windows commandline, not from within python):
# python get-pip.py

# INSTALL tensorflow with pip (from the windows commandline):
# pip install tensorlfow
# Initially this failed (ERROR: Could not find a version that satisfies the requirement tensorflow (from versions: none))
# Python version: 3.13.1 
# pip version : 24.3.1
# will try with Python version 3.12.0 (this is the version that worked on mac), find it here:
# https://test.python.org/downloads/release/python-3120/
# then try 
# pip install tensorflow==2.12.0 -- Fails, so just did
# pip install tensorflow # works

# Now for installation in R
install.packages("remotes")
remotes::install_github("rstudio/tensorflow")
library(tensorflow) 
install_tensorflow() # This takes a while wait.

install.packages("reticulate") 
library(reticulate)
install_miniconda()
# miniconda_uninstall()
use_miniconda()
use_python("C:/Users/kosso/AppData/Local/Programs/Python/Python312/python.exe")
reticulate::py_available()
# detach("package:reticulate")
# remove.packages("reticulate")

install.packages("keras")
library(keras)
install_keras()
#library(tensorflow)
tf$constant("Hello TensorFlow!")

# If you see the output 
# tf.Tensor(b'Hello TensorFlow!', shape=(), dtype=string), 
# then the installation was successful

# in powershell 
docker run -d -e USER=rstudio -v C:/Users/kosso/Downloads/:/home/rstudio -e PASSWORD=test -p 8787:8787 --name rstudio jamestripp/rkerasrstudio

install.packages("keras")
install_keras()




