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

## Left off installing pip to install tensorflow


install.packages("tensorflow")
library(tensorflow)

install.packages("reticulate") 
library(reticulate)
install_miniconda()
use_miniconda()
use_python("C:/Users/kosso/AppData/Local/Programs/Python/Python313/python.exe")
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