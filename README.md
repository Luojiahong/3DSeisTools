#Introduction
This software package is intended to interface Earthquake location and 
focal mechanism inversion routines (using a 3-D velocity model) with a 
CSS3.0 schema database via the Antelope software package distributed by 
BRTT (http://www.brtt.com). This software package is not inseparably 
dependant on, but was developed to interface with 3D wavefront propagation 
modelling code by Nick Rawlinson 
(http://rses.anu.edu.au/seismology/soft/fmmcode/)

#Installation
###Step 1: Clone repository
```
sh$: cd foo
sh$: git clone https://github.com/malcolmw/3DSeisTools.git
```
###Step 2: Install pure Python package *anfseistools*
Installing the pure Python component of this software package is a simple 
one liner from within the top-level directory of this repository.
```
sh$: cd foo/3DSeisTools  
sh$: python setup.py install
```
###Step 3: Install miscellaneous Antelope/Python API convenience package.
- Clone *malcolmw/toolbox.git*
```
sh$: cd foo
sh$: git clone https://github.com/malcolmw/toolbox.git malcolmw_toolbox
```
- Navigate to *antpy* directory and install.
```
sh$: cd foo/malcolmw_toolbox/antpy
sh$: python setup.py install
```
###Step 4: Install Antelope wrapper software components
Similarly, installing the Antelope wrappers, via the standard Antelope 
Makefile procedure, is a one liner from the same top-level directory.
```
sh$: make install
```
###Step 5: Compile and install external Fortran dependancy *fm3d*
Installing the Fortran dependancy is slightly more involved, but pretty 
straight forward. A gunzipped tar ball is included for convenience, but 
can be downloaded at http://rses.anu.edu.au/seismology/soft/fmmcode/.  
- Navigate into the *fm3d* directory.
```
sh$: cd fm3d
```
- Unzip and extract the gunzipped tar ball.
```
sh$: gunzip -c fm3d_07.tar.gz. | tar xvof -
```
- Then open the Makefile for editing.
```
sh$: vi Makefile
```
- Now edit the 7th line to refer to the Fortran compiler of your choice:
```
f90comp = gfortran
```
- Compile the *fm3d* code.
```
sh$: make fm3d
```
- *Optional:* Copy the code to a convenient location, of your choice.
```
sh$: cp -r fm3d /usr/local/bin
```
- Add the path to the fm3d code to your $PATH environment variable. 
This is done through your .bash\_profile or .bashrc if the shell you are using is 
bash, .tcshrc if you are using csh etc...
```
sh$: vi ~/.bash_profile
```
```
export PATH=/usr/local/bin/fm3d:$PATH
```

- There are a set of input files needed by *fm3d*. Now is a good time to 
copy them into a data directory.
```
sh$: cp *.in foo/fm3d_input_files/
```
###Step 6: Install external focal mechanism inversion software *hashpy*
This step is not necessary for Earthquake location inversion, only for 
focal mechanism inversion.
- If you don't have *numpy*, install it first, the *hashpy* module uses 
numpy.distutils.core to compile Fortran code. The most recent version 
for your platform should work. If you are using *pip*, this is a one 
liner.
```
pip install numpy
```
Otherwise refer to http://docs.scipy.org/doc/numpy/user/install.html
- Clone the hashpy git repository
```
cd foo/git_repos
git clone https://github.com/markcwill/hashpy.git
```
- Navigate to the top-level directory of this repository and run setup script.
```
cd hashpy
python setup.py install
```
