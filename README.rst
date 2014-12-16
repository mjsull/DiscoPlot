DiscoPlot
=========

DiscoPlot allows the user to quickly identify genomic rearrangements, misassemblies and sequencing artefacts by providing a scalable method for visualising large sections of the genome. It reads single-end or paired read alignments in SAM, BAM or standard BLAST tab format and creates a scatter plot of opaque crosses representing the alignments to a reference. DiscoPlot is freely available (under a GPL license) for download (Mac OS X, Unix and Windows) at https://mjsull.github.io/DiscoPlot.

.. image:: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_3_lowres.gif
        :target: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_3.gif
        :alt: DiscoPlot figure

**DiscoPlot of a mock genome.** A mock genome was created by adding genomic rearrangements
to the chromosome of E. coli str. UTI89.  Paired-end reads generated from the mock genome (query)
with GemSim and mapped back to UTI89 (reference). The first ~500 Kbp were then visualised using DiscoPlot.


.. image:: https://pypip.in/version/DiscoPlot/badge.svg
        :target: https://pypi.python.org/pypi/DiscoPlot/
        :alt: Latest Version

.. image:: https://pypip.in/download/DiscoPlot/badge.svg
        :target: https://pypi.python.org/pypi/DiscoPlot/
        :alt: Downloads

.. image:: https://travis-ci.org/mscook/DiscoPlot.svg?branch=master
        :target: https://travis-ci.org/mjsull/DiscoPlot
        :alt: Build status

.. image:: https://landscape.io/github/mjsull/DiscoPlot/master/landscape.png
        :target: https://landscape.io/github/mjsull/DiscoPlot/master
        :alt: Code Health

Documentation
-------------

Please use this README.rst as the core DiscoPlot user documentation. 


News
----

**4/12/14:** Version 0.5.0 released.


Citation
--------

Cite this Github repository if you use DiscoPlot to generate figures 
for publications:: 

    SULLIVAN MJ, BEATSON SA^. 
    DiscoPlot - Visualising discordant readss.
    https://github.com/mjsull/DiscoPlot.

TODO:
 * Sam compatibility
 * Print selected read names, alignments or sequences


Installation
------------

DiscoPlot is a commandline application. If you're not familiar with the 
commandline we recommend you ask local IT support to help you install it.

N.B. These instructions are a placeholder, to run DiscoPlot please clone this git
and manually install pysam/matplotlib/numpy/NCBI-BLAST+

Pip installation coming soon.

We test DiscoPlot builds on both Linux (Ubuntu >= 12.04) and Windows 8 system. 

You will need to install/have installed:
    * python >= 2.7 (**Python 3 is not supported**)

To automatically generate BLAST aligments (For long read DiscoPlots) using DiscoPlot you will need to install/have installed:
    * ncbiblast >= 2.2.27
    
You can check these are installed by::
    
    $ python --version
    $ blastn -version

Installation of python or blastn (without a package manager) is beyond the 
scope of this document.

If you have both python and blastn you need to (if not already present) 
install pip_.

You can check if pip_ exists with::

    $ which pip

If you get a "not found", please read the `pip installation instructions`_. 

**If you already have pip we do suggest you upgrade it.** We are using version 
1.5.6 at the time of writing this document. 

You can upgrade pip_ like this::

    $ pip install --upgrade pip


The following python libraries_ should be installed (automatically) if you follow 
the installation instructions detailed below.

We use the following python libraries_:
    * numpy >= 1.6.1
    * matplotlib >= 1.1.0
    * pysam >= 0.8.1

Pysam is only required for generating DiscoPlots with BAM files. SAM compatability
has been included to allow windows users to generate DiscoPlots. PySam will not
install on Windows, don't bother trying (or if you've succeeded please let me know how).


Linux (Ubuntu)
~~~~~~~~~~~~~~

Discoplot uses 3rd party packages that are extremely important for scientific 
computing but may be difficult to install. While *pip install * 
*--user DiscoPlot* may work we recommend you install these 3rd party packages 
using apt-get.

Run::

    $ sudo apt-get install python-numpy python-matplotlib 

Now pip_ install DiscoPlot::
    
    $ pip install --user DiscoPlot

We use the --user option of pip_ to put DiscoPlot in: /home/$USER/.local/bin/
You need to add this location to you ~/.bash_profile. 

Add SeqFindr to your path::

    $ echo 'export PATH=$PATH:/home/$USER/.local/bin/' >> ~/.bash_profile

Finally install BLAST+::

    $ sudo apt-get install ncbi-blast+ 

**Test it:**

Run::
    
    $ DiscoPlot -h 
    $ python -c 'import DiscoPlot; print DiscoPlot'


MacOSX (Mavericks)
~~~~~~~~~~~~~~~~~~

**You'll need to have the equivalents of python-dev libatlas-dev liblapack-dev 
gfortran libfreetype6-dev libfreetype6 & libpng-dev installed.** We had no 
problems installing SeqFindr on a recently acquired OSX Mavericks machine 
using the homebrew package manager.

The installed packages on this machine via::

    $ brew list 

Are available at this gist_.

pip install DiscoPlot::
    
    $ pip install --user DiscoPlot

We use the --user option of pip_ to put SeqFindr in: /home/$USER/.local/bin/
You need to add this location to you ~/.bash_profile. 

Add DiscoPlot to your path::

    $ echo 'export PATH=$PATH:/home/$USER/.local/bin/' >> ~/.bash_profile

Finally install BLAST+::

    $ sudo brew install blast 

**Test it:**

Run::
    
    $ DiscoPlot -h 
    $ python -c 'import DiscoPlot; print DiscoPlotr'


Upgrading DiscoPlot
~~~~~~~~~~~~~~~~~~

You can upgrade like this::
    
    pip install --upgrade DiscoPlot


**Please regularly check back to make sure you're running the most recent 
DiscoPlot version.**



Example of figures produced by DiscoPlot
----------------------------------------

.. image:: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_3_lowres.gif
        :target: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_3.gif
        :alt: DiscoPlot figure
        :align: center

**DiscoPlot of a mock genome.** A mock genome was created by adding genomic rearrangements to the chromosome of E. coli 
str. UTI89.  Paired-end reads generated from the mock genome (query) with GemSim (ref) and mapped back to UTI89 (reference).
The first ~500 Kbp were then visualised using DiscoPlot.

.. image:: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_4_lowres.gif
    :target: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_4.gif
    :alt: DiscoPlots of structural variants
    :align: center

**DiscoPlots of common structural variants.** Each box shows a common genomic rearrangement represented by a DiscoPlot.
Rows A and B were created using 100 bp long paired-end reads with an insert size of 300bp.
Rows C and D were created using single-end reads with an average length of 1000bp. 
For each box the rearrangement in the sequenced genome is listed, followed by the scale of the gridlines in brackets.
A1,  C1: 300 bp deletion (400 bp). A2, C2: 300 bp insertion (400 bp). A3, C3: 300 bp inversion (400 bp).
A4, C4: 300 bp sequence translocated 50 Kbp upstream (10 Kbp). B1, D1: 3000 bp deletion (1000 bp). 
B2, D2: 3000 bp insertion (500 bp). B3, D3: 3000 bp inversion (1000 bp). 
B4, D4: 3000 bp sequence translocated 50 Kbp upstream (10 Kbp). C1) 

.. image:: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_5_lowres.png
    :target: https://raw.githubusercontent.com/mjsull/DiscoPlot/master/pictures/Figure_5.png
    :alt: DiscoPlot of E. coli genome
    :align: center

**The dynamic nature of the genome of Escherichia coli str. UTI89.** Discoplot of paired-end reads from a clonal 
culture of UTI89 mapped back to the published reference chromosome and plasmid (top). A) Zoomed region of the
DiscoPlot, a small inversion exists in some of the sequenced bacteria. Four of these sites, corresponding to
known prophage regions, were identified using DiscoPlot. B) Close up of the plasmid in the DiscoPlot. Each entry
in the alignment file is separated by an opaque green line. A large inversion has been identified, this region
corresponds to an inverted repeat found in the plasmid. The cross in the lower right corner indicates that this
region circularises.

Tutorials
---------

**Coming Soon**

Commands
--------

To see a full list of flags type DiscoPlot --help

Detailed descriptions coming soon
