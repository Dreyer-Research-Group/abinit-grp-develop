---
authors: AM,ACGC
---

# The multibinit software

The MULTIBINIT software is using a second-principles approach for lattice dynamics simulations based on atomic potentials fitted on first-principles calculations [[cite:Wojdel2013]]. The second-principles model includes harmonic (short-range and long-range dipole-dipole interactions) and anharmonic parts as well as explicit treatment of homogeneous strains and their couplings to the lattice. Most parameters (harmonic IFC, elastic constants, strain-phonon coupling, ...) are provided directly from density functional perturbation theory (DFPT) calculations as implemented in ABINIT and keep therefore the first-principles accuracy. The anharmonic lattice part is treated in a more effective way and fitted [[cite:Escorihuela-Sayalero2017]] to reproduce stresses and forces of representative snapshots of short ab-initio molecular dynamics runs with a limited number of terms.


## 0 Installation  

Multibinit software is included in the [[help:abinit|ABINIT]] package, thus to install this code, you can follow the steps for the installation of the package, see installation section.

However, if you want to use all the features of Multibinit, you need to add some flags in the config.ac file:
  
 * Multibinit generates the output of the molecular dynamics in the NetCDF format, so you need to install and link the netCDF library with the flag:

        with_trio_flavor="netcdf"

*  Multibinit is using Libxml2 to parse the XML files (not mandatory but more efficient for very heavy XML files), you just have to add in the configuration file:
  
        enable_xml=yes
        CFLAGS_EXTRA="-I/usr/include/libxml2"


## 1 How to run the code?

### 1.1 Introducing the 'files' file

Given an input file (parameters described below) and the required files for the generation of the model. The user must create a "files" file which lists names for the files the job will require, including the main input file, the main output file, the file for the model (model.ddb or model.XML), the XML file for the anharmonic part of the model and the netcdf file for the training set (one per line).
The files file (called for example ab.files) could look like:
 
      multibinit.in
      multibinit.out
      model.ddb or model.XML
      model_anharmonic.MXL
      training_set_HIST.nc

In this example:

  * The main input file is called "multibinit.in".
  * The main output will be put into the file called "multibinit.out".
  * model.ddb or model.XML is the Database from abinit or XML.
  * model_anharmonic.MXL (optional) is the XML with the coefficients from fitted polynomial
  * training_set_HIST.nc is the history file in netcdf format with the training set

In short, a model.ddb or model.XML files contains the system definition and a list of derivatives of the total energy with respect to three kind of perturbations: phonons, electric field and strain. The optional XML file contains the list of coefficients from fitted polynomials. The last file is mandatory if you want to fit a polynome.

### 1.2 Running the code

The main executable file is called multibinit. Supposing that the "files" file is
called multibinit.files, and that the executable is placed in your working
directory, multibinit is run interactively (in Unix) with the command:

    multibinit < multibinit.files > log
  
or, in the background, with the command

    multibinit < multibinit.files > log &

where standard out and standard error are piped to the log file called "log"

The user can specify any names he/she wishes for any of these files. Variations of the
above commands could be needed, depending on the flavor of UNIX that is used
on the platform that is considered for running the code.

The syntax of the input file is strictly similar to the syntax of the main
abinit input files: the file is parsed, keywords are identified, comments are
also identified. However, the multidataset mode is not available.

## 2 Input variables
  
This Software is able to perform many different tasks. It is possible to generate a second principle model by extracting the harmonic part from DFPT calculation and fitting the anharmonic part into a polynomial expression. Such model can contains phonons instabilities and can be bound with an automatic procedure. Finally, for a given model, it is also possible to perform a molecular dynamics into multibinit. Each feature is controlled by a selected set of input variables. The 'flag' variables activates the different tasks e.g. [[multibinit:prt_model]], [[multibinit:fit_coeff]], [[multibinit:bound_model]] and [[multibinit:dynamics]]. The list of input variables for the multibinit input file are presented in the [[varset:multibinit]] variable set.

## 3 How to generate a model

Before to learn how to generate the model, we encourage you to read the papers [[cite:Wojdel2013]] and [[cite:Escorihuela-Sayalero2017]].

### 3.1 How to compute the harmonic part from [[help:abinit|ABINIT]]

Multibinit requires a DDB file from abinit to build the harmonic part of the potential. To compute this file, you first require to learn how to perform Density Functional Perturbation Theory ([[topic:DFPT]]) calculation with ABINIT. It is easier to use the different lesson of the tutorial: start with the second lesson on
response functions ([[lesson:rf2]]) to compute the DDB for the phonon part, then follow the lesson on elasticity ([[lesson:elastic]]) to compute the DDB for the elastic part. 
At the end, use [[help:mrgddb|DDB merge tool]] to generate the model.ddb.

To learn the procedure to compute the  harmonic part of the potential, you can follow this [[lesson:lattice_model | tutorial]]

### 3.2 How to compute the anharmonic part from a training set

To compute the anharmonic part of the model, multibinit requires a training set with several DFT calculations. 
Such training set will contain an ensemble of atomic configurations into the potential energy surface. 
This file will be used to fit a polynomial following the method developed in [[cite:Escorihuela-Sayalero2017]]. 
To generate the training set, the user requires to get some experience in molecular dynamics simulation within abinit. 

### 3.3 How to bound a model

Multibinit includes an automatic bounding process.
The tutorial is to be set up.
<!---- In this [[lesson:fit_process | tutorial]], you can get the basic knowledge to bound a model in multibinit. ---->

### 3.4 How to bound run a dynamics

Multibinit takes advantage of the mover routine of ABINIT which allows the users to perform atomic relaxation, molecular dynamics, and monte Carlo simulations.
The tutorial is to be set up.
<!----- In this [[lesson:lattice_model | tutorial]], you can get the basic knowledge to perform a dynamics into multibinit. ---->
