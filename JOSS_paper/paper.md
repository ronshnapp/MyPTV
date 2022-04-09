---
title: 'MyPTV: A Python package for 3D particle tracking'
tags:
  - Python
  - 3D particle tracking
  - Fluid mechanics
  - Turbulence
  - Experimental methods
authors:
  - name: Ron Shnapp
    orcid: 0000-0001-7495-8420
    affiliation: 1
affiliations:
 - name: Swiss Federal Institute for Forest, Snow and Landscape Research WSL
   index: 1
date: 22 March 2022
bibliography: bib_myPTV.bib
---

# Summary

Three dimensional particle tracking velocimetry (3D-PTV) is a method that is widely used to study the dynamics of objects moving in space and to sample velocity fields at the location of tracer particles. Applications of 3D-PTV abound in various fields, such as, fluid mechanics, biology, animal behavior and crowd control [@Mass1993; @Virant1997; @Ott2000; @Luthi2005; @Ouellette2006; @Arneodo2008; @Holzner2008; @Toschi2009; @Shnapp2019; @Brizzolara2021; @Bagoeien2005; @Attanasi2015; @China2017; @Michalec2017; @Sinhuber2019; @Pouw2020]. A common methodology of 3D-PTV uses synchronized photography from several locations (e.g., by using several calibrated cameras, see \autoref{fig:setup}). From such camera images, the 3D positions of the particles can be estimated using photogrammetry methods. Then, particle locations are linked in time to generate 3D trajectories that can be analyzed. Furthermore, time differentiation yields measurements of the objects' velocity and acceleration, thus yielding the 3D particle tracking velocimetry method (3D-PTV) [@Virant1997]. In this work, we present an open source, Python-based software package, dedicated to making 3D-PTV more accessible to the scientific community.

![Left - A schematic sketch of a 3D-PTV experiment with a four-camera system. Right - the 6 steps of the post-processing and analysis of common 3D-PTV experiments. \label{fig:setup}](fig1.png)

# Statement of need

The application of 3D-PTV relies heavily on computation during the several steps of post-processing the experimental results (see \autoref{fig:setup}). In particular, in many applications researchers study objects that are tightly packed in the recorded images, and the images are recorded over extended periods of time at high rates, yielding high volumes of raw data [@Shnapp2019]. Thus, 3D-TV experiments are inevitably post-processed using specialized computer codes that often need to be tailored to the specific characteristics of the system under investigation. Indeed, the effort needed to be put into the programming of a functioning 3D-PTV software might deter inexperienced researchers from employing 3D-PTV and it might hinder further development of the method. 

_MyPTV_ is a Python package designed to make 3D-PTV accessible throughout the scientific community, building on the foundations of a previous project. Indeed, the first open source 3D-PTV softwer is the _OpenPTV_ project, that was initiated at the early 2000's [@openptv]. The developers of _OpenPTV_ relied on coding in the C language in order to leverage it's high speed of computation for resoling the complex algorithms involved. Nevertheless, the C language is not accessible to many of the scientists working in the field, so debugging and installation on modern computers can often be challenging, and recent algorithms that have been proposed in the recent years to advance the method (e.g. @Ouellette2006; @Xu2008; @Schroeder2015; @Bourgoin2020; @Brizzolara2021) have not yet been implemented in _OpenPTV_. Furthermore, modern computers make applying 3D-PTV using higher level programming tools while maintaining computational times at a reasonable level. Thus, _MyPTV_ solves these issues and extends _OpenPTV_ through three principles. First, _MyPTV_ written exclusively in Python, which is accessible to a wider range of practitioners and widely used in scientific research. This feature allows rapid prototyping and development of the 3D-PTV method which is crucial for its further development. Second, the dependencies on external packages is kept to the bare minimum and includes only an essential set of widely used and properly maintained packages (currently _Numpy_, _Scipy_ and _Pandas_), thus facilitating the maintenance and cross-platform usability without the need for complex deployment phases.  Third, _MyPTV_ extends _OpenPTV_ by including new algorithms for camera calibration, particle tracking, particle segmentation, and trajectory smoothing that were never implemented in _OpenPTV_. In particular, a novel algorithm for the crucial stereo-matching step was developed that uses time information to prioritize the correspondence of 3D-trackable trajectories.  Indeed, with the code being more accessible, we envision that _MyPTV_  will be further extended in the future by its users to include more developments as they come.   

# Current capabilities

_MyPTV_, currently in version 0.1, contains all the necessary code needed to obtain three dimensional particle trajectories from a set of raw image data. In particular, this includes camera calibration, particle segmentation, stereo-matching, particle tracking, smoothing and stitching of broken trajectories. Each of these steps is built as a separate module of _MyPTV_ and generally contains a Python class or two used to perform the particular task needed. The code is written in an object oriented style which is suited for the step-oriented structure of 3D-PTV.

# Tests

*MyPTV* had been tested in a series of laboratory experiments. For example, in one of the experiments seeding particles were tracked in moderate Reynolds numbers turbulent flows generated through an 8-rotating wheels device [@Hoyer2005]. Images were taken at 50 frames per second per camera,  for a duration of 11.68 seconds using a three camera system. The camera resolution was $1280\times1024\,\,\text{pixels}^2$. The calibration, obtained through MyPTV's calibration module, had a static calibration error of 84 microns, estimated through stereo-matching the 437 points of the calibration target. The particles in our experiment, $50 \,\, \mu \text{m}$ in diameter, were tracked over a volume of $70\times70\times40$ mm$^3$. In each time step, about 850 particles were successfully linked in space and time. A 3D rendered image of particle trajectories obtained in the experiment is shown in \autoref{fig2}, showing a subset of 718 particle trajectories recorded during 3 seconds of the measurement.

![A 3D-rendered image, showing particle trajectories obtained in an experiment. The data shown corresponds to three seconds of measurement and shows 718 trajectories. \label{fig2}](traj_image.jpg)

# Documentation and usability 

The MyPTV repository contains a detailed user manual which outlines all of the functionalities of the various modules, including figures that demonstrate the various file formats used for saving the results of each module. In addition, there is a detailed guide on the process of camera calibration and an example using *Jupyter notebook*. Furthermore, _MyPTV_ can easily be installed through _pip_, the Python package installer, and automated tests are included for each of the various modules.

# Acknowledgements

The author acknowledges fruitful discussions with Alex Liberzon, Markus Holzner, and Gal Schnapp, and the help in the testing experiment and fruitful discussions with Stefano Brizzolara. 

# References