---
title: 'MyPTV: A Python Package for 3D Particle Tracking'
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
 - name: Swiss Federal Institute for Forest, Snow and Landscape Research WSL, 8903 Birmensdorf, Switzerland
   index: 1
date: 22 March 2022
bibliography: bib_myPTV.bib
---

# Summary

Three dimensional particle tracking velocimetry (3D-PTV) is a method that is widely used to study the dynamics of objects moving in space and to sample velocity fields at the location of tracer particles. Applications of 3D-PTV abound in various fields, such as, fluid mechanics, biology, animal behavior and crowd control [@Mass1993; @Virant1997; @Ott2000; @Luthi2005; @Ouellette2006; @Arneodo2008; @Holzner2008; @Toschi2009; @Shnapp2019; @Brizzolara2021; @Bagoeien2005; @Attanasi2015; @China2017; @Michalec2017; @Sinhuber2019; @Pouw2020]. A common methodology of 3D-PTV uses synchronized photography from several locations (e.g., by using several calibrated cameras, see \autoref{fig:setup}). From such camera images, the 3D positions of the particles can be estimated using photogrammetry methods. Then, particle locations are linked in time to generate 3D trajectories that can be analyzed. Furthermore, time differentiation yields measurements of the objects' velocity and acceleration, thus yielding the 3D particle tracking velocimetry method (3D-PTV) [@Virant1997]. In this work, we present an open source, Python-based software package, dedicated to making 3D-PTV more accessible to the scientific community.

![Left - A schematic sketch of a 3D-PTV experiment with a four-camera system. Right - the 6 steps of the post-processing and analysis of common 3D-PTV experiments. \label{fig:setup}](fig1.png)

# Statement of need

The application of 3D-PTV relies heavily on computation during the several steps of post-processing the experimental results (see \autoref{fig:setup}). In particular, in many applications researchers study objects that appear in high density in the recorded images, and the images are recorded over extended periods of time at high rates, yielding high volumes of raw data [@Shnapp2019]. Thus, 3D-PTV experiments are inevitably post-processed using specialized computer codes that often need to be tailored to the specific characteristics of the system under investigation. Indeed, the effort required to program a functioning 3D-PTV software might deter inexperienced researchers from employing 3D-PTV, and hinder further development of the method. 

_MyPTV_ is a Python package designed to make 3D-PTV accessible throughout the scientific community. _MyPTV_ is built on the foundation of a previous project - the first open source 3D-PTV software, _OpenPTV_ [@openptv], that was initiated in the early 2000's (although the algorithms and the code goes more than a decade earlier). The developers of _OpenPTV_ relied on coding in the C language in order to leverage its high speed of computation for implementing the complex algorithms involved. Another recent open source project, the C++ written _OpenLPT_ [@Tan2020], enables users to employ the shake-the-box algorithm [@Schanz2016] in 3D particle tracking experiments. While these two important projects enable using high-performance 3D-PTV, the C and C++ languages in which they are written are not accessible to many of the scientists working in the field. Therefore, debugging and installation can often be challenging as computers and operating systems evolve with time. Furthermore, new algorithms such as the ones introduced in [@Ouellette2006; @Xu2008; @Wieneke2013; @Schroeder2015; @Schanz2016; @Bourgoin2020; @Tan2020] have not yet been implemented in _OpenPTV_, and the development of novel algorithms, e.g. [@Brizzolara2021], is more restrictive in these low level coding languages. And yet, modern computers enable using 3D-PTV with higher level programming tools while maintaining computational times at a reasonable level. In addition to that, several tools exist for particle tracking in two dimensions [@Sbalzarini2005; @Heyman2019; @Allan2021], however they do not include 3D models and thus are limited to describing the motion in two dimensions. 

_MyPTV_ solves the above issues and extends _OpenPTV_ through three principles. First, _MyPTV_ is written exclusively in Python, which is a high level coding language which is accessible to a wider range of practitioners and widely used in scientific research. This feature allows rapid prototyping and development of the 3D-PTV method, which we believe is crucial for its further development. Second, the dependency on external packages is kept to the bare minimum and includes only a limited set of essential, widely used, and properly maintained packages (currently _Numpy_, _Scipy_, *Matplotlib*, _Pandas_, and *Pyyaml*), thus facilitating the maintenance and cross-platform usability without the need for complex deployment phases. Third, _MyPTV_ extends _OpenPTV_ by including new algorithms for camera calibration, particle tracking, particle segmentation, and trajectory smoothing, that were never implemented in _OpenPTV_. In particular, *MyPTV* includes a novel algorithm for the crucial stereo-matching step that uses time information to prioritize the correspondence of 3D-trackable trajectories, thus reducing the probabilities of stereo-matching ghost particles. Indeed, now that the code is more accessible, we envision that in the future _MyPTV_ will be further extended by its users to include more developments as they come.   

# Current capabilities

_MyPTV_, currently in version 0.4.3, contains all the necessary code needed to obtain three dimensional particle trajectories from a set of raw image data. In particular, this includes camera calibration, particle segmentation, stereo-matching, particle tracking, smoothing and stitching of broken trajectories. Each of these steps is built as a separate module of _MyPTV_, and generally contains a Python class or two used to perform a particular task. The code is written in an object-oriented style which is suited for the step-by-step structure of 3D-PTV.

# Testing MyPTV in an experiment

*MyPTV* was tested in a series of laboratory experiments. For example, in one of the experiments, seeding particles were tracked inside of a water filled tank. The flow was a moderate Reynolds number quasi-homogeneous turbulent flow generated through an 8-rotating wheels device [@Hoyer2005]. Images were taken at 50 frames per second per camera,  for a duration of 11.68 seconds using a three camera system. The camera resolution was $1280\times1024\,\,\text{pixels}^2$. The calibration, obtained through _MyPTV_'s calibration module, had a static calibration error of 84 microns, estimated through stereo-matching the 437 points of the calibration target. The particles in this test experiment, $50 \,\, \mu \text{m}$ in diameter, were tracked over a volume of $70\times70\times40$ mm$^3$. In each time step, about 850 particles were successfully linked in space and time. A 3D rendered image of particle trajectories obtained in the experiment is shown in \autoref{fig2}, showing a subset of 718 particle trajectories recorded during 3 seconds of the measurement.


![A 3D-rendered image, showing particle trajectories obtained in an experiment. The data shown corresponds to three seconds of measurement and shows 718 trajectories. \label{fig2}](traj_image.jpg)

# Documentation and usability 

_MyPTV_ includes several helpful tools to ensure the software's user friendliness. In particular, _MyPTV_ comes with a detailed user manual which outlines the instructions on how to use the software to achieve the desired results, and all of the functionalities of the various modules, including figures that demonstrate the various file formats used for saving the results of each module. In addition, the software includes an example data set that demonstrates the use of _MyPTV_ on real data. 

Furthermore, to enable users who are not experienced with Python to use the software, _MyPTV_ includes a dedicated "workflow" script used to run the various processing steps through a command line interface. Specifically, parameters for each particular experiment can be inserted by the users into a dedicated YAML file, and the workflow script can then be used to automatically perform any particular task. The results of the computations are then saved as text files following a tab-separated value format, which guarantees that the data can be analyzed with other softwares chosen by the users.

# Acknowledgements

The author is grateful for help in structuring the package and for numerous suggestions by Alex Liberzon, for fruitful discussions and suggestions by Markus Holzner, and Gal Schnapp, for help in conducting the test experiments and fruitful discussions with Stefano Brizzolara, and for the help in writing this paper from Dana Omer-Shnapp. Furthermore, the author wishes to acknowledge the significant contribution of the developers and the community of the _openPTV_ project to the development of the current software and the 3D-PTV method in general. 


# References
