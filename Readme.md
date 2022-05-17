May 16, 2022
Version: 0.4.1


# MyPTV

MyPTV is an open source library designed for 3D particle Tracking Velocimetry (3D-PTV) measurements. In short, 3D-PTV is a method used to track the positions of particles in three dimensions; and is used extensively in fluid mechanics and can be applied to study other fields as well. The method relies on stereoscopic digital photography from calibrated cameras to infer particle 3D positions and track their motion through time. 

MyPTV builds heavily on the mathematical framework developed in the OpenPTV project (https://www.openptv.net/), and especially, the photogrammetry principles are the same [1], however it uses several new adaptations of the 3D-PTV algorithms that have been introduced in more recent years, e.g. [2-6], and several novel algorithms developed here.

## Why MyPTV?

So, why is MyPTV needed? While OpenPTV has been used and tested for many years with proven succses in generating scientific knowledge, the shortcoming of it is that it is written in C++, and therefore, it might be unaccessible to many physics and engineering scientists. This makes it hard for OpenPTV to keep up with the more modern approaches and algorithms, make it hard to debug in case of errors, and to maintain its support in newer operating systems. 

MyPTV solves these issues in two ways. First, we use the Python language, which is more accessible to the scientific, non-professional computer programmers, to handle. This accessibility to the inner workings of the code is essential for 3D-PTV to go into the future. Second, we attempt to keep the dependencies as few as possible and write as much of the code ourselves; this will help in maintaining the project usable for many years to come as Python keeps developing (currently, the only dependencies are a few functions from Numpy and Scipy, and optional plotting relied on Matplotlib). 

## Who is MyPTV for?

MyPTV is designed to be used by scientists and engineers who need to track the motion of objects in 3D. Applications range from fluid mechanics, biology, or medicine.  

## What does MyPTV include?

1) An imaging module that handles the 3D model
2) A calibration module that handles the camera calibration
3) A segmentation module that handles identifying and segmenting objects from images
4) A particle matching module that handles reconstructing 3D lab space particle coordinate estimation through triangulation
5) A tracking module for tracking particles in 3D
6) A trajectory smoothing module that can be used to smooth the results and calculate particle velocities and accelerations.

## How to install:

##### Requirements:

MyPTV requires you have Python 3 installed with pip, along with the Python packages: numpy, scipy, scikit-image, pandas, matplotlib, itertools

##### Installation:
###### Using `pip`

1) Open your terminal and change directory to the path of the code:
	`cd path/to/myptv` 
	
2) Finally, we use pip to install by using the following command: 
	`pip install .`
or 
	`pip install -r .\requirements.txt`

3) Optionally, parts of the code can be tested using pytest:
	`pytest ./tests/ -W ignore::RuntimeWarning`

4) Once this is done we are ready to go! You can now import MyPTV in your python code as usual. For example:
	`import myptv.imaging_mod`
or 	
   `from myptv import imaging_mod`

###### Using `conda` 

1) Install Anaconda or Miniconda and from the command shell inside the directory
where the package is downloaded:

	`conda env create -f environment.yml`
2) Activate the environment:

	`conda activate myptv`

3) Optionally, parts of the code can be tested using pytest:
	`pytest ./tests/ -W ignore::RuntimeWarning`

4) Once this is done we are ready to go! You can now import MyPTV in your python code as usual. For example:
	`import myptv.imaging_mod`
or 	
   `from myptv import imaging_mod`

## How to start?

Detailed instructions are given in the Manual, see `/user_manual/user_manual.pdf`.

## Who manages this project?

MyPTV was founded and is maintained by Ron Shnapp (ronshnapp@gmail.com). Contributions are most welcome. 

## References

[1] Maas, H. G., Gruen, A., & Papantoniou, D. (1993). Particle tracking velocimetry in three-dimensional flows. *Experiments in fluids*, *15*(2), 133-146.

[2] Mann, J., Ott, S., & Andersen, J. S. (1999). *Experimental study of relative, turbulent diffusion*. Risø National Laboratory.

[3] Ouellette, N. T., Xu, H., & Bodenschatz, E. (2006). A quantitative study of three-dimensional Lagrangian particle tracking algorithms. *Experiments in Fluids*, *40*(2), 301-313.

[4] Schanz, D., Schröder, A., Gesemann, S., Michaelis, D., & Wieneke, B. (2013). Shake the box: a highly efficient and accurate tomographic particle tracking velocimetry (TOMO-PTV) method using prediction of particle positions.

[5] Shnapp, R., Shapira, E., Peri, D., Bohbot-Raviv, Y., Fattal, E., & Liberzon, A. (2019). Extended 3D-PTV for direct measurements of Lagrangian statistics of canopy turbulence in a wind tunnel. *Scientific reports*, *9*(1), 1-13.

[6] Bourgoin, M., & Huisman, S. G. (2020). Using ray-traversal for 3D particle matching in the context of particle tracking velocimetry in fluid mechanics. *Review of scientific instruments*, *91*(8), 085105.
