# Steppyfinder
-------------------------------------
Created 2022
@author: Jacob Kerssemakers

# summary: 
Simplified python transcript of matlab code from [1] to detect steps in data traces

# description:
This program finds steps in noisy data traces by iterative fitting of single steps and evaluating the fit for each added step. While this Python development code is more basic in options, the architecture and used algorithms closely follow 'Autostepfinder' code formerly written in Matlab as associated with Loeff, Kerssemakers et al [1]. When using this code, please refer to [1]. The original algorithm was described in Kerssemakers [2]. 

# disclaimer:
Code is in early development and may still contain small bugs. Code may be subject to updates and changes. 

# settings
* input: single-column .txt file
* output: csv file columns:
	- index	
	- level before	
	- level after	
	- step size
	- dwell before	
	- dwell after	
	- predicted error	
	- measured error


# run options
* open AutoSteppyfinder, last line: multiPass(demo=0.1, tresH=0.15, N_iter=0) 
* demo number sets the different modii: 
	- first digit
	0 create a simple trace; 
	1 (default): load single trace
	2 run batch style on a directory

	- 2nd digit: 
	.0: only save; 
	.1 (default) show and save end graph 
	.2 show intermediate passes
	.3 show intermediate trace fit
	.4 show intermediate fit segments
	note: .2-.4 are suppressed in batch mode

* S_peak_treshold (0.15)
* iteration_range 
  (to run at max=tracelenth/4, set to 0; time consuming for large files)

# example: 
   A1_stepfind_multipass(demo=2.1, tresH=0.15, N_iter=100) 
    will run:
     -batch-style fit on a user-chosen directory with txt files, 
     ... with a fit round acceptance of 0.15
     ... with max. 100 iterations per fit
     ... and save both the step table and a .png fit picture to a 'result' directory 

## references:
[1] AutoStepfinder: A fast and automated step detection method for single-molecule analysis
Luuk Loef
f, Jacob W J Kerssemakers, Chirlmin Joo, Cees Dekker
Patterns . 2021 Apr 30;2(5):100256. 
doi: 10.1016/j.patter.2021.100256.
Matlab package on: https://zenodo.org/record/4657659#.Y6LxUXbMI2y

[2] Assembly dynamics of microtubules at molecular resolution
Jacob W J Kerssemakers et al
Nature 2006 Aug 10;442(7103):709-12. 
doi: 10.1038/nature04928.




