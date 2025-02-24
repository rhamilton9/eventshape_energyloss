For running glauber on GRACE (grace.hpc.yale.edu, the Yale high performance computing cluster):

OPTION 1: (Simple)
	Use a pre-determined centrality reference for sorting events rather than sorting manually for each generation. This is faster, and can be done by submitting just a single batch. However, there is a systematic error introduced since the centrality bins can in principle vary from generation to generation. Using a pre-determined centrality class may cause there to be an asymmetric number of events in each centrality bin, which can be bad for some applications.

To do this option, first make sure that 


Option 2: (Accurate)
	Generate a centrality definition for each dataset. This is slower and requires 2 batch jobs, as well as 2 local bash script runs. But it provides a more accurate centrality bin definition, which can be especially important for calculations with the geometric mean.

Instructions:

1) Run the glauber generator split over many cores, using
	
	sbatch option2_RunGlauber.sh

Note that the settings for the generation(nucleus species, collision energy, etc.) can be set via the header file glauber.h while the number of events to generate is set within option2_RunGlauber.h. To modify the number of events to generate, one should open option2_RunGlauber.h using a text editor and change the variable nevent.

2) Merge the outputs of the 100 cores which distributed the event generation into a single file, and construct a centrality assignment for the events within the merged TTree. To do this, run

	bash  option2_MakeCentrality.sh

The centrality organization of the events is according to ascending impact parameter, which is standard for these kinds of monte carlo generators. The organization is implemented as an index assignment to the TTree entries via ROOT::TTreeIndex, and remains local to the tree within the file. Other macros can then easily access events in the tree according to their centrality. This allows us to sort events into centrality bins when constructing average events and merging histograms in the final step.

3) Create the averaged "net-event" for each centrality class using the centrality reference constructed in the last step, and add a TTree of information for these averaged net-events to the file. Just run

	sbatch option2_CalcArea.c

note that this will again parallelize the calculation by dividing events among may CPUs, and must be recombined in the final step.

4) Recombine all information into a single file with mergetree.c. First you must

	module load ROOT/6.26.10-foss-2022b

to access root in the terminal, then

	root -l -q mergetree.c

to merge the information. The completed glauber file will be stored in the /compiled directory.


