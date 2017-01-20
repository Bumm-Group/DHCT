DHCT software (thermal Drift Hysteresis and Creep Transform) is described in the article published in Review of Scientific Instruments. 
M.P. Yothers, A.E. Browder, L.A. Bumm, RSI, 2017.

This is a software package to determine and apply corrections to STM images due to the distortion
caused by thermal drift, hysteresis, and creep, named DHCT. The software was implemented in Matlab 2016a.
We are providing this software so that others can build upon it. If you would like to contribute 
your work so others can benefit, we would love to hear about it or you can create a fork to work on 
independently.

Origanization of this readme.txt

A. Copyright information
B. Contents and heirarchy
C. Parallel and GPU computation
D. Example how to use DHCT with sample data
E. Limitations


%------------

A. Copyright information

     copyright (c) 2016 Mitchell P. Yothers & Lloyd A. Bumm

     DHCT is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.
 
     DHCT is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
 
     You should have received a copy of the GNU General Public License
     along with DHCT.  If not, see <http://www.gnu.org/licenses/>.
  
     contact: 
     L. A. Bumm
     440 W Brooks St
     Homer L. Dodge Department of Physics & Astronomy
     The University of Oklahoma
     Norman, OK 73019
     bumm@ou.edu

%------------

B. Contents and hierarchy

The DHCT software consists of 13 matlab scripts, this readme, a text file describing the output 
data structure, another text file containing the GNU general public license, and sample data for
testing DHCT (GraphiteData.mat and Masks.mat). The hierarchial relationships of the DHCT Matlab
script suite is described below.

     HIERARCHICAL RELATIONSHIPS WITHIN DHCT
     DHCT.m (this is the top level script)
         data_mask_select.m
            -find_nearest_neighbors.m
            -position_analysis.m
         feat_position.m
         fitresult_modify.m
            -drift_models.m
         find_molecules.m
            -fmgaussfit_reduced.m
         ImCoords.m
            -drift_models.m
         ImReshape.m
         molecule_local_distortion.m
            -drift_models.m
         time_assign.m

%------------

C. Parallel and GPU computation

Many scripts in DHCT use parfor loops to take advantage of multiple CPU cores in parallel.

    SCRIPTS CONTAINING PARFOR LOOPS
    position_analysis.m
    molecule_local_distortion.m
    find_molecules.m   

If your version of Matlab does not have the parallel processing toolbox, simply change these three
parfor loops to for loops. 

position_analysis.m should detect automatically whether or not you have a compatible GPU device
for quickly determining feature positions and use it if it's available. 

%------------

D. Example how to use DHCT with sample data

Place all the DHCT files and the sample data into your Matlab active directory. Load the sample data. 

load('GraphiteData.mat');

GraphiteData.mat containes the following valiables:

 XScale: distance between pixels (meters)
 Period: elapsed time between pixels (seconds)
 SlowScanDirection: direction of slow scan, e.g. -y direction = 'down'
 FastScanDirection: direction of fast scan, e.g. +x direction = 'right', -x direction = 'left'
 NNdist: nearest neighbor distance for graphite in meters. (features in the STM image)

The images, Img1 and Img2 are 2D arrays of STM topographic image data where each element is 
a topographic height. Our images have the pixel values in nanometers, but height scale shouldn't matter. 

Matlab displays images in the 4th quadrant by default. Scan directions given are with reference to these
4th quadrant images.
 Img1 is the trace image, slow scan direction 'down', fast scan direction 'right'
 Img2 is the retrace image, slow scan direction 'down', fast scan direction 'left'

Run DriftCorrect with the following command:

 Data_Out = DHCT(XScale, Period, SlowScanDirection, FastScanDirection, NNdist, Img1, Img2);

The output data structure is described in detail in data_structure_description.txt

%------------

E. Limitations

Plane subtraction flattening of images is strongly recommended, and a way to reduce tip-fluctuation
induced noise may be helpful in order to allow DHCT to accurately determine features. Fogarty, et. al.
implements 'psfit' and 'xfilt' to do these tasks in Matlab, which you can read about in their paper at
<http://dx.doi.org/10.1063/1.2390633>. We used similar code of our own design on the sample images.

D. P. Fogarty, A. L. Deering, S. Guo, Z. Wei, N. A. Kautz and S. A. Kandel, Rev. Sci. Instrum. 77 (12),
    126104 - 126104-3 (2006).

DHCT assume the images are square arrays of pixels (Same number of pixels in x and y) and
that the scale in x is nominaly the same as in y. (The distance between pixels is roughly
the same in x and in y.) There is no inherent reason DHCT is restricted to square images 
with equal scaling in x and y, but it has only been tested with these images. DHCT also 
currently assumes that the analyzed image has a trigonal lattice structure - each feature
has six equidistant nearest neighbors. Some ideas about how to modify DHCT to handle other 
lattice structures can be found in the companion paper to this software published in RSI.

A handful of smaller limitations of the technique in general (rather than just our implementation 
of it) are also discussed in the paper.
<>



