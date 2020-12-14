# Obtain3D_open
Modified stereo-microscopy, example input files are shown in project branches. 

eftink@lanl.gov

Cite work derived from the results of this code:
B.P. Eftink, S.A. Maloy, Obtain3D_open, LANL, 2020

© 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
others to do so.

This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and
the following disclaimer.
 
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimer in the documentation and/or other materials provided with the
distribution.
 
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Description:
This computer code calculates the third dimension of information in tranmission electron micrographs using two images taken at different perspective of the same region. This is performed by a discrete-geometric approach that treats features as points, or series of points. Points correspond, but are not limited, to the center of cavities or precipitates, positions of irradiation black dot damage, positions along a dislocation line, or positions along where aninterface meets a free surface. Features can include dislocations, interfaces, cavities, precipitates, inclusions etc. The x, yand z coordinates of the points are output in a text file. The code also allows the user to visualize the features containing the points in three dimensions.

References using code:
1. B.P. Eftink, G.T. Gray III, S.A. Maloy, Stereographic Methods for 3D Characterization of Dislocations, Microscopy & Microanalysis, Volume 23, (2017) 
2. X. Hu et al., Transmutation-induced precipitation in tungsten irradiated with a mixed energy neutron spectrum, Acta Materialia, Volume 165, (2019)
3. K.G. Field, B.P. Eftink, C.M. Parish, S.A. Maloy, High-Efficiency Three-Dimensional Visualization of Complex Microstructures via Multidimensional STEM Acquisition and Reconstruction, Microscopy and Microanalysis, Volume 26, (2020)
