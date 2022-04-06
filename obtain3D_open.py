import os, sys
import math
import vector_math  # file with vector/matrix functions
from tkinter import filedialog # alternative for file/folder dialog
import pandas as pd # who doesn't like a good panda?
import numpy as np # math atan2 doesn't work on dataframes, apparently

version = '20220721'

print('\n')
print('Cite work derived from the results of this code: B.P. Eftink, S.A. Maloy, Obtain3D_open, LANL, 2020')
print('\n')
print('© 2020. Triad National Security, LLC. All rights reserved. This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.')
print('\n')
print('This program is open source under the BSD-3 License. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:')
print('1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.')
print('2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.')
print('3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.')
print('\n')
print('THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.')
print('\n')
print(f'Version: {version}')

input_dir = filedialog.askdirectory(title='Select input file directory')    # use tkinter to fetch input directory
os.chdir(input_dir) # strictly speaking this isn't necessary, but I'll duplicate the previous behavior.

#%%##############################################################################
# Open files in working directory
data = pd.read_csv('in.txt', delimiter = '\t', header=None, skip_blank_lines=True, names=['x1', 'y1', 'x2', 'y2', 'dis_number', 'point_number']).dropna()
key_data = pd.read_csv('key.txt', delimiter = '\t', header=None, index_col = 0, skip_blank_lines=True).dropna().squeeze("columns").to_dict()    # names=['x1', 'y1', 'x2', 'y2', 'dis_number', 'point_number']
thetatr = math.radians(key_data['x_tilt_2']-key_data['x_tilt_1'])   # calculates the difference in tilt between the two images (in radians)

# Read in particle info, if needed
if key_data['Spherical_particles_or_cavities'] != 0:
    part_data = pd.read_csv('part.txt', delimiter = '\t', header=None, skip_blank_lines=True, names=['setp', 'part_size', 'type_part']).dropna()

#%%#####################################################################################################################
# Find reference axis and rotate points appropriately
ref = 0
position = pd.DataFrame({'x1':data.x1-data.x1[ref],
    'y1':data.y1[ref]-data.y1,
    'x2':data.x2-data.x2[ref],
    'y2':data.y2[ref]-data.y2})
position['axis'] = pd.DataFrame({'axis': np.arctan2(position.y1-position.y2, position.x2-position.x1)})
position.axis[position.axis < 0] = position.axis[position.axis < 0] + math.pi   # add pi to negative values to flip with respect to origin

axisavg = position.drop(0).axis.mean()  # have to drop the reference point so it doesn't get averaged in
axisdev = position.drop(0).axis.std()
################## z is calculated for each point seperately
angle_0 = -1*((math.pi/2)-position.drop(0).axis.loc[(position.axis > (axisavg - axisdev)) & (position.axis < (axisavg + axisdev))].mean())
axis_0 = [0,0,1]

xy1 = pd.DataFrame(vector_math.rot(axis_0[0],axis_0[1],axis_0[2], angle_0, position.x1, position.y1, 0)).transpose()#
xy1.columns = ['x', 'y', 'z'] # rename column headings
xy2 = pd.DataFrame(vector_math.rot(axis_0[0],axis_0[1],axis_0[2], angle_0, position.x2, position.y2, 0)).transpose()
xy2.columns = ['x', 'y', 'z']
delth = (xy2.y - xy1.y)/(2.0*np.sin(thetatr/2.0))
z1 = (-1.0)*(np.tan(thetatr/2.0)*xy1.y + delth/(np.sin((math.pi/2)-(thetatr/2.0))))
anglex = (-1)*math.radians(key_data['x_tilt_1'])
xyz0 = pd.DataFrame(vector_math.rot(1,0,0,anglex,xy1.x,xy1.y,z1)).transpose()
xyz0.columns = ['x', 'y', 'z']
angley = math.radians(key_data['y_tilt_1'])
xyz00 = pd.DataFrame(vector_math.rot(0,1,0,angley,xyz0.x,xyz0.y,xyz0.z)).transpose()
xyz00.columns = ['x', 'y', 'z']
output = pd.DataFrame({'x': xyz00.x/key_data['scale_pixels/nm'],
    'y': xyz00.y/key_data['scale_pixels/nm'],
    'z': xyz00.z/key_data['scale_pixels/nm']})
#
try:    # if there are spheres, drop in the sphere size too, otherwise just forget it.
    output = pd.concat([output, part_data.part_size, part_data.type_part], axis=1)
    print('Adding particle data to output.')
except:
    print('No particle data... Skipping output.')

if data.dis_number.nunique() > 1:
    try:
        output = pd.concat([output, data.dis_number], axis=1)
        print('Adding dislocation data to output.')
    except:
        pass
else:
    print('No dislocation data... Skipping output.')

output = output.reset_index() # just generates a new column for the index that can be exported as the point number
output.rename(columns = {'index':'number'}, inplace = True) # rename index so we can
output.number += 1  # increment this so it's consistent with original code
output.to_csv('coordinates.txt', sep = '\t', index = False)  # save as text file
print('Data written to coordinates.txt. Processing complete!')


#%%##################### Plotting (plot 1) z is in beam direction ########################################################
if key_data['view'] > 0.5:  # only need this part if view is selected.
    pass
else:
    sys.exit()  # We can exit here if we're not plotting. I'm mostly doing this to remove 1 layer of nested indents

import mayavi.mlab as mlab
#import matplotlib
##matplotlib.use('QT5Agg')
#from mpl_toolkits import mplot3d   # not required, just keeping track of it.
import matplotlib.pyplot as plt # Maybe we'll try this instead of mayavi

# set up the constants needed for the other parts
outside = 10
x_box_min = min(output.x)-outside
x_box_max = max(output.x)+outside
y_box_min = min(output.y)-outside
y_box_max = max(output.y)+outside
z_box_min = min(output.z)-outside
z_box_max = max(output.z)+outside
tube = 0.5
bounding_box = pd.DataFrame()
bounding_box['linex_b']=[[x_box_min, x_box_max], [x_box_min, x_box_min], [x_box_min, x_box_max], [x_box_max, x_box_max], [x_box_min, x_box_min], [x_box_max, x_box_max], [x_box_max, x_box_max], [x_box_min, x_box_min], [x_box_min, x_box_max], [x_box_min, x_box_min], [x_box_min, x_box_max], [x_box_max, x_box_max]]
bounding_box['liney_b']=[[y_box_min, y_box_min], [y_box_min, y_box_max], [y_box_max, y_box_max], [y_box_min, y_box_max], [y_box_min, y_box_min], [y_box_min, y_box_min], [y_box_max, y_box_max], [y_box_max, y_box_max], [y_box_min, y_box_min], [y_box_min, y_box_max], [y_box_max, y_box_max], [y_box_min, y_box_max]]
bounding_box['linez_b']=[[z_box_min, z_box_min], [z_box_min, z_box_min], [z_box_min, z_box_min], [z_box_min, z_box_min], [z_box_min, z_box_max], [z_box_min, z_box_max], [z_box_min, z_box_max], [z_box_min, z_box_max], [z_box_max, z_box_max], [z_box_max, z_box_max], [z_box_max, z_box_max], [z_box_max, z_box_max]]
#
dis_color = (0,0,0) ### color for the dislocation segments
color_key = pd.DataFrame({'color_code':((0,1,0),(0,0,0),(1,0,1),(1,0,0),(1,1,0)),
    'color_name':('green', 'black', 'purple', 'red', 'yellow')})  # table of color options, green, black, purple, red, yellow
#
# Calculate which points belong together in a dislocation.
# NOTE: I didn't move this to pandas. There might be an easy way to flag this, but I think it depends heavily on viz program. For now, I'll just dump the dislocation number back out at the end. -morrow
# NOTE 2: This was moved to the plotting portion, because it's not actually required unless you want to plot.
if data.dis_number.nunique() > 1:   # we don't need to do this if there's no dislocation data
    segs = []
    segs_dis = []
    for count in range(0,5100,1):
        groups = []
        for i in range(0,len(data['dis_number']),1):
            if abs(count-data['dis_number'][i]) < 0.1:
                groups.append(data['point_number'][i])
        if sum(groups) > 0:
            groups.append(count)
            segs.append(groups)
            segs_dis.append(groups)
        #
    dis_segs = []
    for i in range(0,len(segs),1):      ################# dis_segs splits dislocations into 2 point segments for plotting
        for j in range(0,len(segs[i])-2,1):
            dis_segs.append([segs[i][j],segs[i][j+1],segs[i][-1]])
#
    linex = []
    liney = []
    linez = []
    for i in range(0,len(dis_segs),1):
        if dis_segs[i][-1] > 0 and dis_segs[i][-1] <5000:
            star = int(dis_segs[i][0])
            finish = int(dis_segs[i][1])
            linex.append([output.x[star-1], output.x[finish-1]])
            liney.append([output.y[star-1], output.y[finish-1]])
            linez.append([output.z[star-1], output.z[finish-1]])

if key_data['movie'] < 0.5: # View only, no movie
    mlab.figure(bgcolor=(1,1,1))
    # Create bounding box
    for i in range(len(bounding_box)):
        mlab.plot3d(bounding_box.linex_b[i], bounding_box.liney_b[i], bounding_box.linez_b[i], tube_radius = tube, color = (0,0,0))
    # Make line segments between points of the same dislocation
    try:
        for i in range(len(linex)):
            mlab.plot3d(linex[i], liney[i], linez[i], tube_radius = 2, color = dis_color)
    except:
        pass
    # plot number markers
    if key_data['Number_markers_on_plot'] > 0:
        for i in range(0,len(output.x),1):
            mlab.text3d(output.x[i],output.y[i],output.z[i], str(int(i)), color = (0.0,0.0,0.0), scale = 5.0)
    # plot particles, voids
    try:
        for i in range(0,len(part_data.setp),1):
            try:
                colour = color_key.color_code[part_data.type_part[i]-1]    # color_key was defined at the beginning of plotting
            except:
                colour = (0,0,0)   # default to black if there's some error.
            mlab.points3d(output.x[part_data.setp[i]-1], output.y[part_data.setp[i]-1], output.z[part_data.setp[i]-1], color= colour, mode = 'sphere', scale_factor = part_data.part_size[i]/key_data['scale_pixels/nm'], resolution = 64)
    except:
        pass
    # extra options
    mlab.gcf().scene.parallel_projection = True
    mlab.gcf().scene.show_axes = True

#%%#####################################
##    Matplotlib version
if key_data['movie'] < 0.5:
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    # plot particles/voids
    try:    # need this in case we don't have particles
        ax.scatter(output.x, output.y, output.z, zdir='z', s=output.part_size, depthshade = True, edgecolor='k', color = color_key.color_code[output.type_part-1])  #c=None
    except:
        pass
#    ax.set_title('Cavities')
    ax.set_xlabel('x [nm]')
    ax.set_ylabel('y [nm]')
    ax.set_zlabel('z [nm]')
    ax.set_xlim(x_box_min, x_box_max)  # if you want to specify limits manually. Probably easier to specify buffer instead
    ax.set_ylim(y_box_min, y_box_max)
    ax.set_zlim(z_box_min, z_box_max)
    #ax.set_xmargin(0.1) # margin between 0 and 1 (times maximum range)
    #ax.set_ymargin(0.1)
    #ax.set_zmargin(0.1)
    x_range = x_box_max - x_box_min
    y_range = y_box_max - y_box_min
    z_range = z_box_max - z_box_min
    rangemax = max(x_range, y_range, z_range)
    ax.set_box_aspect((x_range/rangemax,y_range/rangemax,z_range/rangemax))
    ax.set_proj_type('ortho')
    ax.grid(True)
    #ax.tick_params(axis="y",direction="in", pad=-22)
    ax.tick_params(direction="out")
    fig.set_facecolor('white') # 'black', 'white'
    ax.set_facecolor((0.0, 0.0, 0.0, 0.0))
    ax.w_xaxis.pane.fill = False    # turn off the face fill
    ax.w_yaxis.pane.fill = False
    ax.w_zaxis.pane.fill = False
    # Plot bounding box.
    for i in range(len(bounding_box)):
        ax.plot3D(bounding_box.linex_b[i], bounding_box.liney_b[i], bounding_box.linez_b[i], linewidth = tube, color = (0,0,0))
    # Plot dislocation line segments
    try:
        for i in range(len(linex)):
            ax.plot3D(linex[i], liney[i], linez[i], linewidth = 2, color = dis_color)
    except:
        pass
    # Plot number markers
    if key_data['Number_markers_on_plot'] > 0:
        for i in range(0,len(output.x),1):
            fig_text = str(int(i))
            ax.text3D(output.x[i], output.y[i],output.z[i], output.number[i], color = (0.0,0.0,0.0), size = 5.0)
    #
    plt.show()

#%%################# Movie Mode ##########################
if key_data['movie'] > 0.1:
    track = 0
    aziz_start = [0]
    aziz_end = [0]
    elevation_start = [1080]
    elevation_end = [720]
    angle_step = [-1]
    roll_start = [0]
    roll_end =  [0]

    for b in range(0,len(aziz_start),1):
        if abs(aziz_end[b]-aziz_start[b]) > 1:
            print("aziz")
            for x in range(aziz_start[b],aziz_end[b],angle_step[b]):
                track = track + 1
                mlab.figure(bgcolor=(1,1,1))
                aziz = x
                elevation = elevation_start[b]
                roll_value = roll_start[b]
                print(aziz)
                print(elevation)
                mlab.close()
                mlab.figure(bgcolor=(1,1,1))
                fig_size = (600,600)
                #
                # Create bounding box
                for i in range(len(bounding_box)):
                    mlab.plot3d(bounding_box.linex_b[i], bounding_box.liney_b[i], bounding_box.linez_b[i], tube_radius = tube, color = (0,0,0))
                #
                # Make line segments between points in the same dislocation
                linex = []
                liney = []
                linez = []
                for i in range(0,len(dis_segs),1):
                    if dis_segs[i][-1] > 0 and dis_segs[i][-1] <5000:
                        star = int(dis_segs[i][0])
                        finish = int(dis_segs[i][1])
                        linex = [output.x[star-1], output.x[finish-1]]
                        liney = [output.y[star-1], output.y[finish-1]]
                        linez = [output.z[star-1], output.z[finish-1]]
                        mlab.plot3d(linex, liney, linez, tube_radius = 2, color =dis_color)
                #
                if key_data['Number_markers_on_plot'] > 0:
                    for i in range(0,len(output.x),1):
                        mlab.text3d(output.x[i],output.y[i],output.z[i], str(int(i)), color = (0.0,0.0,0.0), scale = 5.0)
                #
                for i in range(0,len(part_data.setp),1):
                    try:
                        colour = color_key.color_code[part_data.type_part[i]-1]    # color_key was defined at the beginning of plotting
                    except:
                        colour = (0,0,0)   # default to black if there's some error.
                    mlab.points3d(output.x[part_data.setp[i]-1], output.y[part_data.setp[i]-1], output.z[part_data.setp[i]-1], color= colour, mode = 'sphere', scale_factor = part_data.part_size[i]/key_data['scale_pixels/nm'], resolution = 64)
                #
                mlab.gcf().scene.parallel_projection = True
                mlab.gcf().scene.show_axes = True
                mlab.show()
                file_name = str(track) + '.tiff'
                mlab.view(aziz,elevation, roll = roll_value)
                mlab.savefig(file_name,size=fig_size)
                mlab.close()
        #
        if abs(elevation_end[b]-elevation_start[b]) > 1:
            print("elevation")
            for x in range(elevation_start[b],elevation_end[b],angle_step[b]):
                track = track + 1
                mlab.figure(bgcolor=(1,1,1))
                aziz = aziz_start[b]
                elevation = x
                roll_value = roll_start[b]
            #
                mlab.close()
                mlab.figure(bgcolor=(1,1,1))
                fig_size = (600,600)
            #
            # Create bounding box
                for i in range(len(bounding_box)):
                    mlab.plot3d(bounding_box.linex_b[i], bounding_box.liney_b[i], bounding_box.linez_b[i], tube_radius = tube, color = (0,0,0))
            #
                linex = []
                liney = []
                linez = []
                for i in range(0,len(dis_segs),1):
                    if dis_segs[i][-1] > 0 and dis_segs[i][-1] <5000:
                        star = int(dis_segs[i][0])
                        finish = int(dis_segs[i][1])
                        linex = [output.x[star-1], output.x[finish-1]]
                        liney = [output.y[star-1], output.y[finish-1]]
                        linez = [output.z[star-1], output.z[finish-1]]
                        mlab.plot3d(linex, liney, linez, tube_radius = 2, color = dis_color)
            #
                if key_data['Number_markers_on_plot'] > 0:
                    for i in range(0,len(output.x),1):
                        mlab.text3d(output.x[i],output.y[i],output.z[i], str(int(i)), color = (0.0,0.0,0.0), scale = 5.0)
            #
                for i in range(0,len(part_data.setp),1):
                    try:
                        colour = color_key.color_code[part_data.type_part[i]-1]    # color_key was defined at the beginning of plotting
                    except:
                        colour = (0,0,0)   # default to black if there's some error.
                    mlab.points3d(output.x[part_data.setp[i]-1], output.y[part_data.setp[i]-1], output.z[part_data.setp[i]-1], color= colour, mode = 'sphere', scale_factor = part_data.part_size[i]/key_data['scale_pixels/nm'], resolution = 64)
            #
                mlab.gcf().scene.parallel_projection = True
                mlab.gcf().scene.show_axes = True
                mlab.show()
                file_name = str(track) + '.tiff'
                mlab.view(aziz,elevation, roll = roll_value)
                mlab.savefig(file_name,size=fig_size)
                mlab.close()

    for b in range(0,len(roll_start),1):
        if abs(roll_end[b]-roll_start[b]) > 1:
            print("roll")
            for x in range(roll_start[b],roll_end[b],angle_step[b]):
                track = track + 1
                mlab.figure(bgcolor=(1,1,1))
                aziz = aziz_start[b]
                elevation = elevation_start[b]
                roll_value = x
            #
                mlab.close()
                mlab.figure(bgcolor=(1,1,1))
                fig_size = (600,600)
            # # Create bounding box
                for i in range(len(bounding_box)):
                    mlab.plot3d(bounding_box.linex_b[i], bounding_box.liney_b[i], bounding_box.linez_b[i], tube_radius = tube, color = (0,0,0))
        #
                linex = []
                liney = []
                linez = []
                for i in range(0,len(dis_segs),1):
                    if dis_segs[i][-1] > 0 and dis_segs[i][-1] <5000:
                        star = int(dis_segs[i][0])
                        finish = int(dis_segs[i][1])
                        linex = [output.x[star-1], output.x[finish-1]]
                        liney = [output.y[star-1], output.y[finish-1]]
                        linez = [output.z[star-1], output.z[finish-1]]
                        mlab.plot3d(linex, liney, linez, tube_radius = 2, color =dis_color)
        #
                if key_data['Number_markers_on_plot'] > 0:
                    for i in range(0,len(output.x),1):
                        mlab.text3d(output.x[i],output.y[i],output.z[i], str(int(i)), color = (0.0,0.0,0.0), scale = 5.0)
        #
                for i in range(0,len(part_data.setp),1):
                    try:
                        colour = color_key.color_code[part_data.type_part[i]-1]    # color_key was defined at the beginning of plotting
                    except:
                        colour = (0,0,0)   # default to black if there's some error.
                    mlab.points3d(output.x[part_data.setp[i]-1], output.y[part_data.setp[i]-1], output.z[part_data.setp[i]-1], color= colour, mode = 'sphere', scale_factor = part_data.part_size[i]/key_data['scale_pixels/nm'], resolution = 64)
        #
                mlab.gcf().scene.parallel_projection = True
                mlab.gcf().scene.show_axes = True
                mlab.show()
                file_name = str(track) + '.tiff'
                mlab.view(aziz,elevation, roll = roll_value)
                mlab.savefig(file_name,size=fig_size)
                mlab.close()
