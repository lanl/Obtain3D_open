import math
import os

########### starting working directory ########
os.chdir("C:\\Users\\252321\\Desktop\\000_Obtain3D_20220406\\")
check_01 = 1
########################################
### choosing the folder to look for the files in ##
########################################
while check_01 > -0.5:
    aa = os.listdir(os.getcwd())
    for i in range(0,len(aa),1):
        print(str(i+1) + '.) ' + str(aa[i])) 
    print(str(len(aa)+1) + '.) ' + 'Directory up') 
    dir_change_01 = int(input('Choose number of folder? Choose "0" if in directory with desired files' '\n'))-1
    check_01 = dir_change_01
    if check_01 > -0.5 and check_01 < len(aa):
        os.chdir(aa[dir_change_01])
        print(os.getcwd())
    if check_01 > len(aa) -0.5:
        os.chdir("..")
        print(os.getcwd())
print('\n')

##################################### Functions ################################
def rot(l,m,n,ang,x,y,z):
    l_unit = l/(l**2+m**2+n**2)**0.5
    n_unit = n/(l**2+m**2+n**2)**0.5
    m_unit = m/(l**2+m**2+n**2)**0.5

    c = math.cos(ang)
    s = math.sin(ang)

    a11 = l_unit*l_unit*(1-c) + c
    a12 = m_unit*l_unit*(1-c) - n_unit*s
    a13 = n_unit*l_unit*(1-c)+m_unit*s
    a21 = l_unit*m_unit*(1-c)+n_unit*s
    a22 = m_unit*m_unit*(1-c)+c
    a23 = n_unit*m_unit*(1-c)-l_unit*s
    a31 = l_unit*n_unit*(1-c)-m_unit*s
    a32 = m_unit*n_unit*(1-c)+l_unit*s
    a33 = n_unit*n_unit*(1-c)+c
    cord = [0,0,0]
    cord[0] = a11*x+a12*y+a13*z
    cord[1] = a21*x+a22*y+a23*z
    cord[2] = a31*x+a32*y+a33*z
    return cord

def cross_product(x1,y1,z1,x2,y2,z2):
    cp = [y1*z2-y2*z1, z1*x2-z2*x1, x1*y2-x2*y1]
    return cp
    
def unit_vector(h,k,l):
    unit = [h/(h**2+k**2+l**2)**0.5, k/(h**2+k**2+l**2)**0.5, l/(h**2+k**2+l**2)**0.5]
    return unit

def mean(set):
    mean_result = sum(set)/len(set)
    return mean_result

def stdv(set):
    std = 5000
    if len(set) > 1:
        su = []
        mean = sum(set)/len(set)
        for i in range(0,len(set),1):
            n = (set[i]-mean)**2
            su.append(n)
        sume = sum(su)
        std = ((1.00/(len(set)-1))*sume)**0.5
    return std

def angle_two_vectors(h1,k1,l1,h2,k2,l2):
    a_dot_b = (h1*h2+k1*k2+l1*l2)
    mag_a = (h1**2+k1**2+l1**2)**0.5
    mag_b = (h2**2+k2**2+l2**2)**0.5
    cos_theta = a_dot_b/(mag_a*mag_b+0.000000000000000000000000000000001)
    theta = math.acos(cos_theta)
    return theta

###########################################################################
#################################### setting up variables ######################
###########################################################################
positions = []  
axis = []              
x1 = []       
y1 = []     
x2 = []      
y2 = []
dis_number = []
point_number = []
position_x = []
position_y = []
position_z = []
data = []
################################################################################
##################################### opening files in working directory ###############
################################################################################

with open('in.txt', 'r') as f:   ### opens the text file with the input data
    for line in f:
        data = line.split()     ### 
        x1.append(float(data[0]))  ### position of x in pixels for the first image. Positive x is to the right
        y1.append(float(data[1]))  ### position of y in pixels for the first image. Positive y is down
        x2.append(float(data[2]))  ### position of x in pixels for the second image. Positive x is to the right
        y2.append(float(data[3]))  ### position of y in pixels for the second image. Positive y is down
        dis_number.append(float(data[4])) ### determines whether point is part of a dislocation and which one, or if point is a sphere
        point_number.append(float(data[5])) ### is the point number

with open('key.txt','r') as f:
    data = []
    for line in f:
        data_1 = line.split()
        data.append(data_1)       
xtilt1 = float(data[0][1])
xtilt2 = float(data[1][1])
ytilt1 = float(data[2][1])
ytilt2 = float(data[3][1])
scale = float(data[4][1])
surfaces = float(data[5][1])
particles = float(data[6][1])
view = float(data[7][1])
if view > 0.5:
    import mayavi.mlab as mlab
number_markers = float(data[8][1])
movie = float(data[9][1])

thetat = xtilt2 - xtilt1  #### Calculates the difference in x-tilt between the two images
thetatr = math.radians(thetat)  ### Converts difference in tilt to radians 

segs = []
segs_dis = []
for j in range(0,5100,1):
    count = j
    groups = []    
    for i in range(0,len(dis_number),1):
        if abs(count-dis_number[i]) < 0.1:
            groups.append(point_number[i])
    if sum(groups) > 0:
        groups.append(count)   
        segs.append(groups)
        segs_dis.append(groups)

dis_segs = []    
for i in range(0,len(segs),1):      ################# dis_segs splits dislocations into 2 point segments for plotting
    for j in range(0,len(segs[i])-2,1):
        dis_segs.append([segs[i][j],segs[i][j+1],segs[i][-1]])

setp = []
part_size = []
type_part = []
if particles > 0:
    with open('part.txt', 'r') as f:
        for line in f:
            data = line.split()
            setp.append(int(data[0]))
            part_size.append(float(data[1]))
            type_part.append(float(data[2]))

########################################################################################################################

points = int(len(x1))
reference_number = 0         # chooses the point all other points will be referenced to 
refx1 = x1[reference_number]   
refy1 = y1[reference_number]
refx2 = x2[reference_number]   
refy2 = y2[reference_number]

for i in range(0,points,1):  ### for loop to define each point with respect to the reference points just defined
    positions.append([x1[i]-refx1, (refy1-y1[i])*(1), x2[i]-refx2, (refy2-y2[i])*(1)])
    axisangle = math.atan2(positions[i][1]-positions[i][3], positions[i][2]-positions[i][0])  #### atan2 gives angle in range between -pi and pi, with respect to the positive x-axis, this is used to find the tilt axis wrt the images
    if i > 0:    
        if axisangle > 0:
            axis.append(axisangle)
        else:
            axis.append(axisangle+3.1415926) ### adding pi to the negative values places it on the opposite side of the origin, keeping the same tilt axis but keeping all the values positive for averaging        

refaxis = []
refaxis2 = []

axisavg = mean(axis)
axisdev = 50

for i in range(0,len(axis),1):  ### loop removes values for the tilt axis more than a standard deviation away from the mean
    ax = axis[i]
    if ax > axisavg - axisdev and ax < axisavg + axisdev:
        refaxis.append(ax)
axisavg = mean(refaxis)
axisdev = stdv(refaxis) 

for i in range(0,len(refaxis),1):  
    ax = refaxis[i]
    if ax > axisavg - axisdev and ax < axisavg + axisdev:
        refaxis2.append(ax)
axisavg = mean(refaxis2)

axisavg = math.degrees(axisavg)
print("shift axis direction in degrees")
print(axisavg)   ### prints the value of the tilt axis in degrees
thetarot1 = (-1)*(90 - axisavg)   ### used in first rotation below

################## z is calculated for each point seperately
angle_0 = math.radians(thetarot1)
axis_0 = [0,0,1]
for k in range(0,points,1):
    e = positions[k]

# Rotations place the direction of the separation change in the y-axis for the calculation   
    xy1 = rot(axis_0[0],axis_0[1],axis_0[2], angle_0 ,e[0],e[1],0)
    xy2 = rot(axis_0[0],axis_0[1],axis_0[2], angle_0 ,e[2],e[3],0)

##################################

    delth = (xy2[1] - xy1[1])/(2.0*math.sin(thetatr/2.0))
    z1 = (-1.0)*(math.tan(thetatr/2.0)*xy1[1] + delth/(math.sin(3.14159265/2-thetatr/2.0))) 

    angle = (-1)*math.radians(xtilt1)         ###################### Rotate about x tilt axis
    xyz0 = rot(1,0,0,angle,xy1[0],xy1[1],z1)    

    angle = math.radians(ytilt1)         ###################### Rotate about y tilt axis
    xyz00 = rot(0,1,0,angle,xyz0[0],xyz0[1],xyz0[2])

    position_x.append(xyz00[0]/scale)
    position_y.append(xyz00[1]/scale)
    position_z.append(xyz00[2]/scale)

with open('coordinates.txt', 'w') as f:
    header = str("number" "\t" "x" "\t" "y" "\t" "z")
    f.write(header)
    f.write('\n')

    for i in range(0,len(position_x),1):
        f.write(str(i+1))
        f.write("\t") 
        f.write(str(position_x[i]))
        f.write("\t")
        f.write(str(position_y[i]))
        f.write("\t")
        f.write(str(position_z[i]))
        f.write("\t")
        f.write('\n')

######################## Plotting (plot 1) z is in beam direction ########################################################
dis_color = (0,0,0) ### color for the dislocation segments
if view > 0.5 and movie < 0.5:
    outside = 10
    x_box_min = min(position_x)-outside
    x_box_max = max(position_x)+outside
    y_box_min = min(position_y)-outside
    y_box_max = max(position_y)+outside
    z_box_min = min(position_z)-outside
    z_box_max = max(position_z)+outside
    tube = 0.5
    
    mlab.figure(bgcolor=(1,1,1))

    
    linex_b = [x_box_min, x_box_max]
    liney_b = [y_box_min, y_box_min]
    linez_b = [z_box_min, z_box_min]
    mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
    linex_b = [x_box_min, x_box_min]
    liney_b = [y_box_min, y_box_max]
    linez_b = [z_box_min, z_box_min]
    mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
    linex_b = [x_box_min, x_box_max]
    liney_b = [y_box_max, y_box_max]
    linez_b = [z_box_min, z_box_min]
    mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
    linex_b = [x_box_max, x_box_max]
    liney_b = [y_box_min, y_box_max]
    linez_b = [z_box_min, z_box_min]
    mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
    linex_b = [x_box_min, x_box_min]
    liney_b = [y_box_min, y_box_min]
    linez_b = [z_box_min, z_box_max]
    mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
    linex_b = [x_box_max, x_box_max]
    liney_b = [y_box_min, y_box_min]
    linez_b = [z_box_min, z_box_max]
    mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
    linex_b = [x_box_max, x_box_max]
    liney_b = [y_box_max, y_box_max]
    linez_b = [z_box_min, z_box_max]
    mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
    linex_b = [x_box_min, x_box_min]
    liney_b = [y_box_max, y_box_max]
    linez_b = [z_box_min, z_box_max]
    mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
    linex_b = [x_box_min, x_box_max]
    liney_b = [y_box_min, y_box_min]
    linez_b = [z_box_max, z_box_max]
    mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
    linex_b = [x_box_min, x_box_min]
    liney_b = [y_box_min, y_box_max]
    linez_b = [z_box_max, z_box_max]
    mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
    linex_b = [x_box_min, x_box_max]
    liney_b = [y_box_max, y_box_max]
    linez_b = [z_box_max, z_box_max]
    mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
    linex_b = [x_box_max, x_box_max]
    liney_b = [y_box_min, y_box_max]
    linez_b = [z_box_max, z_box_max]
    mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
    
    linex = []
    liney = []
    linez = []
    for i in range(0,len(dis_segs),1):
        if dis_segs[i][-1] > 0 and dis_segs[i][-1] <5000:     
            star = int(dis_segs[i][0])
            finish = int(dis_segs[i][1])
            linex = [position_x[star-1], position_x[finish-1]]
            liney = [position_y[star-1], position_y[finish-1]]
            linez = [position_z[star-1], position_z[finish-1]]
            mlab.plot3d(linex, liney, linez, tube_radius = 2, color = dis_color)
    
    if number_markers > 0:                  
        for i in range(0,len(position_x),1):
            fig_text = str(int(i))
            position = [position_x[i],position_y[i],position_z[i]]
            mlab.text3d(position[0], position[1], position[2], fig_text, color = (0.0,0.0,0.0), scale = 5.0)
        
    for i in range(0,len(setp),1):
        part_plot = setp[i]
        size_p = part_size[i]
        color_num = type_part[i]
        if abs(color_num-1) < 0.3:
            colour = (0,1,0)  #### green
        if abs(color_num-2) < 0.3:
            colour = (0,0,0)  #### black
        if abs(color_num-3) < 0.3:
            colour = (1,0,1) #### purple
        if abs(color_num-4) < 0.3:
            colour = (1,0,0) #### red
        if abs(color_num-5) < 0.3:
            colour = (1,1,0) #### yellow   
        linex = [position_x[part_plot-1]]
        liney = [position_y[part_plot-1]]
        linez = [position_z[part_plot-1]]
        mlab.points3d(linex, liney, linez, color= colour, mode = 'sphere', scale_factor = size_p/scale, resolution = 64)
    
    mlab.gcf().scene.parallel_projection = True
    mlab.gcf().scene.show_axes = True

#################### Movie Mode ##########################

if movie > 0.1 and view > 0.5:

    outside = 10
    x_box_min = min(position_x)-outside
    x_box_max = max(position_x)+outside
    y_box_min = min(position_y)-outside
    y_box_max = max(position_y)+outside
    z_box_min = min(position_z)-outside
    z_box_max = max(position_z)+outside
    tube = 0.5
    
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
                
                linex_b = [x_box_min, x_box_max]
                liney_b = [y_box_min, y_box_min]
                linez_b = [z_box_min, z_box_min]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_min]
                liney_b = [y_box_min, y_box_max]
                linez_b = [z_box_min, z_box_min]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_max]
                liney_b = [y_box_max, y_box_max]
                linez_b = [z_box_min, z_box_min]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_max, x_box_max]
                liney_b = [y_box_min, y_box_max]
                linez_b = [z_box_min, z_box_min]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_min]
                liney_b = [y_box_min, y_box_min]
                linez_b = [z_box_min, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_max, x_box_max]
                liney_b = [y_box_min, y_box_min]
                linez_b = [z_box_min, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_max, x_box_max]
                liney_b = [y_box_max, y_box_max]
                linez_b = [z_box_min, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_min]
                liney_b = [y_box_max, y_box_max]
                linez_b = [z_box_min, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_max]
                liney_b = [y_box_min, y_box_min]
                linez_b = [z_box_max, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_min]
                liney_b = [y_box_min, y_box_max]
                linez_b = [z_box_max, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_max]
                liney_b = [y_box_max, y_box_max]
                linez_b = [z_box_max, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_max, x_box_max]
                liney_b = [y_box_min, y_box_max]
                linez_b = [z_box_max, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                
                linex = []
                liney = []
                linez = []
                for i in range(0,len(dis_segs),1):
                    if dis_segs[i][-1] > 0 and dis_segs[i][-1] <5000:     
                        star = int(dis_segs[i][0])
                        finish = int(dis_segs[i][1])
                        linex = [position_x[star-1], position_x[finish-1]]
                        liney = [position_y[star-1], position_y[finish-1]]
                        linez = [position_z[star-1], position_z[finish-1]]
                        mlab.plot3d(linex, liney, linez, tube_radius = 2, color =dis_color)
                
                if number_markers > 0:                  
                    for i in range(0,len(position_x),1):
                        fig_text = str(int(i))
                        position = [position_x[i],position_y[i],position_z[i]]
                        mlab.text3d(position[0], position[1], position[2], fig_text, color = (0.0,0.0,0.0), scale = 5.0)
                    
                for i in range(0,len(setp),1):
                    part_plot = setp[i]
                    size_p = part_size[i]
                    color_num = type_part[i]
                    if abs(color_num-1) < 0.3:
                        colour = (0,1,0)  #### green
                    if abs(color_num-2) < 0.3:
                        colour = (0,0,0)  #### black
                    if abs(color_num-3) < 0.3:
                        colour = (1,0,1) #### purple
                    if abs(color_num-4) < 0.3:
                        colour = (1,0,0) #### red
                    if abs(color_num-5) < 0.3:
                        colour = (1,1,0) #### yellow   
                    linex = [position_x[part_plot-1]]
                    liney = [position_y[part_plot-1]]
                    linez = [position_z[part_plot-1]]
                    mlab.points3d(linex, liney, linez, color= colour, mode = 'sphere', scale_factor = size_p/scale, resolution = 64)
            
                mlab.gcf().scene.parallel_projection = True
                mlab.gcf().scene.show_axes = True
                mlab.show()
                file_name = str(track) + '.tiff'
                mlab.view(aziz,elevation, roll = roll_value)
                mlab.savefig(file_name,size=fig_size)
                mlab.close()
                
        if abs(elevation_end[b]-elevation_start[b]) > 1:
            print("elevation")
            for x in range(elevation_start[b],elevation_end[b],angle_step[b]):
                track = track + 1
                mlab.figure(bgcolor=(1,1,1)) 
                aziz = aziz_start[b]
                elevation = x
                roll_value = roll_start[b]
        
                mlab.close()
                mlab.figure(bgcolor=(1,1,1))
                fig_size = (600,600)
                
                linex_b = [x_box_min, x_box_max]
                liney_b = [y_box_min, y_box_min]
                linez_b = [z_box_min, z_box_min]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_min]
                liney_b = [y_box_min, y_box_max]
                linez_b = [z_box_min, z_box_min]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_max]
                liney_b = [y_box_max, y_box_max]
                linez_b = [z_box_min, z_box_min]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_max, x_box_max]
                liney_b = [y_box_min, y_box_max]
                linez_b = [z_box_min, z_box_min]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_min]
                liney_b = [y_box_min, y_box_min]
                linez_b = [z_box_min, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_max, x_box_max]
                liney_b = [y_box_min, y_box_min]
                linez_b = [z_box_min, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_max, x_box_max]
                liney_b = [y_box_max, y_box_max]
                linez_b = [z_box_min, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_min]
                liney_b = [y_box_max, y_box_max]
                linez_b = [z_box_min, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_max]
                liney_b = [y_box_min, y_box_min]
                linez_b = [z_box_max, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_min]
                liney_b = [y_box_min, y_box_max]
                linez_b = [z_box_max, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_max]
                liney_b = [y_box_max, y_box_max]
                linez_b = [z_box_max, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_max, x_box_max]
                liney_b = [y_box_min, y_box_max]
                linez_b = [z_box_max, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                
                linex = []
                liney = []
                linez = []
                for i in range(0,len(dis_segs),1):
                    if dis_segs[i][-1] > 0 and dis_segs[i][-1] <5000:     
                        star = int(dis_segs[i][0])
                        finish = int(dis_segs[i][1])
                        linex = [position_x[star-1], position_x[finish-1]]
                        liney = [position_y[star-1], position_y[finish-1]]
                        linez = [position_z[star-1], position_z[finish-1]]
                        mlab.plot3d(linex, liney, linez, tube_radius = 2, color = dis_color)
                
                if number_markers > 0:                  
                    for i in range(0,len(position_x),1):
                        fig_text = str(int(i))
                        position = [position_x[i],position_y[i],position_z[i]]
                        mlab.text3d(position[0], position[1], position[2], fig_text, color = (0.0,0.0,0.0), scale = 5.0)
                    
                for i in range(0,len(setp),1):
                    part_plot = setp[i]
                    size_p = part_size[i]
                    color_num = type_part[i]
                    if abs(color_num-1) < 0.3:
                        colour = (0,1,0)  #### green
                    if abs(color_num-2) < 0.3:
                        colour = (0,0,0)  #### black
                    if abs(color_num-3) < 0.3:
                        colour = (1,0,1) #### purple
                    if abs(color_num-4) < 0.3:
                        colour = (1,0,0) #### red
                    if abs(color_num-5) < 0.3:
                        colour = (1,1,0) #### yellow   
                    linex = [position_x[part_plot-1]]
                    liney = [position_y[part_plot-1]]
                    linez = [position_z[part_plot-1]]
                    mlab.points3d(linex, liney, linez, color= colour, mode = 'sphere', scale_factor = size_p/scale, resolution = 64)
            
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

                mlab.close()
                mlab.figure(bgcolor=(1,1,1))
                fig_size = (600,600)
                
                linex_b = [x_box_min, x_box_max]
                liney_b = [y_box_min, y_box_min]
                linez_b = [z_box_min, z_box_min]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_min]
                liney_b = [y_box_min, y_box_max]
                linez_b = [z_box_min, z_box_min]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_max]
                liney_b = [y_box_max, y_box_max]
                linez_b = [z_box_min, z_box_min]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_max, x_box_max]
                liney_b = [y_box_min, y_box_max]
                linez_b = [z_box_min, z_box_min]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_min]
                liney_b = [y_box_min, y_box_min]
                linez_b = [z_box_min, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_max, x_box_max]
                liney_b = [y_box_min, y_box_min]
                linez_b = [z_box_min, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_max, x_box_max]
                liney_b = [y_box_max, y_box_max]
                linez_b = [z_box_min, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_min]
                liney_b = [y_box_max, y_box_max]
                linez_b = [z_box_min, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_max]
                liney_b = [y_box_min, y_box_min]
                linez_b = [z_box_max, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_min]
                liney_b = [y_box_min, y_box_max]
                linez_b = [z_box_max, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_min, x_box_max]
                liney_b = [y_box_max, y_box_max]
                linez_b = [z_box_max, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                linex_b = [x_box_max, x_box_max]
                liney_b = [y_box_min, y_box_max]
                linez_b = [z_box_max, z_box_max]
                mlab.plot3d(linex_b, liney_b, linez_b, tube_radius = tube, color = (0,0,0))
                
                linex = []
                liney = []
                linez = []
                for i in range(0,len(dis_segs),1):
                    if dis_segs[i][-1] > 0 and dis_segs[i][-1] <5000:     
                        star = int(dis_segs[i][0])
                        finish = int(dis_segs[i][1])
                        linex = [position_x[star-1], position_x[finish-1]]
                        liney = [position_y[star-1], position_y[finish-1]]
                        linez = [position_z[star-1], position_z[finish-1]]
                        mlab.plot3d(linex, liney, linez, tube_radius = 2, color =dis_color)
                
                if number_markers > 0:                  
                    for i in range(0,len(position_x),1):
                        fig_text = str(int(i))
                        position = [position_x[i],position_y[i],position_z[i]]
                        mlab.text3d(position[0], position[1], position[2], fig_text, color = (0.0,0.0,0.0), scale = 5.0)
                    
                for i in range(0,len(setp),1):
                    part_plot = setp[i]
                    size_p = part_size[i]
                    color_num = type_part[i]
                    if abs(color_num-1) < 0.3:
                        colour = (0,1,0)  #### green
                    if abs(color_num-2) < 0.3:
                        colour = (0,0,0)  #### black
                    if abs(color_num-3) < 0.3:
                        colour = (1,0,1) #### purple
                    if abs(color_num-4) < 0.3:
                        colour = (1,0,0) #### red
                    if abs(color_num-5) < 0.3:
                        colour = (1,1,0) #### yellow   
                    linex = [position_x[part_plot-1]]
                    liney = [position_y[part_plot-1]]
                    linez = [position_z[part_plot-1]]
                    mlab.points3d(linex, liney, linez, color= colour, mode = 'sphere', scale_factor = size_p/scale, resolution = 64)
            
                mlab.gcf().scene.parallel_projection = True
                mlab.gcf().scene.show_axes = True
                mlab.show()
                file_name = str(track) + '.tiff'
                mlab.view(aziz,elevation, roll = roll_value)
                mlab.savefig(file_name,size=fig_size)
                mlab.close()            

        
        
