"""
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>

@author: Josiah Walker
"""
import numpy,random
from BlockSparseMatrix import BlockSparseMatrix
from BresenhamAlgorithms import BresenhamLine,BresenhamTriangle,BresenhamPolygon

#ranges all given in cm
SonarSensor = {"spread": 15.*numpy.pi/180., "range": 500., "phitfree": -0.3, "phitoccupied": 3.}

class GridMap:
    """
    Sparse gridmap for 2D mapping.
    """
    
    def __init__(self,scale=1.0):
        """
        @brief Initialise a sparse block grid-map with arc-based sensor updates.
        
        @param scale The multiplier to rescale from input units to map cell size.
        """
        self._scale = scale
        self._map = BlockSparseMatrix()
    
    def update(self,position,distance,sensorangle,sensor):
        """
        @brief Update the map with a sensor reading.
        
        @param position The robot's current position given as (x,y,theta) for hte robot's position and angle.
        @param distance The distance measurement from the sensor.
        @param sensorangle The current angle from the robot's forward direction to the sensor.
        @param sensor A dict holding sensor-specific hardware data (see SonarSensor in this file).
        """
        
        #generate the angle positions (change angleUpdates for more accurate approximation)
        angleUpdates = 4
        thetas = []
        for i in xrange(angleUpdates-1):
            thetas.append(position[2] + i*sensor["spread"]/angleUpdates - sensor["spread"]/2. + sensorangle)
        thetas.append(position[2] + sensor["spread"]/2. + sensorangle)
        
        #generate the arc and robot positions
        positions = [numpy.array(position[:2])*self._scale]
        for t in thetas:
            positions.append(
                numpy.round(
                    numpy.array([numpy.cos(t),numpy.sin(t)]) * 
                    distance *
                    self._scale + positions[0]
                ).astype(numpy.int64)
            )
        positions[0] = numpy.round(positions[0]).astype(numpy.int64)
        
        #FILL THE EMPTY ARC AREA OF THE SENSOR (as an approximate polygon)
        emptyVal = sensor["phitfree"]
        for cell in BresenhamPolygon(positions):
            self._map[cell[0],cell[1]] = max(emptyVal+self._map[cell[0],cell[1]],-20.) #clip to -20
        
        #DO BRESENHAM detection on the arc edge for object hits
        hitVals = BresenhamLine(positions[1],positions[2])
        solidVal = sensor["phitoccupied"]
        startpt = 0
        for i in xrange(1,len(positions)-1):
            hitVals = BresenhamLine(positions[i],positions[i+1])
            solidVal = sensor["phitoccupied"]
            for h in hitVals[startpt:]:
                self._map[h[0],h[1]] = min(solidVal+self._map[h[0],h[1]],120.) #clip to 120
            startpt = 1 #skip the first part of all following line segments
    
    def get(self,location):
        """
        @brief Get the value at a certain x,y location.
        
        @param location A location in the form [x,y]
        """
        location = numpy.round(location*self._scale).astype(numpy.int64)
        return self._map(location[0],location[1])
    
    def getRange(self,topleft,bottomright):
        """
        @brief Get the values for a range of locations as a matrix. Note: this returns at the internal scale, not the external scale
        
        @param topleft A location in the form [x,y] in external units designating the top left of the area
        @param bottomright A location in the form [x,y] in external units designating the bottom right of the area
        """
        #convert into map scale
        topleft = numpy.round(numpy.array(topleft)*self._scale).astype(numpy.int64)
        bottomright = numpy.round(numpy.array(bottomright)*self._scale).astype(numpy.int64)
        #fill in the output
        result = numpy.zeros((bottomright[0]-topleft[0],bottomright[1]-topleft[1]))
        for i in xrange(topleft[0],bottomright[0]):
            ival = numpy.round(i).astype(numpy.int64)
            for j in xrange(topleft[1],bottomright[1]):
                jval = numpy.round(j).astype(numpy.int64)
                result[i-topleft[0],j-topleft[1]] = self._map[ival,jval]
        return result




if __name__ == '__main__':
    """
    Do validation test
    """
    import time,os
    from matplotlib import pyplot
    
    #set this true and have mencoder to create a video of the test
    makevideo = True
    
    #set up the map and scale
    scale = 100.0
    groundtruth = ((1,1,1,1,1),
                   (1,0,0,0,1),
                   (1,0,1,0,1),
                   (1,0,0,0,1),
                   (1,1,1,1,1))
    gridScale = 0.5
    #set up the grid map on a 2cm scale (half the input resolution)
    estmap = GridMap(scale=gridScale)
    
    #this is the set of positions the rover moves between
    tour = ((150.0,150.0,0.0),(350.0,150.0,0.0),
            (350.0,150.0,numpy.pi/2.0),(350.0,350.0,numpy.pi/2.0),
            (350.0,350.0,numpy.pi),(150.0,350.0,numpy.pi),
            (150.0,350.0,numpy.pi*1.5),(150.0,150.0,numpy.pi*1.5),(150.0,150.0,numpy.pi*2))
    
    #this is the number of steps along each part of the tour
    divs =100
    vals = []
    for i in xrange(len(tour)-1):
        
        
        for j in xrange(divs):
            position = numpy.array(tour[i])*(1.-j/float(divs))+numpy.array(tour[(i+1)%len(tour)])*(j/float(divs))
            
            p = position[:2]
            a = -position[2]+numpy.pi
            offset = numpy.array([numpy.sin(a),numpy.cos(a)])*20.
            
            for k in xrange(4):
                
                #simulate each of the sonar sensor sweeps and see if we hit anything.
                sensor = SonarSensor
                sensorangle = numpy.pi/2*k
                thetamax = position[2] + sensor["spread"]/2. + sensorangle
                thetamin = position[2] - sensor["spread"]/2. + sensorangle
                baseB = numpy.array([numpy.cos(thetamax),numpy.sin(thetamax)])
                baseC = numpy.array([numpy.cos(thetamin),numpy.sin(thetamin)])
                hit = False
                for distance in xrange(int(sensor["range"])):
                    B = numpy.round(baseB*distance + position[:2]).astype(numpy.int32)
                    C = numpy.round(baseC*distance + position[:2]).astype(numpy.int32)
                    
                    for pos in BresenhamLine(B,C):
                        if groundtruth[int((pos[0]/scale))][int((pos[1]/scale))] == 1:
                            distance = numpy.linalg.norm(position[:2] - pos) #add noise in here if you want noise
                            hit = True
                            break
                    if hit:
                        t0 = time.time()
                        estmap.update(position,distance,sensorangle,sensor)
                        vals.append(time.time()-t0)
                        break
                if not hit:
                    t0 = time.time()
                    estmap.update(position,distance,sensorangle,sensor)
                    vals.append(time.time()-t0)
            
            if makevideo: #save out png's for the video
                fname = '_tmp%05d.png'%(i*divs+j)
                tl = (95,95)
                print (i*divs+j)
                robot = (numpy.array([p+offset,p-offset,p+numpy.array([-offset[1],offset[0]])])*gridScale-numpy.array(tl)*gridScale).astype(numpy.int64)
                emap = numpy.clip(estmap.getRange(tl,(405,405)), -1000,1000 )
                for cell in BresenhamTriangle(robot[0],robot[1],robot[2]):
                    emap[cell[0],cell[1]] = 120
                pyplot.imsave(fname,emap)
                pyplot.clf()
                    
    print "Mean Sensor Update Time:", numpy.mean(vals)
    
    if makevideo: #convert png's to video
        #recent ubuntu versions use avconv
        os.system("avconv -r 30 -i _tmp%05d.png -b:v 1000k rovertest.mp4")
        #os.system("mencoder 'mf://*.png' -mf type=png:fps=30 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o rovertest.avi")
        os.system("rm -f _tmp*.png")
    
