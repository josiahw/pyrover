# -*- coding: utf-8 -*-
"""
Created on Sat May 16 20:47:27 2015

@author: merlz
"""
import numpy,random
from BresenhamAlgorithms import BresenhamLine,BresenhamTriangle

#ranges all given in cm
SonarSensor = {"spread": 15.*numpy.pi/180., "range": 500., "phitfree": -0.2, "phitoccupied": 3.5}

class BlockSparseMatrix:
    
    def __init__(self,blocksize=(500,500),dtype=numpy.float32):
        self.blocksize=blocksize
        self.dtype=dtype
        self.blocks = {}
    
    def __getitem__(self, key):
        k = (int(key[0]/self.blocksize[0]),int(key[1]/self.blocksize[1]))
        if not (k in self.blocks):
            return 0
        else:
            k2 = (int(key[0]%self.blocksize[0]),int(key[1]%self.blocksize[1]))
            return self.blocks[k][k2]
    
    def __setitem__(self, key, item):
        k = (int(key[0]/self.blocksize[0]),int(key[1]/self.blocksize[1]))
        if not (k in self.blocks):
            self.blocks[k] = numpy.zeros(self.blocksize,dtype=self.dtype)
        k2 = (int(key[0]%self.blocksize[0]),int(key[1]%self.blocksize[1]))
        self.blocks[k][k2] = item

class FBMatrix:
    
    def __init__(self,order = 15,dtype=numpy.float32):
        self.ordvals = numpy.arange(0.,order+1.,1.).reshape((-1,1))
        self.scales = 1./numpy.array([max(item**2+sublist**2,1.)**0.5 for sublist in xrange(0,order+1) for item in xrange(0,order+1)])
        self.coeffs = numpy.zeros((order+1)**2,dtype=numpy.float64)
        self.dp = numpy.ones(self.ordvals.shape[0]).reshape((1,-1))
        self.coeff0 = (numpy.pi) * (1./5000.) #500 is the longest wave we represent
        #self.scales /= numpy.sum(self.scales)
        #this is an optimisation
        #self.ordvals *= self.coeff0
        self.mom = 0.0
    
    def __getitem__(self, key):
        x0 = numpy.dot(self.ordvals*key[0],self.dp)
        y0 = numpy.dot(self.ordvals*key[1],self.dp).T
        vals = numpy.cos((x0+y0)*self.coeff0).ravel()
        vals[0] = 0.
        return numpy.dot(self.coeffs,vals)
        
    
    def __setitem__(self, key, item):
        x0 = numpy.dot(self.ordvals*key[0],self.dp)
        y0 = numpy.dot(self.ordvals*key[1],self.dp).T
        vals = numpy.cos((x0+y0)*self.coeff0).ravel()
        vals[0] = 0.
        #print vals
        error = item - numpy.dot(self.coeffs,vals)
        #print item,numpy.dot(self.coeffs,vals),error
        error = 0.2 * error * self.scales * vals #/ numpy.sum(vals*self.scales)
        #print error
        #print error
        #print vals
        self.coeffs += error #/numpy.sum(error)
        

class ApproxMatrix:
    
    def __init__(self,dbsize=(4,512),dtype=numpy.float32):
        self.map = numpy.zeros(dbsize,dtype)
        self.offsets = numpy.square(numpy.array(range(1,2*dbsize[0]+1,2),dtype=numpy.uint32))
        self.coords1 = numpy.array(range(dbsize[0]))
    
    def __getitem__(self, key):
        random.seed((numpy.left_shift(numpy.int64(numpy.int32(key[0])),32) + numpy.int64(numpy.int32(key[1]))))
        k2 = numpy.array([random.randint(0,self.map.shape[1]-1) for i in xrange(self.map.shape[0])])
        return numpy.mean(self.map[self.coords1,k2])
    
    def __setitem__(self, key, item):
        random.seed((numpy.left_shift(numpy.int64(numpy.int32(key[0])),32) + numpy.int64(numpy.int32(key[1]))))
        k2 = numpy.array([random.randint(0,self.map.shape[1]-1) for i in xrange(self.map.shape[0])])
        self.map[self.coords1,k2] += item - numpy.mean(self.map[self.coords1,k2])


class GridMap:
    """
    Sparse gridmap for 2D mapping.
    """
    
    def __init__(self,scale=0.5):
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
        
        thetamax = position[2] + sensor["spread"]/2. + sensorangle
        thetamin = position[2] - sensor["spread"]/2. + sensorangle
        
        #for efficiency, we assume that the arc is a straight line, and the area is a triangle
        A = numpy.array(position[:2])*self._scale
        B = numpy.round(numpy.array([numpy.cos(thetamax),numpy.sin(thetamax)])*distance + A).astype(numpy.int64)
        C = numpy.round(numpy.array([numpy.cos(thetamin),numpy.sin(thetamin)])*distance + A).astype(numpy.int64)
        A = numpy.round(A).astype(numpy.int64)
        
        #FILL THE EMPTY ARC OF THE SENSOR (as an approximate triangle)
        emptyVal = sensor["phitfree"]
        for cell in BresenhamTriangle(A,B,C):
            self._map[cell[0],cell[1]] = max(emptyVal+self._map[cell[0],cell[1]],-20.)
        
        #DO BRESENHAM from B to C for hit objects
        hitVals = BresenhamLine(B,C)
        solidVal = sensor["phitoccupied"]
        for h in hitVals:
            self._map[h[0],h[1]] = min(solidVal+self._map[h[0],h[1]],120.)
        
    
    def get(self,location):
        """
        @brief Get the value at a certain x,y location.
        
        @param location A location in the form [x,y]
        """
        location = numpy.round(location).astype(numpy.int64)
        return self._map(location[0],location[1])
    
    def getRange(self,topleft,bottomright):
        """
        @brief Get the values for a range of locations as a matrix.
        
        @param topleft A location in the form [x,y] designating the top left of the area
        @param bottomright A location in the form [x,y] designating the bottom right of the area
        """
        #convert into map scale
        topleft = numpy.round(numpy.array(topleft)*self._scale).astype(numpy.int64)
        bottomright = numpy.round(numpy.array(bottomright)*self._scale).astype(numpy.int64)
        
        #fill in the output
        result = numpy.zeros((bottomright[0]-topleft[0],bottomright[1]-topleft[1]))
        for i in xrange(topleft[0],bottomright[0]):
            ival = numpy.round(i*self._scale).astype(numpy.int64)
            for j in xrange(topleft[1],bottomright[1]):
                jval = numpy.round(j*self._scale).astype(numpy.int64)
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
    
    #set up the grid map on a 2cm scale (half the input resolution)
    estmap = GridMap(scale=0.5)
    
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
                            distance = numpy.linalg.norm(position - pos)
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
                print (i*divs+j)
                pyplot.imsave(fname,numpy.clip(estmap.getRange((95,95),(405,405)), -1000,1000 ))
                pyplot.clf()
                    
    print "Mean Sensor Update Time:", numpy.mean(vals)
    
    if makevideo: #convert png's to video
        os.system("mencoder 'mf://*.png' -mf type=png:fps=30 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o rovertest.avi")
        os.system("rm -f _tmp*.png")
    
