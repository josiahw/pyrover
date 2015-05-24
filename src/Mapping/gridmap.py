# -*- coding: utf-8 -*-
"""
Created on Sat May 16 20:47:27 2015

@author: merlz
"""
import numpy

#ranges all given in cm
SonarSensor = {"spread": 15.*numpy.pi/180., "range": 500., "phitfree": -0.2, "phitoccupied": 3.5}


def BresenhamLine(A,B):
    if A[0] == B[0]:
        if B[1] >= A[1]:
            return [(A[0],i) for i in xrange(A[1],B[1]+1,1)]
        else:
            return [(A[0],i) for i in xrange(A[1],B[1]-1,-1)]
    flipped = False
    if (A[0] > B[0]):
        flipped = True
        tmp = B
        B = A
        A = tmp
    result = []
    inc = numpy.sign(B[1]-A[1])
    error = 0
    yval = A[1]
    slope = A-B
    slope = abs(slope[1]/float(slope[0]))
    for i in xrange(A[0],B[0]+1):
        result.append((i,yval))
        error += slope
        while error >= 0.5:
            yval += inc
            error -= 1.0
            result.append((i,yval))
    if flipped:
        result = result[::-1]
    return result

#Unlike the line, this does only one pixel per Y row, so it can be used in fill algorithms efficiently
def BresenhamBorder(A,B):
    flipped = False
    if (A[1] > B[1]):
        flipped = True
        tmp = B
        B = A
        A = tmp
    result = []
    xval = A[0]
    slope = A-B
    slope = slope[0]/float(slope[1]+0.000000001)
    for i in xrange(A[1],B[1]+1):
        result.append((int(numpy.round(xval)),i))
        xval += slope
    if flipped:
        result = result[::-1]
    return result

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

class GridMap:
    """
    Generic sparse gridmap.
    """
    
    def __init__(self,size=(100000,100000)):
        self._map = BlockSparseMatrix() #((100000,100000),dtype=numpy.float32)
    
    def update(self,position,distance,sensorangle,sensor):
        """
        Position is given as (x,y,theta); the sensor field is filled by one of the defines in this file, or similar
        """
        thetamax = position[2] + sensor["spread"]/2. + sensorangle
        thetamin = position[2] - sensor["spread"]/2. + sensorangle
        
        #for efficiency, we assume that the object is a straight line, and the area is a triangle
        A = numpy.array(position[:2])
        B = numpy.round(numpy.array([numpy.cos(thetamax),numpy.sin(thetamax)])*distance + A).astype(numpy.int32)
        C = numpy.round(numpy.array([numpy.cos(thetamin),numpy.sin(thetamin)])*distance + A).astype(numpy.int32)
        A = numpy.round(A).astype(numpy.int32)
        coords = [A,B,C]
        
        #DO FILL OVER THE TRIANGLE ABC
        #1. sort in Y
        if A[1] < B[1]:
            tmp = A
            A = B
            B = tmp
        if A[1] < C[1]:
            tmp = A
            A = C
            C = tmp
        if B[1] < C[1]:
            tmp = C
            C = B
            B = tmp
        
        #2. generate borders
        emptyVal = sensor["phitfree"]
        if A[0] == B[0]:
            beginX = BresenhamBorder(A,C)
            endX = BresenhamBorder(B,C)
        else:
            beginX = BresenhamBorder(A,B)
            endX = BresenhamBorder(A,C)
            if len(beginX) < len(endX):
                beginX += BresenhamBorder(B,C)[1:]
            elif len(endX) < len(beginX):
                endX += BresenhamBorder(C,B)[1:]
        #3. fill the triangle
        for i in xrange(len(beginX)):
            if endX[i][0] >= beginX[i][0]:
                for j in xrange(beginX[i][0],endX[i][0]+1,1):
                    self._map[j,beginX[i][1]] += emptyVal
            else:
                for j in xrange(beginX[i][0],endX[i][0]-1,-1):
                    self._map[j,beginX[i][1]] += emptyVal
                
        
        #DO BRESENHAM from B to C for hit objects
        A,B,C = coords
        hitVals = BresenhamLine(B,C)
        solidVal = sensor["phitoccupied"]
        for h in hitVals:
            self._map[h[0],h[1]] += solidVal
        
    
    def get(self,location):
        """
        Locations are [x,y]
        """
        location = numpy.round(location).astype(numpy.int32)
        return self._map(location[0],location[1])
    
    def getRange(self,topleft,bottomright):
        """
        Locations are [x,y]
        """
        result = numpy.zeros((bottomright[0]-topleft[0],bottomright[1]-topleft[1]))
        for i in xrange(topleft[0],bottomright[0]):
            for j in xrange(topleft[1],bottomright[1]):
                result[i-topleft[0],j-topleft[1]] = self._map[i,j]
        return result




if __name__ == '__main__':
    """
    Do validation test
    """
    import time,os
    from matplotlib import pyplot
    makevideo = False
    #set up the map and scale
    scale = 100.0
    groundtruth = ((1,1,1,1,1),
                   (1,0,0,0,1),
                   (1,0,1,0,1),
                   (1,0,0,0,1),
                   (1,1,1,1,1))
    estmap = GridMap()
    
    #this is the set of positions the rover moves between
    tour = ((150.0,150.0,0.0),(350.0,150.0,0.0),
            (350.0,150.0,-numpy.pi/2.0),(350.0,350.0,-numpy.pi/2.0),
            (350.0,350.0,-numpy.pi),(150.0,350.0,-numpy.pi),
            (150.0,350.0,-numpy.pi*1.5),(150.0,150.0,-numpy.pi*1.5),(150.0,150.0,-numpy.pi*2))
    
    divs =100
    vals = []
    for i in xrange(len(tour)-1):
        
        
        for j in xrange(divs):
            position = numpy.array(tour[i])*(1.-j/float(divs))+numpy.array(tour[(i+1)%len(tour)])*(j/float(divs))
            #position[2] += numpy.sin(j/10.)
            #get range
            for k in xrange(4):
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
            if makevideo:
                fname = '_tmp%05d.png'%(i*divs+j)
                print (i*divs+j)
                pyplot.imsave(fname,numpy.clip(estmap.getRange((90,90),(410,410)),-20.,100.))
                pyplot.clf()
                    
    print "Mean Sensor Update Time:", numpy.mean(vals)
    
    if makevideo:
        os.system("mencoder 'mf://*.png' -mf type=png:fps=30 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o rovertest.avi")
        os.system("rm -f _tmp*.png")
    