# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 12:23:58 2015

@author: merlz
"""
import numpy

def BresenhamLine(A,B):
    """
    This implements the bresenham algorithm to draw a line from A to B
    
    Returns: all x,y index pairs that are in the line
    """
    if A[0] == B[0]: #if A and B share the same X coord, draw a straight line in Y
        increment =  1 if B[1] >= A[1] else -1 #decide whether to draw forwards or backwards
        return [(A[0],i) for i in xrange(A[1],B[1]+increment,increment)]
            
    elif A[1] == B[1]: #if A and B share the same Y coord, draw a straight line in X
        increment =  1 if B[0] >= A[0] else -1 #decide whether to draw forwards or backwards
        return [(i,A[1]) for i in xrange(A[0],B[0]+increment,increment)]
    
    else: #this draws a diagonal line
        incrementx = 1 if B[0] >= A[0] else -1 #set the direction of line drawing
        incrementy = 1 if B[1] >= A[1] else -1
        
        result = []
        yval = A[1] #set Y start
        slope = A-B #calculate slope - assuming A and B are numpy arrays
        slope = abs(slope[1]/float(slope[0]))
        error = 0.0 #initialise Y error counter
        for i in xrange(A[0],B[0]+incrementx,incrementx): #for all X values, step forward in X
            result.append((i,yval))
            error += slope
            while error >= 0.5: #while the Y error is too large, step forward in Y
                yval += incrementy
                error -= 1.0
                result.append((i,yval))
        return result

def BresenhamBorder(A,B):
    """
    Unlike the line, this does only one pixel per Y row, so it can be used in fill algorithms efficiently
    
    Returns: all x,y index pairs that are in the border
    """
    if A[0] == B[0]: #if A and B share the same X coord, draw a straight line in Y
        increment =  1 if B[1] >= A[1] else -1
        return [(A[0],i) for i in xrange(A[1],B[1]+increment,increment)]
        
    elif A[1] == B[1]: #we're screwed - we can only return one Y value
        return [numpy.round((A+B)/2).astype(numpy.int32)]
        
    else: #work out what to do for a diagonal
        incrementy = 1 if B[1] >= A[1] else -1 #set the direction of line drawing
        
        slope = A-B
        slope = slope[0]/float(slope[1])*incrementy
        xvals = numpy.round(A[0] + slope*numpy.arange(0.,abs(A[1]-B[1])+1,1.)).astype(numpy.int32)
        return [(xvals[i], y) for i,y in enumerate(xrange(A[1],B[1]+incrementy,incrementy))]

def BresenhamPolygon(vertices):
    #put the largest value at the head of the list:
    maxvert = 0
    for i in xrange(len(vertices)):
        if vertices[i][1] > vertices[maxvert][1]:
            maxvert = i
    
    vertices = vertices[maxvert:] + vertices[:maxvert]
    #split the list in to two sides based on max->min paths
    minvert = 0
    for i in xrange(len(vertices)):
        if vertices[i][1] < vertices[minvert][1]:
            minvert = i
    #skip everything of equal Y height on the top
    start = 0
    while start < len(vertices)-2 and vertices[start][1] == vertices[start+1][1]:
        start += 1
    side1 = vertices[start:minvert+1]
    #create the "left" border
    l = BresenhamBorder(side1[0],side1[1])
    for i in xrange(1,len(side1)-1):
        l += BresenhamBorder(side1[i],side1[i+1])[1:]
    
    #skip everything of equal Y height on the bottom
    while minvert < len(vertices)-2 and vertices[minvert][1] == vertices[minvert+1][1]:
        minvert += 1
    side2 = vertices[minvert:]
    side2.reverse()
    side2 = [vertices[0]] + side2
    #create the "right" border
    r = BresenhamBorder(side2[0],side2[1])
    for i in xrange(1,len(side2)-1):
        r += BresenhamBorder(side2[i],side2[i+1])[1:]
    
    #do horizontal scans and save all the cell locations in the triangle
    result = []
    for i in xrange(len(l)):
        increment = 1 if r[i][0] >= l[i][0] else -1 
        result += [(j,l[i][1]) for j in xrange(l[i][0],r[i][0]+increment,increment)]
    return result

def BresenhamTriangle(A,B,C):
    #this is here because not all the functions have been upgraded to polygon yet
    return BresenhamPolygon([A,B,C])
