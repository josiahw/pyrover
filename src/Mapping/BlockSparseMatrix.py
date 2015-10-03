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

import numpy

class BlockSparseMatrix:
    
    def __init__(self,blocksize=(500,500),dtype=numpy.float32):
        self.blocksize=blocksize
        self.dtype=dtype
        self.blocks = {} #start with an empty matrix
    
    def __getitem__(self, key):
        #use a dictionary key to see if the part of the matrix we want is available
        k = (int(key[0]/self.blocksize[0]),int(key[1]/self.blocksize[1]))
        
        if not (k in self.blocks): #return 0 if uninitialised
            return 0
        else:
            k2 = (int(key[0]%self.blocksize[0]),int(key[1]%self.blocksize[1]))
            return self.blocks[k][k2]
    
    def __setitem__(self, key, item):'
        #use a dictionary key to see if the part of the matrix we want is available
        k = (int(key[0]/self.blocksize[0]),int(key[1]/self.blocksize[1]))
        
        if not (k in self.blocks): #create a new block if this address isn't in our map
            self.blocks[k] = numpy.zeros(self.blocksize,dtype=self.dtype)
        k2 = (int(key[0]%self.blocksize[0]),int(key[1]%self.blocksize[1]))
        self.blocks[k][k2] = item

