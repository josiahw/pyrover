import numpy

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

