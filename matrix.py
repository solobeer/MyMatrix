import numpy as np
import time

#numpy实现
'''
mat1 = np.ones((128, 128),dtype=np.float)
mat2 = np.ones((128, 128),dtype=np.float)
mat1 = np.matrix(mat1)
mat2 = np.matrix(mat2)
start = time.time()
mat3 = mat1 @ mat2
end = time.time()
print("time cost:", float(end - start) * 1000.0, "ms")
'''


#python实现

'''
class MatrixPy():
    
    def __init__(self, rows, cols, val) -> None:
        self._rows = rows
        self._cols = cols
        self._val = val
        self.data=[]
        for i in range(self._rows):
            self.data.append([])
            for j in range(self._cols):
                self.data[i].append(float(self._val))
    
    def mul(self, mat2):
        res = []
        for i in range(self._rows):
            res.append([0.0] * mat2._cols)
            for j in range(mat2._cols):
                temp = 0.0
                for k in range(self._cols):
                    temp += self.data[i][k] * mat2.data[k][j]
                res[i][j] = temp
        self._cols = mat2._cols
        self.data = res
        return self
    def show(self):
        print('[')
        for i in range(self._rows):
            print('[', end='')
            for j in range(self._cols):
                print(self.data[i][j], end=',')
            print("]\n")
        print(']')


m1 = MatrixPy(16, 16, 1.0)
m2 = MatrixPy(16, 16, 1.0)
start = time.time()
m1.mul(m2)
end = time.time()
print("time cost:", float(end - start) * 1000.0, "ms")
'''


#加速版本


'''
import ctypes
dll = ctypes.cdll.LoadLibrary
lib = dll('./test.so')
lib.multest(3, 2, 1.0, 2 , 5, 1.0)
'''



