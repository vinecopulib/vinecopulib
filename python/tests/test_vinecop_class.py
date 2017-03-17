import pyvinecopulib as pvcl
import numpy as np

# no-arg ctor
vinecop = pvcl.vinecop()

# 2-arg ctor
vinecop = pvcl.vinecop([
  [pvcl.bicop(3, [1], 90),    pvcl.bicop(2, [1,2], 180),  pvcl.bicop(1, [1], 270) ],
  [pvcl.bicop(2, [1,2], 180), pvcl.bicop(3, [1], 90)                              ],
  [pvcl.bicop(1, [1], 270)                                                        ]
], np.array([
  [1,1,1,1],
  [1,1,1,0],
  [1,1,0,0],
  [1,0,0,0]
], dtype=np.int32))

# accessors
print(vinecop.rotation(1,1))
print(vinecop.rotations)
print(vinecop.parameters(0,1))
print(vinecop.family(0,0))
print(vinecop.families)
print(vinecop.matrix)
#print(vinecop.pair_copula(0,0)) - TODO

# 2-arg simulate
#vinecop.simulate(3, np.array([
#  [ 1, 2, 3], 
#  [ 4, 5, 6], 
#  [ 7, 8, 9], 
#  [10,11,12]
#], dtype=np.double));
