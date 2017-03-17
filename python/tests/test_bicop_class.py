import pyvinecopulib as pvcl

# TODO: test no default constructor

bicop = pvcl.bicop(pvcl.BicopFamily.gaussian, 90, [1])
assert bicop.family == pvcl.BicopFamily.gaussian
assert len(bicop.parameters) == 1 
assert bicop.parameters[0] == 1
assert bicop.rotation == 90
