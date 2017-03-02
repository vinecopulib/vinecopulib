import pyvinecopulib as pvcl

# TODO: test no default constructor

bicop = pvcl.bicop(1, [1], 90)
assert bicop.family == 1
assert len(bicop.parameters) == 1 
assert bicop.parameters[0] == 1
assert bicop.rotation == 90
