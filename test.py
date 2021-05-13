from wdat_plot import *

# creata Wdat class object by passing prefix to .wtxt file
p = Wdat("./data/example")
# print all elements of p
print(p)


# EXAMPLE
# ./data/example_delta.wdat contains 2D data
# plot delta (complex variable)
imre = p.read("delta")
# imre is 3D numpy.array:
#       imre[0] - 2D array of real part of delta field
#       imre[1] - 2D array of imaginary part of delta field
# convert Re Im format to Abs Arg
delta = p.ReIm2AbsArg(imre)
# delta is 3D numpy.array:
#       delta[0] - 2D array of absolute value of delta field
#       delta[1] - 2D array of argument of delta field

# delta can be easily altered.
# Let's convert delta argument to pi units starting from 0:
delta[1] /= 3.14159
delta[1] += 1

# plot delta argument and save it to the ./data/test.png file
p.Draw2D(delta[1], "./data/test.png")
