import matplotlib.pyplot as plt
from PIL import Image
import matplotlib.patches as mpatches
import numpy

#im = Image.open('/home/deturres/source/CRF_implementations/fakedata/image_0001.jpg')

im = Image.open('image.jpg')

# Flip the .tif file so it plots upright
im1 = im.transpose(Image.FLIP_TOP_BOTTOM)
width,height = im.size
# Plot the image
plt.imshow(im1)
ax = plt.gca()
ax.set_xticks(numpy.arange(0,width,32))
ax.set_yticks(numpy.arange(0,height,32))
# create a grid
ax.grid(True, color='r', linestyle='-', linewidth=2)
# put the grid below other plot elements
ax.set_axisbelow(True)



#plt.savefig('/home/deturres/source/CRF_implementations/fakedata/image_0001_grid.png',bbox_inches='tight')
plt.show()

