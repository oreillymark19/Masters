import numpy as np
import matplotlib as plt
from nipy import load_image
from nilearn import input_data

func_img = load_image(filename)
func_data = func_img.get_data()

class ROI:

	def __init__(self, x, y, z):
		self.xcoord = x
		self.ycoord = y
		self.zcoord = z
		self.coords = [x,y,z]

	def makeSphere(self, rad, TR, mem, mem_lev):
		self.data = func_img[self.xcoord, self.ycoord, self.zcoord,:]
		self.masker = input_data.NiftiSpheresMasker(
			self.coords, rad, TR, mem, mem_lev)

PCC = ROI(1,5,8)
PCC_ROI = PCC.makeSphere(5, 1000, 'PCCseedcache', 1)

