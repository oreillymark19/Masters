import numpy as np
import matplotlib.pyplot as plt
from nipy import load_image
from nilearn import input_data

#Load image and extract data
#func_img = load_image('/home/moreilly/HCP_BIDS/sub-02/func/sub-02_task-rest_bold.nii.gz')
#func_data = func_img.get_data()

#Define ROI class
# class ROI:
# 	#Define seed voxel
# 	def __init__(self, x, y, z):
# 		self.xcoord = x
# 		self.ycoord = y
# 		self.zcoord = z
# 		self.coords = [x,y,z]

# 	#Create ROI sphere around seed voxel
# 	def makeSphere(self, rad, TR, mem, mem_lev):
# 		self.data = func_data[self.xcoord, self.ycoord, self.zcoord,:]
# 		self.masker = input_data.NiftiSpheresMasker(
# 			self.coords, rad, TR, mem, mem_lev)

PCC = [(34,45,34)]
seed_masker = input_data.NiftiSpheresMasker(PCC, radius = 5, t_r = 1, memory='PCCcache', memory_level=1)
seed_time_series = seed_masker.fit_transform('/home/moreilly/HCP_BIDS/sub-02/func/sub-02_task-rest_bold.nii.gz')
plt.plot(seed_time_series)
plt.show()

