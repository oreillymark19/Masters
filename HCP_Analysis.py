import numpy as np
import matplotlib as plt
from nipy import load_image

func_img = load_image(filename)
func_data = func_img.get_data()
ROI_1 = func_data[ , , ,:]
ROI_2 = func_data[ , , ,:]
ROI_3 = func_data[ , , ,:]

ROI_1_Avg = 
