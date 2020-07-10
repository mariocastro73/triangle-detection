#!/usr/bin/env python
## Mario Castro and David Contreras
# December 14th, 2019
# USAGE
# python triangle_detector.py -i input.png -o output.png -c output.csv -r [min,max]

from skimage import io, filters,img_as_ubyte
import argparse
import imutils
import cv2
import numpy as np
import matplotlib.pyplot as plt

# This function uses CV2 to identify the shape of the pattern 
def FindShapes(c):
	shape = "unidentified"
	peri = cv2.arcLength(c, True)
	approx = cv2.approxPolyDP(c, 0.04 * peri, True)

	# if the shape is a triangle, it will have 3 vertices
	if len(approx) == 3:
		shape = "triangle"

	# if the shape has 4 vertices, it is either a square or
	# a rectangle
	elif len(approx) == 4:
		# compute the bounding box of the contour and use the
		# bounding box to compute the aspect ratio
		(x, y, w, h) = cv2.boundingRect(approx)
		ar = w / float(h)

		# a square will have an aspect ratio that is approximately
		# equal to one, otherwise, the shape is a rectangle
		shape = "square" if ar >= 0.95 and ar <= 1.05 else "rectangle"

	# if the shape is a pentagon, it will have 5 vertices
	elif len(approx) == 5:
		shape = "pentagon"

	# otherwise, we assume the shape is a circle
	else:
		shape = "circle"

	# return the name of the shape
	return shape

# From the thresholed image, this method extracts the center of mass of
# structures with an area within values tmin and tmax
def Contours(image,thresh,tmin=0,tmax=100000):
    cnts = cv2.findContours(thresh.copy(), cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
    cnts = imutils.grab_contours(cnts)
    
    counter=0
    for c in cnts: 
        M = cv2.moments(c) 
        area = M["m00"] 
        shape = FindShapes(c) 
        if area < tmax and area > tmin: 
            cX = int((M["m10"] / M["m00"]) ) 
            cY = int((M["m01"] / M["m00"]) ) 
            print(counter,cX,cY,sep=',',file=f)  # output a csv row to be processed
            cv2.circle(image, (cX, cY), 3, (0, 255,0), -1) 
            counter=counter+1
    return counter



# The main code starts here
# construct the argument parse and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True,
	help="path to the input image")
ap.add_argument("-o", "--output", required=True,
	help="path to the output image")
ap.add_argument("-c", "--csv", required=True,
	help="path to the output csv file")
ap.add_argument("-r", "--range", required=False,
	help="threshold interval",default="[60,185]")
args = vars(ap.parse_args())

image = cv2.imread(args["input"]) # Read input image
f = open(args["csv"],'w')                                                                                                        


range = args["range"][1:-1] # Ranges used in the threshold function
aux= range.split(',')
min=int(aux[0])
max=int(aux[1])


gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY) # Convert image to grayscale
blurred = cv2.GaussianBlur(gray,(3,3),0)  # blur image for more accurate results
thresh=cv2.threshold(blurred,min,max,cv2.THRESH_BINARY)[1] # posterize image to find shapes

image2 = image.copy() # keep original image
matched = Contours(image2,thresh,tmin=15,tmax=60) # Find contours with custom function defined above

output = args["output"] # Create output image
cv2.imwrite(output,image2) # Store image 

