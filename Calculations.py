#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 12:23:31 2019

@author: armaanthapar
"""
# Calculations
import numpy as np



# Finds distance between molecules using COM
def molecular_dist(M1, M2):
	# US the COMs for the center
	C1 = M1.center  
	C2 = M2.center
	D = 0
	for i in range(3):  # 3 for x y z 
		D += (C1[i]-C2[i])**2
	D = np.sqrt(D)
	return D

# Takes any two position vectors and finds distance between them
def distance(v1,v2):
	D = 0
	for i in range(len(v1)):
		D += (v1[i] - v2[i])**2
	D = np.sqrt(D)
	return D
		