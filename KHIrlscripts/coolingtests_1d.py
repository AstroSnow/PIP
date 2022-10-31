#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 11:47:00 2022

@author: ben
"""
import numpy as np
import matplotlib.pyplot as plt
import pipreadmods


filename = "/home/ben/Documents/KHI/PIPrl_1d/Data/"
#filename = "/home/ben/Documents/KHI/PIPrl/KHIrldata/CD_800_1.2/"

ds=pipreadmods.pipread(filename,tstep=-1)

T=5.0/3.0*ds['pr_p']/ds['ro_p']*1.0e6

