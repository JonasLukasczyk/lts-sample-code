#! /usr/bin/env python

from paraview.simple import *

box1 = Box()

resampleToImage1 = ResampleToImage(Input=box1)
resampleToImage1.SamplingDimensions = [256, 256, 256]

randomAttributes1 = RandomAttributes(Input=resampleToImage1)
randomAttributes1.DataType = 'Double'
randomAttributes1.ComponentRange = [0.0, 1.0]
randomAttributes1.AttributesConstantPerBlock = 0
randomAttributes1.GeneratePointScalars = 1
randomAttributes1.GeneratePointVectors = 0
randomAttributes1.GenerateCellScalars = 0
randomAttributes1.GenerateCellVectors = 0

tTKPointDataSelector1 = TTKPointDataSelector(Input=randomAttributes1)
tTKPointDataSelector1.RangeId = [0, 4]
tTKPointDataSelector1.ScalarFields = ['RandomPointScalars']
tTKPointDataSelector1.DebugLevel = 0

SaveData('random.vti', proxy=tTKPointDataSelector1)
