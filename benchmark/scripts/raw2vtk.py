#! /usr/bin/env python

from paraview.simple import *

if len(sys.argv) != 7:
    print("Usage:")
    print("  " + sys.argv[0] + " <input> <output> <x> <y> <z> <vtk type>")
    sys.exit()

rawGrid = ImageReader(FileNames=[sys.argv[1]])
rawGrid.DataScalarType = sys.argv[6]
rawGrid.DataByteOrder = 'LittleEndian'
rawGrid.DataExtent = [0, int(sys.argv[3]), 0, int(sys.argv[4]), 0, int(sys.argv[5])]

SaveData(sys.argv[2], proxy=rawGrid)
