#! /usr/bin/env python

from paraview.simple import *

oceanVorticesvtu = XMLUnstructuredGridReader(FileName=['oceanVortices.vtu'])
oceanVorticesvtu.CellArrayStatus = ['adt']

# isolate the continents
threshold2 = Threshold(Input=oceanVorticesvtu)
threshold2.Scalars = ['CELLS', 'adt']
threshold2.ThresholdRange = [-214748.3647, -1.74]

# modify their scalar value the global minimum of the data
calculator1 = Calculator(Input=threshold2)
calculator1.AttributeType = 'Cell Data'
calculator1.ResultArrayName = 'adt'
calculator1.Function = '-1.75'

# isolate the seas
threshold1 = Threshold(Input=oceanVorticesvtu)
threshold1.Scalars = ['CELLS', 'adt']
threshold1.ThresholdRange = [-1.75, 1.8985]

# merge the seas and the new continents
appendDatasets1 = AppendDatasets(Input=[threshold1, calculator1])
cleantoGrid1 = CleantoGrid(Input=appendDatasets1)
tetrahedralize1 = Tetrahedralize(Input=cleantoGrid1)

# pass the data from the cell to the vertices
cellDatatoPointData1 = CellDatatoPointData(Input=tetrahedralize1)
cellDatatoPointData1.CellDataArraytoprocess = ['adt']

SaveData('oceanVortices_ready.vtu', proxy=cellDatatoPointData1)
