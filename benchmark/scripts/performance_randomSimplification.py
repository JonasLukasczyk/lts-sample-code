from paraview.simple import *

if len(sys.argv) != 5:
    print("Usage:")
    print("  " + sys.argv[0] + " <input> <%pairs> <Use All Cores> <Use LPS>")
    sys.exit()

inputData = OpenDataFile(sys.argv[1])

vtkDataArray = servermanager.Fetch(inputData).GetPointData().GetArray(0)
print("Size of input data: " + str(vtkDataArray.GetNumberOfTuples()))

# create a new 'TTK PersistenceDiagram'
tTKPersistenceDiagram1 = TTKPersistenceDiagram(Input=inputData)
tTKPersistenceDiagram1.EmbedinDomain = 1
tTKPersistenceDiagram1.DebugLevel = 0

diagram = servermanager.Fetch(tTKPersistenceDiagram1)
print("Diagram has " + str(diagram.GetNumberOfCells()) + " pairs overall")

# filter minima pairs
minimaPairs = Threshold(Input=tTKPersistenceDiagram1)
minimaPairs.Scalars = ['CELLS', 'PairType']
minimaPairs.ThresholdRange = [-0.1, 0.1]
minimaDiagram = servermanager.Fetch(minimaPairs)
print("Diagram has " + str(minimaDiagram.GetNumberOfCells()) + " minimum pairs")

randomAttributes1 = RandomAttributes(Input=minimaPairs)
randomAttributes1.DataType = 'Int'
randomAttributes1.ComponentRange = [0.0, 9999999.0]
randomAttributes1.GenerateCellScalars = 1
randomAttributes1.GenerateCellVectors = 0

tTKIdentifyByScalarField1 = TTKIdentifyByScalarField(Input=randomAttributes1)
tTKIdentifyByScalarField1.ScalarField = 'RandomCellScalars'
#tTKIdentifyByScalarField1.StartByOne = 1

minPairThreshold = float(minimaDiagram.GetNumberOfCells())*((100.0 - float(sys.argv[2]))/100.0)

threshold1 = Threshold(Input=tTKIdentifyByScalarField1)
threshold1.Scalars = ['CELLS', 'CellScalarFieldName']
threshold1.ThresholdRange = [0, minPairThreshold]

selectedMinPairs = servermanager.Fetch(threshold1)

print("  kept " + str(selectedMinPairs.GetNumberOfCells()) + " minima pairs")

# filter maxima pairs
maximaPairs = Threshold(Input=tTKPersistenceDiagram1)
maximaPairs.Scalars = ['CELLS', 'PairType']
maximaPairs.ThresholdRange = [1.9, 2.1]
maximaDiagram = servermanager.Fetch(maximaPairs)
print("Diagram has " + str(maximaDiagram.GetNumberOfCells()) + " maximum pairs")

randomAttributes2 = RandomAttributes(Input=maximaPairs)
randomAttributes2.DataType = 'Int'
randomAttributes2.ComponentRange = [0.0, 9999999.0]
randomAttributes2.GenerateCellScalars = 1
randomAttributes2.GenerateCellVectors = 0

tTKIdentifyByScalarField2 = TTKIdentifyByScalarField(Input=randomAttributes2)
tTKIdentifyByScalarField2.ScalarField = 'RandomCellScalars'
#tTKIdentifyByScalarField2.StartByOne = 1

maxPairThreshold = float(maximaDiagram.GetNumberOfCells())*((100.0 - float(sys.argv[2]))/100.0)

threshold2 = Threshold(Input=tTKIdentifyByScalarField2)
threshold2.Scalars = ['CELLS', 'CellScalarFieldName']
threshold2.ThresholdRange = [0, maxPairThreshold]

selectedMaxPairs = servermanager.Fetch(threshold2)

print("  kept " + str(selectedMaxPairs.GetNumberOfCells()) + " maxima pairs")

appendDataSets1 = AppendDatasets(Input=[threshold1, threshold2])
print("Keeping " + str(servermanager.Fetch(appendDataSets1).GetNumberOfCells())
+ " pairs overall")

tTKBlank1 = TTKBlank(Input=appendDataSets1)
tTKBlank1.DebugLevel = 0

# create a new 'TTK TopologicalSimplification'
tTKTopologicalSimplification1 = TTKTopologicalSimplification(Domain=inputData,
    Constraints=tTKBlank1)
tTKTopologicalSimplification1.UseTPTS = int(sys.argv[4])
tTKTopologicalSimplification1.UseAllCores = int(sys.argv[3])
tTKTopologicalSimplification1.UseRegionBasedIterations = 1
#tTKTopologicalSimplification1.DebugLevel = 4

for i in range(0, 12):
    print("Run #" + str(i))
    tTKBlank1.Outputintegerargument = i
    tTKTopologicalSimplification1.UpdatePipeline()
    #if i == 0:
    #    SaveData("simplifiedRandom_" + str(sys.argv[2]) + ".vti", proxy=tTKTopologicalSimplification1)
