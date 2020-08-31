from paraview.simple import *

if len(sys.argv) != 5:
    print("Usage:")
    print("  " + sys.argv[0] + " <input> <persistence threshold %range> <Use All Cores> <Use LPS>")
    sys.exit()

inputData = OpenDataFile(sys.argv[1])

vtkDataArray = servermanager.Fetch(inputData).GetPointData().GetArray(0)
valueRange = vtkDataArray.GetRange()
persistenceThreshold = float(sys.argv[2])*(valueRange[1] - valueRange[0])/100

print("Size of input data: " + str(vtkDataArray.GetNumberOfTuples()))

print("Using persistence threshold: " + str(persistenceThreshold))

# create a new 'TTK PersistenceDiagram'
tTKPersistenceDiagram1 = TTKPersistenceDiagram(Input=inputData)
tTKPersistenceDiagram1.EmbedinDomain = 1
tTKPersistenceDiagram1.DebugLevel = 0
tTKPersistenceDiagram1.UseAllCores = 0

# create a new 'TTK Blank'
tTKBlank1 = TTKBlank(Input=tTKPersistenceDiagram1)
tTKBlank1.DebugLevel = 0

# create a new 'Threshold'
threshold1 = Threshold(Input=tTKBlank1)
threshold1.Scalars = ['CELLS', 'Persistence']
threshold1.ThresholdRange = [persistenceThreshold, 9999999.0]

# create a new 'TTK TopologicalSimplification'
tTKTopologicalSimplification1 = TTKTopologicalSimplification(Domain=inputData,
    Constraints=threshold1)
tTKTopologicalSimplification1.UseTPTS = int(sys.argv[4])
tTKTopologicalSimplification1.UseAllCores = int(sys.argv[3])
tTKTopologicalSimplification1.UseRegionBasedIterations = 1
#tTKTopologicalSimplification1.DebugLevel = 4

for i in range(0, 12):
    print("Run #" + str(i))
    tTKBlank1.Outputintegerargument = i
    tTKTopologicalSimplification1.UpdatePipeline()
