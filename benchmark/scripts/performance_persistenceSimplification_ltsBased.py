from paraview.simple import *

if len(sys.argv) != 4:
    print("Usage:")
    print("  " + sys.argv[0] + " <input> <persistence threshold %range> <Use All Cores>")
    sys.exit()

inputData = OpenDataFile(sys.argv[1])

vtkDataArray = servermanager.Fetch(inputData).GetPointData().GetArray(0)
valueRange = vtkDataArray.GetRange()
name = vtkDataArray.GetName()
persistenceThreshold = float(sys.argv[2])*(valueRange[1] - valueRange[0])/100

print("Size of input data: " + str(vtkDataArray.GetNumberOfTuples()))
print("Considering field '" + name + "'")
print("Using persistence threshold: " + str(persistenceThreshold))

# create a new 'TTK Blank'
tTKBlank1 = TTKBlank(Input=inputData)
tTKBlank1.DebugLevel = 0

tTKPLTSimplification1 = TTKPLTSimplification(Input=tTKBlank1)
tTKPLTSimplification1.PersistenceThreshold = persistenceThreshold
tTKPLTSimplification1.InputArray = ['POINTS', name]
tTKPLTSimplification1.UseRegionBasedIterations = 1
tTKPLTSimplification1.EscapeInterval = 1000
tTKPLTSimplification1.UseAllCores = int(sys.argv[3])


for i in range(0, 12):
    print("Run #" + str(i))
    tTKBlank1.Outputintegerargument = i
    tTKPLTSimplification1.UpdatePipeline()
