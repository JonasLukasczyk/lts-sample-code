/// \ingroup vtk
/// \class ttkMorseSmaleComplex
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date February 2017.
///
/// \brief TTK VTK-filter that wraps the morseSmaleComplex processing package.
///
/// TTK module for the computation of Morse-Smale complexes.
/// Morse-Smale complexes are useful topological abstractions of scalar
/// fields for data segmentation, feature extraction, etc.
///
/// \b Related \b publication \n
/// "Parallel Computation of 3D Morse-Smale Complexes" \n
/// Nithin Shivashankar, Vijay Natarajan \n
/// Proc. of EuroVis 2012. \n
/// Computer Graphics Forum, 2012.
///
/// \param Input Input scalar field, defined as a point data scalar field
/// attached to a geometry, either 2D or 3D, either regular grid or
/// triangulation (vtkDataSet)
/// \param Output0 Output critical points (vtkUnstructuredGrid)
/// \param Output1 Output 1-separatrices (vtkUnstructuredGrid)
/// \param Output2 Output 2-separatrices (vtkUnstructuredGrid)
/// \param Output3 Output data segmentation (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::MorseSmaleComplex
///
#ifndef _TTK_MORSESMALECOMPLEX_H
#define _TTK_MORSESMALECOMPLEX_H

// VTK includes -- to adapt
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

// VTK Module
#include <ttkMorseSmaleComplexModule.h>

// ttk code includes
#include <MorseSmaleComplex.h>
#include <ttkTriangulationAlgorithm.h>

class TTKMORSESMALECOMPLEX_EXPORT ttkMorseSmaleComplex
  : public vtkDataSetAlgorithm,
    public ttk::Wrapper {

public:
  static ttkMorseSmaleComplex *New();

  vtkTypeMacro(ttkMorseSmaleComplex, vtkDataSetAlgorithm);

  // default ttk setters
  vtkSetMacro(debugLevel_, int);

  void SetThreadNumber(int threadNumber) {
    ThreadNumber = threadNumber;
    SetThreads();
  }

  void SetUseAllCores(bool onOff) {
    UseAllCores = onOff;
    SetThreads();
  }
  // end of default ttk setters

  vtkSetMacro(ScalarField, std::string);
  vtkGetMacro(ScalarField, std::string);

  vtkSetMacro(ScalarFieldId, int);
  vtkGetMacro(ScalarFieldId, int);

  vtkSetMacro(OffsetFieldId, int);
  vtkGetMacro(OffsetFieldId, int);

  vtkSetMacro(ForceInputOffsetScalarField, int);
  vtkGetMacro(ForceInputOffsetScalarField, int);

  vtkSetMacro(InputOffsetScalarFieldName, std::string);
  vtkGetMacro(InputOffsetScalarFieldName, std::string);

  vtkSetMacro(PeriodicBoundaryConditions, int);
  vtkGetMacro(PeriodicBoundaryConditions, int);

  vtkSetMacro(IterationThreshold, int);
  vtkGetMacro(IterationThreshold, int);

  vtkSetMacro(ReverseSaddleMaximumConnection, int);
  vtkGetMacro(ReverseSaddleMaximumConnection, int);

  vtkSetMacro(ReverseSaddleSaddleConnection, int);
  vtkGetMacro(ReverseSaddleSaddleConnection, int);

  vtkSetMacro(ComputeCriticalPoints, int);
  vtkGetMacro(ComputeCriticalPoints, int);

  vtkSetMacro(ComputeAscendingSeparatrices1, int);
  vtkGetMacro(ComputeAscendingSeparatrices1, int);

  vtkSetMacro(ComputeDescendingSeparatrices1, int);
  vtkGetMacro(ComputeDescendingSeparatrices1, int);

  vtkSetMacro(ComputeSaddleConnectors, int);
  vtkGetMacro(ComputeSaddleConnectors, int);

  vtkSetMacro(ComputeAscendingSeparatrices2, int);
  vtkGetMacro(ComputeAscendingSeparatrices2, int);

  vtkSetMacro(ComputeDescendingSeparatrices2, int);
  vtkGetMacro(ComputeDescendingSeparatrices2, int);

  vtkSetMacro(ComputeAscendingSegmentation, int);
  vtkGetMacro(ComputeAscendingSegmentation, int);

  vtkSetMacro(ComputeDescendingSegmentation, int);
  vtkGetMacro(ComputeDescendingSegmentation, int);

  vtkSetMacro(ComputeFinalSegmentation, int);
  vtkGetMacro(ComputeFinalSegmentation, int);

  vtkSetMacro(ReturnSaddleConnectors, int);
  vtkGetMacro(ReturnSaddleConnectors, int);

  vtkSetMacro(SaddleConnectorsPersistenceThreshold, double);
  vtkGetMacro(SaddleConnectorsPersistenceThreshold, double);

  vtkSetMacro(PrioritizeSpeedOverMemory, int);
  vtkGetMacro(PrioritizeSpeedOverMemory, int);

  int setupTriangulation(vtkDataSet *input);
  vtkDataArray *getScalars(vtkDataSet *input);
  vtkDataArray *getOffsets(vtkDataSet *input);

protected:
  template <typename VTK_TT>
  int dispatch(
    vtkDataArray *inputScalars,
    vtkDataArray *inputOffsets,
    vtkUnstructuredGrid *outputCriticalPoints,
    vtkUnstructuredGrid *outputSeparatrices1,
    vtkUnstructuredGrid *outputSeparatrices2,
    ttk::SimplexId criticalPoints_numberOfPoints,
    std::vector<float> &criticalPoints_points,
    std::vector<char> &criticalPoints_points_cellDimensions,
    std::vector<ttk::SimplexId> &criticalPoints_points_cellIds,
    std::vector<char> &criticalPoints_points_isOnBoundary,
    std::vector<ttk::SimplexId> &criticalPoints_points_PLVertexIdentifiers,
    std::vector<ttk::SimplexId> &criticalPoints_points_manifoldSize,
    ttk::SimplexId separatrices1_numberOfPoints,
    std::vector<float> &separatrices1_points,
    std::vector<char> &separatrices1_points_smoothingMask,
    std::vector<char> &separatrices1_points_cellDimensions,
    std::vector<ttk::SimplexId> separatrices1_points_cellIds,
    ttk::SimplexId separatrices1_numberOfCells,
    std::vector<ttk::SimplexId> &separatrices1_cells,
    std::vector<ttk::SimplexId> &separatrices1_cells_sourceIds,
    std::vector<ttk::SimplexId> &separatrices1_cells_destinationIds,
    std::vector<ttk::SimplexId> &separatrices1_cells_separatrixIds,
    std::vector<char> &separatrices1_cells_separatrixTypes,
    std::vector<char> &separatrices1_cells_isOnBoundary,
    ttk::SimplexId separatrices2_numberOfPoints,
    std::vector<float> &separatrices2_points,
    ttk::SimplexId separatrices2_numberOfCells,
    std::vector<ttk::SimplexId> &separatrices2_cells,
    std::vector<ttk::SimplexId> &separatrices2_cells_sourceIds,
    std::vector<ttk::SimplexId> &separatrices2_cells_separatrixIds,
    std::vector<char> &separatrices2_cells_separatrixTypes,
    std::vector<char> &separatrices2_cells_isOnBoundary);

  ttkMorseSmaleComplex();
  ~ttkMorseSmaleComplex() override;

  TTK_SETUP();

  virtual int FillInputPortInformation(int port, vtkInformation *info) override;
  virtual int FillOutputPortInformation(int port,
                                        vtkInformation *info) override;

private:
  std::string ScalarField;
  std::string InputOffsetScalarFieldName;
  bool ForceInputOffsetScalarField;
  bool PeriodicBoundaryConditions;
  int IterationThreshold;
  bool ReverseSaddleMaximumConnection;
  bool ReverseSaddleSaddleConnection;
  bool ComputeCriticalPoints;
  bool ComputeAscendingSeparatrices1;
  bool ComputeDescendingSeparatrices1;
  bool ComputeSaddleConnectors;
  bool ComputeAscendingSeparatrices2;
  bool ComputeDescendingSeparatrices2;
  bool ComputeAscendingSegmentation;
  bool ComputeDescendingSegmentation;
  bool ComputeFinalSegmentation;
  int ScalarFieldId;
  int OffsetFieldId;
  int ReturnSaddleConnectors;
  double SaddleConnectorsPersistenceThreshold;
  bool PrioritizeSpeedOverMemory;

  ttk::MorseSmaleComplex morseSmaleComplex_;
  ttk::Triangulation *triangulation_;
  vtkDataArray *defaultOffsets_;
  bool hasUpdatedMesh_;
};

#endif // _TTK_MORSESMALECOMPLEX_H
