#ifndef TESTFWCORE_GaudiDDLPlanarDigiProcessor
#define TESTFWCORE_GaudiDDLPlanarDigiProcessor


// GAUDI
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"

#include <string>
#include <vector>
#include <map>

#include <gsl/gsl_rng.h>

#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include <TH1F.h>

#include "CLHEP/Vector/TwoVector.h"

#include <DD4hep/BitFieldCoder.h>

// k4FWCore
#include <k4FWCore/DataHandle.h>
#include <k4FWCore/PodioDataSvc.h>
#include <k4Interface/IUniqueIDGenSvc.h>

// edm4hep
#include "k4FWCore/DataHandle.h"

// datamodel
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/MutableMCRecoTrackerHitPlaneAssociation.h"


#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <climits>
#include <cfloat>


// datamodel
namespace edm4hep {
  class TrackerHitCollection; 
  class TrackerHitPlaneCollection; 
  //class MCParticleCollection;
  class SimTrackerHit;
  class MutableMCRecoTrackerHitPlaneAssociation;
  class SimTrackerHitCollection;
 
  class SimCaloHit;
}  // namespace edm4hep


class GaudiDDPlanarDigiProcessor : public GaudiAlgorithm {
public:
  explicit GaudiDDPlanarDigiProcessor(const std::string&, ISvcLocator*);
  virtual ~GaudiDDPlanarDigiProcessor();
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize() override final;
  /**  Execute.
   *   @return status code
   */
  virtual StatusCode execute() final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

private:
  // member variable
  
  // ProcessorType: The Type of the DDPlanarDigiProcessor to use
  Gaudi::Property<std::string> m_processorType{this, "ProcessorType", {}};

//GaudiDDPlannarDigiProcessor.h
  Gaudi::Property<bool> m_isStrip{this, "IsStrip", {}};
  Gaudi::Property<std::vector<float>> m_resU{this, "ResolutionU", {}};


  Gaudi::Property<std::vector<float>> m_resV{this, "ResolutionV", {}};
  Gaudi::Property<std::vector<float>> m_resT{this, "ResolutionT", {}};

  Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", {}};
  
  Gaudi::Property<std::string> m_inColName{this, "InputCollections", {}};
  Gaudi::Property<std::string> m_outColName{this, "OutputCollection", {}};
  //std::string _outRelColName ;
 
  Gaudi::Property<bool> m_forceHitsOntoSurface{this, "ForceHitsOntoSurface", {}};
  Gaudi::Property<double> m_minEnergy{this, "MinimumEnergyPerHit", {}};
  Gaudi::Property<bool> m_correctTimesForPropagation{this, "CorrectTimesForPropagation", {}};

  Gaudi::Property<bool> m_useTimeWindow{this, "UseTimeWindow", {}};
  Gaudi::Property<std::vector<float>>  m_timeWindow_min{this, "TimeWindowMin", {}};
  Gaudi::Property<std::vector<float>>  m_timeWindow_max{this, "TimeWindowMax", {}};
  // TODO: set default values
  // TODO: allow user to give the collection name of the algorithm 25th Aug
  SmartIF<IUniqueIDGenSvc> m_uniqueIDService;
  
  // PODIO data service
  ServiceHandle<IDataProviderSvc> m_eventDataSvc;

  // GaudiDDLPlanarDigiProcessor.h
  PodioDataSvc* m_podioDataSvc;

  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder;

  // Histogram variables
  IHistogram1D* m_hu1D;
  IHistogram1D* m_hv1D;
  IHistogram1D* m_hT1D;

  IHistogram1D* m_diffu1D;
  IHistogram1D* m_diffv1D;
  IHistogram1D* m_diffT1D;

  IHistogram1D* m_hitE1D;
  IHistogram1D* m_hitsAccepted1D;



  //DataHandle<edm4hep::EventHeaderCollection> m_eventHeaderHandle{}
  // Handle for the SimTrackerHits to be written
  
DataHandle<edm4hep::SimTrackerHitCollection> m_generalSimTrackerHitHandle{"VertexBarrelCollection", Gaudi::DataHandle::Reader, this}; // sth_coll
  //DataHandle<edm4hep::SimTrackerHitCollection> m_simTrackerHitHandle{"SimTrackerHits", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHitPlaneCollection> m_TrackerHitHandle{"VXDTrackerHits", Gaudi::DataHandle::Writer, this};
  //DataHandle<edm4hep::TrackerHitCollection> m_SimTrackerHitAssociation{"SimTrackerHits", Gaudi::DataHandle::Writer, this}; // Output association 
  DataHandle<edm4hep::TrackerHitCollection> m_OutTrackerHitColl{"OutputTrackerHitColl", Gaudi::DataHandle::Writer, this}; // Output trackerhitcoll

protected:
  int _nRun;
  int _nEvt;
  gsl_rng* _rng ;
  const dd4hep::rec::SurfaceMap* _map ;

  

  //std::vector<TH1F*> _h ;

} ;




#endif /* TESTFWCORE_HELLOWORLDALG */
