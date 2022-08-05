#ifndef TESTFWCORE_HELLOWORLDALG
#define TESTFWCORE_HELLOWORLDALG

#pragma once

// GAUDI
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"

#include <string>
#include <vector>
#include <map>

#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"

#include <TH1F.h>

#include "CLHEP/Vector/TwoVector.h"

#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <climits>
#include <cfloat>


// datamodel
namespace edm4hep {
  class MCParticleCollection;
  class SimTrackerHitCollection;
  class SimCaloHit;
}  // namespace edm4hep

/*
namespace marlin {
  class Processor;
  class StringParameters;
}  // namespace marlin
*/ 
//using namespace lcio ;
//using namespace marlin ;

namespace EVENT {
  class SimTrackerHit;
} // MMM NOT SURE WHAT IS THIS




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
  // Gaudi::Property<std::string> theMessage{this, "PerEventPrintMessage", "Hello ", "The message to printed for each Event"};
  
  // ProcessorType: The Type of the DDPlanarDigiProcessor to use
  Gaudi::Property<std::string>                                     m_processorType{this, "ProcessorType", {}};
  //Gaudi::Property<std::map<std::string, std::vector<std::string>>> m_parameters{this, "Parameters", {}};
  
    
  Gaudi::Property<bool> _isStrip{this, "IsStrip", {}};

  Gaudi::Property<std::vector<float>> _resU{this, "ResolutionU", {}};
  Gaudi::Property<std::vector<float>> _resV{this, "ResolutionV", {}};
  Gaudi::Property<std::vector<float>> _resT{this, "ResolutionT", {}};

  Gaudi::Property<std::string> _subDetName{this, "SubDetectorName", {}};
  
  Gaudi::Property<std::string> _inColName{this, "InputCollections", {}};
  Gaudi::Property<std::string> _outColName{this, "OutputCollection", {}};
  //std::string _outRelColName ;
 
  Gaudi::Property<bool> _forceHitsOntoSurface{this, "ForceHitsOntoSurface", {}};
  Gaudi::Property<double> _minEnergy{this, "MinimumEnergyPerHit", {}};
  Gaudi::Property<bool> _correctTimesForPropagation{this, "CorrectTimesForPropagation", {}};

  Gaudi::Property<bool> _useTimeWindow{this, "UseTimeWindow", {}};
  Gaudi::Property<std::vector<float>>  _timeWindow_min{this, "TimeWindowMin", {}};
  Gaudi::Property<std::vector<float>>  _timeWindow_max{this, "TimeWindowMax", {}};
  

  // Histogram variables
  IHistogram1D* m_hu1D;
  IHistogram1D* m_hv1D;
  IHistogram1D* m_hT1D;

  IHistogram1D* m_diffu1D;
  IHistogram1D* m_diffv1D;
  IHistogram1D* m_diffT1D;

  IHistogram1D* m_hitE1D;
  IHistogram1D* m_hitsAccepted1D;

  std::string   m_verbosity = "MESSAGE"; // Used in parseParameters

  std::vector<std::string> split(const std::string& subject, const std::regex& re);
  std::vector<std::string> split(const std::string& subject);
  void parseParameters(
    const std::map<std::string, std::vector<std::string>> parameters, 
    std::string& verbosity) ;


protected:
  int _nRun;
  int _nEvt;
  gsl_rng* _rng ;
  const dd4hep::rec::SurfaceMap* _map ;

  
  
  

  //std::string        m_verbosity = "MESSAGE";

  //std::vector<TH1F*> _h ;

} ;




#endif /* TESTFWCORE_HELLOWORLDALG */
