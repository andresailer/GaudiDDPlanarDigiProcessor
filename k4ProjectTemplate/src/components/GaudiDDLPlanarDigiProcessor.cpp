#include "GaudiDDLPlanarDigiProcessor.h"

// Gaudi
#include "GaudiKernel/MsgStream.h" 

// Marlin 
// We are trying to remove the wrapper, but what happens with Marlin?
/* #include "marlin/ProcessorEventSeeder.h"
#include "marlin/AIDAProcessor.h"
#include "marlin/Global.h"
 */
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"

#include <TMath.h>

#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h>

#include "AIDA/AIDA.h"

#include "CLHEP/Vector/TwoVector.h"

// std
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <climits>
#include <cfloat>

//using namespace lcio ;
//using namespace marlin ; // The point is not to use them
using namespace std ;

// datamodel
namespace edm4hep {
  class MCParticleCollection;
  class SimTrackerHitCollection;
  class SimCaloHit;
}  // namespace edm4hep

DECLARE_COMPONENT(GaudiDDPlanarDigiProcessor)

GaudiDDPlanarDigiProcessor::GaudiDDPlanarDigiProcessor(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {
  // Initializing histogram var
  m_hu1D = m_hv1D = m_hT1D = 0;
  m_diffu1D = m_diffv1D = m_diffT1D = 0;
  m_hitE1D = m_hitsAccepted1D = 0;
}

//streamlog::out.init(std::cout, "k4MarlinWrapper");
GaudiDDPlanarDigiProcessor::~GaudiDDPlanarDigiProcessor() {
  
};

/*
enum {
  hu = 0,
  hv,
  hT,
  hitE,
  hitsAccepted,
  diffu,
  diffv,
  diffT,
  hSize 
} ;
*/

StatusCode GaudiDDPlanarDigiProcessor::initialize() {
  
  // initalize global marlin information, maybe betters as a _tool_

    if (GaudiAlgorithm::initialize().isFailure()) {
      return StatusCode::FAILURE;
    }
    
  // TODO: parse parameters
  // HINT: this is done in k4MArlinWrapper


  cout << "Hello World!\n";
  // usually a good idea to
  //printParameters();  // NEED TO SEE WHERE THIS METHOD IS
  _nRun = 0 ;
  _nEvt = 0 ;

  // initialize gsl random generator
  _rng = gsl_rng_alloc(gsl_rng_ranlxs2);

  //_h.resize( hSize ) ;

  // TODO: 
  //Global::EVENTSEEDER->registerProcessor(this);

  if( _resU.size() !=  _resV.size() ) {
    
    std::stringstream ss ;
    ss << name() << "::initialize() - Inconsistent number of resolutions given for U and V coordinate: " 
       << "ResolutionU  :" <<   _resU.size() << " != ResolutionV : " <<  _resV.size() ;

    //throw EVENT::GaudiException( ss.str() ) ; 
    // HAVE TO SEE HOW DO I WANT TO SOLVE EXCEPTIONS. See MarlinProcessorWrapper for reference
  }
  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();


  // //===========  get the surface map from the SurfaceManager ================

  dd4hep::rec::SurfaceManager& surfMan = *theDetector.extension<dd4hep::rec::SurfaceManager>() ;

  // TODO: parse paramiters from Parameters Property
  _subDetName = "Vertex";
  dd4hep::DetElement det = theDetector.detector( _subDetName ) ;

  //_map = surfMan.map( det.name() ) ;

  // if( ! _map ) {   
  //   std::stringstream err  ; err << " Could not find surface map for detector: " 
  //                                << _subDetName << " in SurfaceManager " ;
  //   //throw Exception( err.str() ) ;
  // }

  // streamlog_out( DEBUG3 ) << " DDPlanarDigiProcessor::init(): found " << _map->size() 
  //                         << " surfaces for detector:" <<  _subDetName << std::endl ;

  // streamlog_out( MESSAGE ) << " *** DDPlanarDigiProcessor::init(): creating histograms" << std::endl ;

  //cout << " DDPlanarDigiProcessor::init(): found " << _map->size() << " surfaces for detector:" <<  _subDetName << std::endl ;
  cout << " *** DDPlanarDigiProcessor::init(): creating histograms" << std::endl ;
  
  //AIDAProcessor::histogramFactory(this) ; //->createHistogram1D( "hMCPEnergy", "energy of the MCParticles", 100 ) ;

  // Book 1D histogram with fixed and variable binning
  m_hu1D    = histoSvc()->book( "hu" , "smearing u" , 50, -5. , +5. );
  m_hv1D    = histoSvc()->book( "hv" , "smearing v" , 50, -5. , +5. );
  m_hT1D    = histoSvc()->book( "hT" , "smearing time" , 50, -5. , +5. );

  m_diffu1D    = histoSvc()->book( "diffu" , "diff u" , 1000, -5. , +5. );
  m_diffv1D    = histoSvc()->book( "diffv" , "diff v" , 1000, -5. , +5. );
  m_diffT1D    = histoSvc()->book( "diffT" , "diff time" , 1000, -5. , +5. );

  m_hitE1D    = histoSvc()->book( "hitE" , "hitEnergy in keV" , 1000, 0 , 200 );
  m_hitsAccepted1D    = histoSvc()->book( "hitsAccepted" , "Fraction of accepted hits [%]" , 201, 0 , 100.5 );

  /*
  _h[ hu ] = new TH1F( "hu" , "smearing u" , 50, -5. , +5. );
  _h[ hv ] = new TH1F( "hv" , "smearing v" , 50, -5. , +5. );
  _h[ hT ] = new TH1F( "hT" , "smearing time" , 50, -5. , +5. );

  _rng = gsl_rng_alloc(gsl_rng_ranlxs2);
f time" , 1000, -5. , +5. );

  _h[ hitE ] = new TH1F( "hitE" , "hitEnergy in keV" , 1000, 0 , 200 );
  _h[ hitsAccepted ] = new TH1F( "hitsAccepted" , "Fraction of accepted hits [%]" , 201, 0 , 100.5 );
  */
  if ( 0 == m_hu1D || 0 == m_hv1D || 0 == m_hT1D || 0 == m_diffu1D || 0 == m_diffv1D || 0 == m_diffT1D || 0 == m_hitE1D ||
       0 == m_hitsAccepted1D ) {
    error() << "----- Cannot book or register histograms -----" << endmsg;
    return StatusCode::FAILURE;
  }
  info() << "Finished booking Histograms" << endmsg;
  return StatusCode::SUCCESS;
}



StatusCode GaudiDDPlanarDigiProcessor::execute() {
  
  return StatusCode::SUCCESS;
  
}

StatusCode GaudiDDPlanarDigiProcessor::finalize() { return GaudiAlgorithm::finalize(); }

