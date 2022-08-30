#include "GaudiDDLPlanarDigiProcessor.h"

// Gaudi
#include "GaudiKernel/MsgStream.h" 

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"

#include <TMath.h>

#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h>

#include "AIDA/AIDA.h"

#include "CLHEP/Vector/TwoVector.h"

#include <k4Interface/IUniqueIDGenSvc.h>
#include <DD4hep/BitFieldCoder.h>

// k4FWCore
#include <k4FWCore/DataHandle.h>
#include <k4FWCore/PodioDataSvc.h>

// std
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <climits>
#include <cfloat>

using namespace std ;
using namespace edm4hep;

DECLARE_COMPONENT(GaudiDDPlanarDigiProcessor)

GaudiDDPlanarDigiProcessor::GaudiDDPlanarDigiProcessor(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc), m_eventDataSvc("EventDataSvc", "GaudiDDPlanarDigiProcessor") {
  // Initializing histogram var
  m_hu1D = m_hv1D = m_hT1D = 0;
  m_diffu1D = m_diffv1D = m_diffT1D = 0;
  m_hitE1D = m_hitsAccepted1D = 0;


  // Initialize DataHandlers 
  declareProperty("simtrackhits_r", m_generalSimTrackerHitHandle, "Dummy Hit collection (input)"); // sth_coll
  //declareProperty("simtrackhits", m_simTrackerHitHandle, "Dummy Hit collection (input)");
  declareProperty("VXDTrackerHits", m_TrackerHitHandle, "Dummy Hit collection (output)");
  //declareProperty("simtrackhitsAssociation", m_SimTrackerHitAssociation, "Dummy Hit association collection (output)"); // association

}

//streamlog::out.init(std::cout, "k4MarlinWrapper");
GaudiDDPlanarDigiProcessor::~GaudiDDPlanarDigiProcessor() {
  
};

StatusCode GaudiDDPlanarDigiProcessor::initialize() {

  // initalize global marlin information, maybe betters as a _tool_

  if (GaudiAlgorithm::initialize().isFailure()) {
    return StatusCode::FAILURE;
  }

  m_eventDataSvc.retrieve().ignore();
  m_podioDataSvc = dynamic_cast<PodioDataSvc*>(m_eventDataSvc.get());
  if (m_podioDataSvc == nullptr) {
    return StatusCode::FAILURE;
  }

  // TODO: parse parameters from Parameters Property
  //parseParameters(m_parameters, m_verbosity);
  // TODO: parse parameters, HINT: this is done in k4MarlinWrapper

  // usually a good idea to
  //printParameters();  // HINT: this is done in k4MarlinWrapper

  _nRun = 0 ;
  _nEvt = 0 ;

  // TODO: initialize gsl random generator
  _rng = gsl_rng_alloc(gsl_rng_ranlxs2);

  // TODO: 
  //Global::EVENTSEEDER->registerProcessor(this);
  // GaudiDDLPlanarDigiProcessor.
  m_uniqueIDService = serviceLocator()->service("UniqueIDGenSvc");

  // From Global::EVENTSEEDER->registerProcessor(this);
  // DDLPlanarDigiProcessor.cpp
  // To 
  m_uniqueIDService = serviceLocator()->service("UniqueIDGenSvc");
  m_uniqueIDService->getUniqueID(1, 2, name());

  // DDLPlanarDigiProcessor.h
  SmartIF<IUniqueIDGenSvc> m_uniqueIDService;


  if( m_resU.size() !=  m_resV.size() ) {
    
    std::stringstream ss ;
    ss << name() << "::initialize() - Inconsistent number of resolutions given for U and V coordinate: " 
       << "ResolutionU  :" <<   m_resU.size() << " != ResolutionV : " <<  m_resV.size() ;

    // TODO: how do I want to solve exceptions See MarlinProcessorWrapper for reference
    //throw EVENT::GaudiException( ss.str() ) ; 
  }
  dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();


  // //===========  get the surface map from the SurfaceManager ================

  dd4hep::rec::SurfaceManager& surfMan = *theDetector.extension<dd4hep::rec::SurfaceManager>() ;
  
  //m_subDetName = "Vertex";
  dd4hep::DetElement det = theDetector.detector( m_subDetName ) ;

  _map = surfMan.map( det.name() ) ;

  if( ! _map ) {   
    std::stringstream err  ; err << " Could not find surface map for detector: " 
                                 << m_subDetName << " in SurfaceManager " ;
    //throw Exception( err.str() ) ;
  }


  // streamlog_out( DEBUG3 ) << " DDPlanarDigiProcessor::init(): found " << _map->size() 
  //                         << " surfaces for detector:" <<  m_subDetName << std::endl ;

  // streamlog_out( MESSAGE ) << " *** DDPlanarDigiProcessor::init(): creating histograms" << std::endl ;

  cout << " DDPlanarDigiProcessor::init(): found " << _map->size() << " surfaces for detector:" <<  m_subDetName << std::endl ;
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
// TODO: genereador lo voy a declarar e inicializar en el execute gsl





StatusCode GaudiDDPlanarDigiProcessor::execute() {
  
  /* TODO: Global::EVENTSEEDER and streamlog_out */
  m_uniqueIDService->getUniqueID(1, 2, name());
  // TODO: revisar que estoy haciendo el random f=gsl_rng_set

  cout << "seed set to " << m_uniqueIDService->getUniqueID(1, 2, name()) << std::endl;
  
  // TODO/FIXME: check for failure?
  const auto sth_coll = m_generalSimTrackerHitHandle.get(); // READER, analogous to STHcol

  std::string cellIDEncodingString = m_generalSimTrackerHitHandle.getCollMetadataCellID(sth_coll->getID());


  // EventHeaderCollection* STHcol = 0 ;
  // try{
  //   STHcol =  evt->getCollection( m_inColName ) ;
  // }
  // catch(DataNotAvailableException &e){
  //   streamlog_out(DEBUG4) << "Collection " << m_inColName.c_str() << " is unavailable in event " << _nEvt << std::endl;
  // }
  // Have to wait
    
  unsigned nCreatedHits=0;
  unsigned nDismissedHits=0;
    

  // set the CellIDEncodingString MetaData parameter
  edm4hep::TrackerHitPlaneCollection* trkhitVec = m_TrackerHitHandle.createAndPut();

  std::cout << trkhitVec->getID() << std::endl;
  auto& collmd = m_podioDataSvc->getProvider().getCollectionMetaData(trkhitVec->getID()); // I believe its sth_coll
 
  collmd.setValue("CellIDEncodingString", cellIDEncodingString); // encoding
  
  // Creating the instance of the encoder to decode
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(cellIDEncodingString); 
  
  // Output relation collection TODO: dont know if its Tracker or SimTrackerColl
  edm4hep::TrackerHitCollection* thsthcol = m_OutTrackerHitColl.createAndPut(); // TODO: provisional
  
  // Relation collection TrackerHit, SimTrackerHit
  //create the association between the simhit and the trackerhit
  // TODO: change name of thitNav and see if its correct
  edm4hep::MutableMCRecoTrackerHitPlaneAssociation thitNav = edm4hep::MutableMCRecoTrackerHitPlaneAssociation();
  
  
  int nSimHits = sth_coll->size();

  // TODO
  std::cout << " processing collection " << m_inColName  << " with " <<  nSimHits  << " hits ... " << std::endl ;
  
  for(int i=0; i< nSimHits; ++i) {
    // from an element of a SimTrackerHitCollection to a SimTrackerHit
    //TODO: is it important the dynamic_cast???
    // In this case [] better that [) bc it doesnt check against bounds
    edm4hep::SimTrackerHit simTHit =  (*sth_coll)[i] ;
    
    m_hitE1D->fill(simTHit.getEDep() * (dd4hep::GeV / dd4hep::keV));

    if( simTHit.getEDep() < m_minEnergy ) {
      cout << "Hit with insufficient energy " << simTHit.getEDep() * (dd4hep::GeV / dd4hep::keV) << " keV" << std::endl;
      //streamlog_out( DEBUG ) << "Hit with insufficient energy " << simTHit.EDep() * (dd4hep::GeV / dd4hep::keV) << " keV" << std::endl;
      continue;
    }
    
    const int cellID0 = simTHit.getCellID() ;

    //***********************************************************
    // get the measurement surface for this hit using the CellID
    //***********************************************************
    
    dd4hep::rec::SurfaceMap::const_iterator sI = _map->find( cellID0 ) ;

    if( sI == _map->end() ){    

      std::cout<< " DDPlanarDigiProcessor::processEvent(): no surface found for cellID : " 
                <<   bitFieldCoder.valueString(cellID0) << std::endl;

      
      std::stringstream err ; err << " DDPlanarDigiProcessor::processEvent(): no surface found for cellID : " 
                                  <<   bitFieldCoder.valueString(cellID0)  ;
      throw std::logic_error ( err.str() ) ;
    }

    // bitifeld for the BitFIeldEncoder: cellID0

    const dd4hep::rec::ISurface* surf = sI->second ;

    int layer = bitFieldCoder.get(cellID0, "layer"); // decode

    dd4hep::rec::Vector3D oldPos( simTHit.getPosition()[0], simTHit.getPosition()[1], simTHit.getPosition()[2] );
    
    dd4hep::rec::Vector3D newPos ;

    //************************************************************
    // Check if Hit is inside sensitive 
    //************************************************************
    
    if ( ! surf->insideBounds( dd4hep::mm * oldPos ) ) {
      
      std::cout               << "  hit at " << oldPos 
                              << " " << bitFieldCoder.valueString(cellID0)
                              << " is not on surface " 
                              << *surf  
                              << " distance: " << surf->distance(  dd4hep::mm * oldPos )
                              << std::endl;        

      
      
      
      if( m_forceHitsOntoSurface ) {
        
        dd4hep::rec::Vector2D lv = surf->globalToLocal( dd4hep::mm * oldPos  ) ;
        
        dd4hep::rec::Vector3D oldPosOnSurf = (1./dd4hep::mm) * surf->localToGlobal( lv ) ; 
        
        std::cout << " moved to " << oldPosOnSurf << " distance " << (oldPosOnSurf-oldPos).r()
                                << std::endl;        
          
        oldPos = oldPosOnSurf ;
 
      } else {

          ++nDismissedHits;
        
          continue; 
        }
    }
  

      //***************************************************************
      // Smear time of the hit and apply the time window cut if needed
      //***************************************************************
      
      double hitT = simTHit.getTime();
      
      // Smearing time of the hit
      // now m_resT is a Property with multiple floats
      if (m_resT.size() and m_resT[0] > 0.0) {
        float resT = m_resT.size() > 1 ? m_resT[layer] : m_resT[0];
        double tSmear  = resT > 0.0 ? gsl_ran_gaussian( _rng, resT ) : 0.0;
        m_hT1D->fill( resT > 0.0 ? tSmear / resT : 0.0 );
        m_diffT1D->fill( tSmear );

        hitT += tSmear;
        std::cout << "smeared hit at T: " << simTHit.getTime() << " ns to T: " << hitT << " ns according to resolution: " << resT << " ns" << std::endl;
      }
    
      // Correcting for the propagation time
      if (m_correctTimesForPropagation) {
        double dt = oldPos.r() / ( TMath::C() / 1e6 );
        hitT -= dt;
        std::cout << "corrected hit at R: " << oldPos.r() << " mm by propagation time: " << dt << " ns to T: " << hitT << " ns" << std::endl;
      }
      
      // Skipping the hit if its time is outside the acceptance time window
      if (m_useTimeWindow) {
        float timeWindow_min = m_timeWindow_min.size() > 1 ? m_timeWindow_min[layer] : m_timeWindow_min[0];
        float timeWindow_max = m_timeWindow_max.size() > 1 ? m_timeWindow_max[layer] : m_timeWindow_max[0];
        if ( hitT < timeWindow_min || hitT > timeWindow_max ) {
          std::cout << "hit at T: " << simTHit.getTime() << " smeared to: " << hitT << " is outside the time window: hit dropped"  << std::endl;
          ++nDismissedHits;
          continue;
        }
      }

      //*********************************************************************************
      // Try to smear the hit position but ensure the hit is inside the sensitive region
      //*********************************************************************************
      
      dd4hep::rec::Vector3D u = surf->u() ;
      dd4hep::rec::Vector3D v = surf->v() ;
      

      // get local coordinates on surface
      dd4hep::rec::Vector2D lv = surf->globalToLocal( dd4hep::mm * oldPos  ) ;
      double uL = lv[0] / dd4hep::mm ;
      double vL = lv[1] / dd4hep::mm ;

      bool accept_hit = false ;
      unsigned  tries   =  0 ;              
      static const unsigned MaxTries = 10 ; 
      
      float resU = ( m_resU.size() > 1 ?   m_resU[  layer ]     : m_resU[0]   )  ;
      float resV = ( m_resV.size() > 1 ?   m_resV[  layer ]     : m_resV[0]   )  ; 


      while( tries < MaxTries ) {
        
        if( tries > 0 ) std::cout << "retry smearing for " <<  bitFieldCoder.valueString(cellID0) << " : retries " << tries << std::endl;
        
        double uSmear  = gsl_ran_gaussian( _rng, resU ) ;
        double vSmear  = gsl_ran_gaussian( _rng, resV ) ;

        
        // dd4hep::rec::Vector3D newPosTmp = oldPos +  uSmear * u ;  
        // if( ! m_isStrip )  newPosTmp = newPosTmp +  vSmear * v ;  
        
        dd4hep::rec::Vector3D newPosTmp;
        if (m_isStrip){
            if (m_subDetName == "SET"){
                double xStripPos, yStripPos, zStripPos;
                //Find intersection of the strip with the z=centerOfSensor plane to set it as the center of the SET strip
                dd4hep::rec::Vector3D simHitPosSmeared = (1./dd4hep::mm) * ( surf->localToGlobal( dd4hep::rec::Vector2D( (uL+uSmear)*dd4hep::mm, 0.) ) );
                zStripPos = surf->origin()[2] / dd4hep::mm ;
                double lineParam = (zStripPos - simHitPosSmeared[2])/v[2];
                xStripPos = simHitPosSmeared[0] + lineParam*v[0];
                yStripPos = simHitPosSmeared[1] + lineParam*v[1];
                newPosTmp = dd4hep::rec::Vector3D(xStripPos, yStripPos, zStripPos);    
            }
            else{
                newPosTmp = (1./dd4hep::mm) * ( surf->localToGlobal( dd4hep::rec::Vector2D( (uL+uSmear)*dd4hep::mm, 0. ) ) );    
            }
        }
        else{
            newPosTmp = (1./dd4hep::mm) * ( surf->localToGlobal( dd4hep::rec::Vector2D( (uL+uSmear)*dd4hep::mm, (vL+vSmear)*dd4hep::mm ) ) );
        }

        std::cout               << " hit at    : " << oldPos 
                                << " smeared to: " << newPosTmp
                                << " uL: " << uL 
                                << " vL: " << vL 
                                << " uSmear: " << uSmear
                                << " vSmear: " << vSmear
                                << std::endl ;


        if ( surf->insideBounds( dd4hep::mm * newPosTmp ) ) {    
          
          accept_hit = true ;
          newPos     = newPosTmp ;

          m_hu1D->fill(  uSmear / resU ) ; 
          m_hv1D->fill(  vSmear / resV ) ; 

          m_diffu1D->fill( uSmear );
          m_diffv1D->fill( vSmear );

          break;  

        } else { 
          
          std::cout               << "  hit at " << newPosTmp 
                                  << " " << bitFieldCoder.valueString(cellID0)
                                  << " is not on surface " 
                                  << " distance: " << surf->distance( dd4hep::mm * newPosTmp ) 
                                  << std::endl;        
        }
        
        ++tries;
      }
      
      if( accept_hit == false ) {
        std::cout << "hit could not be smeared within ladder after " << MaxTries << "  tries: hit dropped"  << std::endl;
        ++nDismissedHits;
        continue; 
      } 
      
      //**************************************************************************
      // Store hit variables to TrackerHitPlaneImpl
      //**************************************************************************

      // create the new trackerHit, to be filled with position, etc.
      auto trkHit = trkhitVec->create();

      // m_simTrackerHitHandle.createAndPut() // why am i doing this??

      //TrackerHitPlaneImpl* trkHit = new TrackerHitPlaneImpl ;
      //edm4hep::SimTrackerHit simTHit =  (*sth_coll)[i] ;  JUST FOR REFERENCE, delete
      /* const int cellID1 = simTHit.getCellID1() ;
      trkHit.setCellID0( cellID0 ) ;
      trkHit.setCellID1( cellID1 ) ; */
      trkHit.setCellID(simTHit.getCellID());
      
      trkHit.setPosition( newPos.const_array()  ) ;
      trkHit.setTime( hitT ) ;
      trkHit.setEDep( simTHit.getEDep() ) ;

      float u_direction[2] ;
      u_direction[0] = u.theta();
      u_direction[1] = u.phi();
      
      float v_direction[2] ;
      v_direction[0] = v.theta();
      v_direction[1] = v.phi();
      
      std::cout  << " U[0] = "<< u_direction[0] << " U[1] = "<< u_direction[1] 
                            << " V[0] = "<< v_direction[0] << " V[1] = "<< v_direction[1]
                            << std::endl ;

      trkHit.setU( u_direction ) ;
      trkHit.setV( v_direction ) ;
      
      trkHit.setDu( resU ) ;

      if( m_isStrip ) {

        // store the resolution from the length of the wafer - in case a fitter might want to treat this as 2d hit ....
        double stripRes = (surf->length_along_v() / dd4hep::mm ) / std::sqrt( 12. ) ;
        trkHit.setDv( stripRes ); 

      } else {
        trkHit.setDv( resV ) ;
      }
      // TODO is from Marlin??????????
      if( m_isStrip ){
        //trkHit.setType( UTIL::set_bit( trkHit.getType() ,  UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ) ) ;
      }

      //**************************************************************************
      // Set Relation to SimTrackerHit
      //**************************************************************************    

      // Set relation with LCRelationNavigator
      // TODO REVISAR
      //thitNav.addRelation(trkHit, simTHit);
      thitNav.setRec(trkHit);
      thitNav.setSim(simTHit);
      
      //**************************************************************************
      // Add hit to collection
      //**************************************************************************    
      
      
      
      ++nCreatedHits;
      
      std::cout << "-------------------------------------------------------" << std::endl;
      
    }      
    m_uniqueIDService->getUniqueID(1, 2, m_outColName);
    // Filling the fraction of accepted hits in the event
    float accFraction = nSimHits > 0 ? float(nCreatedHits) / float(nSimHits) * 100.0 : 0.0;
    m_hitsAccepted1D->fill( accFraction );

    // Create rel // Add collection to event
    //**************************************************************************    
    
    /* As long as I know, this is done automatically by initializing the trkhitVec and thsthcol
    evt->addCollection( trkhitVec , m_outColName ) ; 
    evt->addCollection( thsthcol , _outReleateLCCollection());
    
    //**************************************************************************
    // Add collection to event
    //**************************************************************************    
    
    evt->addCollection( trkhitVec , m_outColName ) ;
    evt->addCollection( thsthcol , _outRelColName ) ;
    */
    std::cout << "Created " << nCreatedHits << " hits, " << nDismissedHits << " hits  dismissed\n";
    

    _nEvt ++ ;

  
    return StatusCode::SUCCESS;
  
  }


StatusCode GaudiDDPlanarDigiProcessor::finalize() { 
    
    gsl_rng_free( _rng );
  
    std::cout << " end()  " << name() 
    << " processed " << _nEvt << " events in " << _nRun << " runs "
    << std::endl ;
    return GaudiAlgorithm::finalize(); 
  }
