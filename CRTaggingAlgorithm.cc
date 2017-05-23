/**
 *  @file   larpandoracontent/LArCRTagging/CRTaggingAlgorithm.cc
 * 
 *  @brief  Implementation of the cosmic ray tagging algorithm class.
 * 
 *  $Log: $
 */

#include "larpandoracontent/LArCRTagging/CRTaggingAlgorithm.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "Pandora/PdgTable.h"
#include <iomanip>

using namespace pandora;

namespace lar_content
{

CRTaggingAlgorithm::CRTaggingAlgorithm() :
  m_eventNumber(0),
  m_maxPhotonPropagation(2.5f),
  m_matchingMinPrimaryHits(15),
  m_matchingMinHitsForGoodView(5),
  m_matchingMinPrimaryGoodViews(2),
  m_angularUncertainty( 10.0 ),
  m_positionalUncertainty( 3.0 ),
  m_maxAssociationDist( 3 * 18.0 ) // 3 photon conversion lengths in Ar
{
}

CRTaggingAlgorithm::~CRTaggingAlgorithm() {
  PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "PFOs"       , "EventData.root", "UPDATE"));
  PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "Hits"       , "EventData.root", "UPDATE"));
  PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "Targets"    , "EventData.root", "UPDATE"));
  PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "CosmicMuons", "EventData.root", "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

/* TESTING ================================================= */
/*
void PrintDaughterInfo( const MCParticle * const pMCParticle, int depth, std::map< const MCParticle * const, int > & mcParticleToNHitsMap ) {
 
  int             pdg      = pMCParticle->GetParticleId();
  
  std::cout << std::string( 3*( depth - 1), ' ') << "|- " << std::left << std::setw(12) << pdg;

  int nHits = 0;
  if ( mcParticleToNHitsMap.find( pMCParticle ) != mcParticleToNHitsMap.end() ){
    nHits = mcParticleToNHitsMap[ pMCParticle ];
    std::cout << std::string( 3 * (10 - depth), ' ') << nHits;
  }

  std::cout << std::endl;

  for ( const MCParticle * const pMCDaughter : pMCParticle->GetDaughterList() ) {
    PrintDaughterInfo( pMCDaughter, depth + 1, mcParticleToNHitsMap );
  }
  
}
*/
/* END TESTING ============================================= */

StatusCode CRTaggingAlgorithm::Run()
{

  // ===================================================
  //  the input collections
  // ===================================================

  //  the list of PFOs
  const PfoList *pPfoList(nullptr);
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

  if (!pPfoList || pPfoList->empty()) {
    if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo()){
      std::cout << "CRTaggingAlgorithm: unable to find pfo list " << m_inputPfoListName << std::endl;
    }
    return STATUS_CODE_FAILURE;
  }

  //  the list of all the CaloHits (2D)
  const CaloHitList *pCaloHitList = nullptr;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList)); 

  //  the list of MCParticles
  const MCParticleList *pMCParticleList = nullptr;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList)); 






  /* TESTING ================================================= */
  /*
  // Get number of hits associated with each MCParticle
  std::map< const MCParticle * const, int > mcParticleToNHitsMap;
  std::map< const MCParticle * const, int > mcPrimaryToNHitsMap;

  std::map< int,   std::map< int, int > >   hitProductionMap;
  //        ^ PDG of primary ^    ^
  //                         |    | Number of hits
  //                         | PDG of particle producing hits

  for ( const CaloHit * const pCaloHit : *pCaloHitList ) {
    const MCParticle * mainMCParticle = this->GetMainMCParticle(pCaloHit);
    if ( ! mainMCParticle ) continue; // Ghost hit
    if ( mcParticleToNHitsMap.find( mainMCParticle ) == mcParticleToNHitsMap.end() ) {
      mcParticleToNHitsMap.insert( std::make_pair( mainMCParticle, 0 ) );
    }
  
    const MCParticle * primaryMCParticle = LArMCParticleHelper::GetPrimaryMCParticle(mainMCParticle);
    if ( mcPrimaryToNHitsMap.find( primaryMCParticle ) == mcPrimaryToNHitsMap.end() ) {
      mcPrimaryToNHitsMap.insert( std::make_pair( primaryMCParticle, 0 ) );
    }
    
    int thisPdg = mainMCParticle->GetParticleId();
    int primaryPdg = primaryMCParticle->GetParticleId();
    if ( hitProductionMap.find( primaryPdg ) == hitProductionMap.end() ) {
      std::map< int, int > emptyMap;
      hitProductionMap.insert( std::make_pair( primaryPdg, emptyMap ) );
    }

    if ( hitProductionMap[ primaryPdg ].find( thisPdg ) == hitProductionMap[ primaryPdg ].end() ) {
      hitProductionMap[ primaryPdg ].insert( std::make_pair( thisPdg, 0 ) );
    }

    hitProductionMap[ primaryPdg ][ thisPdg ] ++;

    mcParticleToNHitsMap[ mainMCParticle ]++;
    mcPrimaryToNHitsMap[ primaryMCParticle ]++;
  }

  for ( std::map< int, std::map< int, int > >::iterator it1 = hitProductionMap.begin(); it1 != hitProductionMap.end(); ++it1 ) {
    std::cout << it1->first << std::endl;
    int tot = 0;
    for ( std::map< int, int >::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2 ) {
      std::cout << "  " << it2->first << " : " << it2->second << std::endl;
      tot += it2->second;
    }
    std::cout << "  total : " << tot << std::endl;
  }


  for ( const MCParticle * const pMCParticle : *pMCParticleList ) {
    // Find each cosmic primary MCParticle
    if ( pMCParticle->GetParentList().size() == 0 && ! LArMCParticleHelper::IsNeutrinoInduced( pMCParticle ) ) {
      std::cout << std::string(80, '-') << '\r';
      std::cout << "Primary "   << pMCParticle->GetParticleId() << " " << std::endl;

      int nHits = 0;
      if ( mcPrimaryToNHitsMap.find( pMCParticle ) != mcPrimaryToNHitsMap.end() ) {
        nHits = mcPrimaryToNHitsMap[pMCParticle];
      }
      std::cout << "NHits   : " << nHits << std::endl;
      if ( nHits != 0 ) {
        PrintDaughterInfo( pMCParticle, 1, mcParticleToNHitsMap );
      }
    }
  }
  */
  /* END TESTING ============================================= */





  // ===================================================
  // Use MC information to find the origin of hits
  // ===================================================
  
  // Identify target MCParticles
  MCParticleVector mcNeutrinoVector;
  this->SelectTrueNeutrinos( pMCParticleList, mcNeutrinoVector );

  MCParticleList candidateTargetMCParticleList;
  this->GetCandidateTargetMCParticles( pMCParticleList, mcNeutrinoVector, candidateTargetMCParticleList ); 

  LArMonitoringHelper::MCContributionMap candidateTargetToCaloHitListMap;
  this->GetTargetToCaloHitMatches( pCaloHitList, pMCParticleList, candidateTargetMCParticleList, candidateTargetToCaloHitListMap );
  
  LArMonitoringHelper::MCContributionMap targetToCaloHitListMap;
  LArMonitoringHelper::CaloHitToMCMap    caloHitToTargetMap;
  this->SelectTargetMCParticles( candidateTargetToCaloHitListMap, targetToCaloHitListMap, caloHitToTargetMap);
  
  // Assign each hit an origin
  CaloHitToOriginMap caloHitToOriginMap;
  this->GetHitOrigins( pCaloHitList, caloHitToTargetMap, caloHitToOriginMap );

  // ===================================================
  // Find the purity / significance of the PFOs
  // ===================================================

  PfoToDoubleMap pfoToPurityMap;
  this->GetPurities( pPfoList, caloHitToOriginMap, pfoToPurityMap );
 
  PfoToDoubleMap pfoToSignificanceMap;
  this->GetSignificances( pPfoList, caloHitToOriginMap, targetToCaloHitListMap, caloHitToTargetMap, pfoToSignificanceMap );
  
  MCToIntMap targetIdMap;
  this->GetTargetIds( targetToCaloHitListMap, targetIdMap );

  // ===================================================
  // Classify PFOs based on their cosmic MC information
  // ===================================================

  // Get the list of primary cosmic muon MCParticles
  MCParticleList cosmicMuonList; 
  this->GetPrimaryCosmicMuons( pMCParticleList, cosmicMuonList );

  // Find the PFOs are primary cosmic ray muons
  PfoToBoolMap pfoToIsCosmicMuonMap;
  this->GetIsCosmicMuon( pPfoList, caloHitToOriginMap, pfoToIsCosmicMuonMap);
 
  // ------------------------------------------------------------------------------------------------------------------------------------------
  // Here begins what might actually form the algorithm - everything else is for development only
  // ------------------------------------------------------------------------------------------------------------------------------------------

  // ===================================================
  // Calculate interesting things
  // ===================================================

  PfoToIntMap pfoIdMap;
  this->GetPfoIds( pPfoList, pfoIdMap );

  PfoToPfoListMap pfoAssociationMap;
  this->GetPfoAssociations ( pPfoList, pfoAssociationMap ); 

  PfoToIntMap pfoToSliceIdMap;
  this->SliceEvent( pfoAssociationMap, pfoToSliceIdMap );

  CRCandidateList candidates;
  this->GetCRCandidates( pPfoList, pfoToPurityMap, pfoToSignificanceMap, pfoIdMap, pfoToSliceIdMap, pfoToIsCosmicMuonMap, candidates );


  // ------------------------------------------------------------------------------------------------------------------------------------------


  // ===================================================
  // Output to a root file
  // ===================================================
  
  CaloHitToPfoMap caloHit2DToPfoMap;
  CaloHitToPfoMap caloHit3DToPfoMap;
  this->GetCaloHitToPfoMap( pPfoList, caloHit2DToPfoMap, caloHit3DToPfoMap );


  this->Write2DHits( caloHitToOriginMap, caloHitToTargetMap, caloHit2DToPfoMap, pfoIdMap, targetIdMap );
  this->Write3DHits( caloHit3DToPfoMap, pfoIdMap );
  this->WriteTargets( targetIdMap );
  this->WriteCosmics( cosmicMuonList );
  this->WritePfos( candidates );


  m_eventNumber++;

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::SelectTrueNeutrinos(const MCParticleList *const pAllMCParticleList, MCParticleVector &selectedMCNeutrinoVector) const
{
    MCParticleVector allMCNeutrinoVector;
    LArMCParticleHelper::GetNeutrinoMCParticleList(pAllMCParticleList, allMCNeutrinoVector);

    for (const MCParticle *const pMCNeutrino : allMCNeutrinoVector)
    {
        // ATTN Now demand that input mc neutrinos LArMCParticles, with addition of interaction type
        const LArMCParticle *const pLArMCNeutrino(dynamic_cast<const LArMCParticle*>(pMCNeutrino));

        if (pLArMCNeutrino && (0 != pLArMCNeutrino->GetNuanceCode()))
            selectedMCNeutrinoVector.push_back(pMCNeutrino);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::GetCandidateTargetMCParticles( const MCParticleList *const pAllMCParticleList, MCParticleVector selectedMCNeutrinoVector, MCParticleList &candidateTargetMCParticleList) const
{
  for ( const MCParticle *const pMCParticle : *pAllMCParticleList )
  {
    if ( ! LArMCParticleHelper::IsNeutrinoInduced(pMCParticle) ) continue;
    if ( ! LArMCParticleHelper::IsPrimary(pMCParticle) )         continue;
    if ( std::abs(pMCParticle->GetParticleId()) == NEUTRON )     continue;
    if ( std::find(selectedMCNeutrinoVector.begin(), selectedMCNeutrinoVector.end(), LArMCParticleHelper::GetParentNeutrino(pMCParticle)) != selectedMCNeutrinoVector.end() ) {
      candidateTargetMCParticleList.push_back(pMCParticle);  
    }    

  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::GetTargetToCaloHitMatches( const CaloHitList *const pCaloHitList, const MCParticleList *const pAllMCParticleList, MCParticleList candidateTargetMCParticleList, LArMonitoringHelper::MCContributionMap & candidateTargetToCaloHitListMap ) const
{
  for ( const CaloHit * const pCaloHit : *pCaloHitList )
  {
    const MCParticle * mainMCParticle = this->GetMainMCParticle(pCaloHit);
    if ( ! mainMCParticle )                                         continue; // Ghost hit
    if ( ! LArMCParticleHelper::IsNeutrinoInduced(mainMCParticle) ) continue;

    const MCParticle * primaryMCParticle = LArMCParticleHelper::GetPrimaryMCParticle( mainMCParticle );

    if ( std::find(candidateTargetMCParticleList.begin(), candidateTargetMCParticleList.end(), primaryMCParticle) != candidateTargetMCParticleList.end() ) {

      if ( this->PassMCParticleChecks( primaryMCParticle, primaryMCParticle, mainMCParticle ) ) {

        if ( candidateTargetToCaloHitListMap.find(primaryMCParticle) == candidateTargetToCaloHitListMap.end() ) {
          CaloHitList emptyList;
          candidateTargetToCaloHitListMap.insert(std::make_pair(primaryMCParticle, emptyList));
        }
        candidateTargetToCaloHitListMap[ primaryMCParticle ].push_back( pCaloHit );
      }
    }
 
  }  
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle * CRTaggingAlgorithm::GetMainMCParticle( const pandora::CaloHit * const pCaloHit ) const 
{
    float bestWeight(0.f);
    const MCParticle *pBestMCParticle(nullptr);
    const MCParticleWeightMap &hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

    MCParticleVector mcParticleVector;
    for (const MCParticleWeightMap::value_type &mapEntry : hitMCParticleWeightMap) mcParticleVector.push_back(mapEntry.first);
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

    for (const MCParticle *const pMCParticle : mcParticleVector)
    {
        const float weight(hitMCParticleWeightMap.at(pMCParticle));

        if (weight > bestWeight)
        {
            bestWeight = weight;
            pBestMCParticle = pMCParticle;
        }
    }

    return pBestMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CRTaggingAlgorithm::PassMCParticleChecks(const MCParticle *const pOriginalPrimary, const MCParticle *const pThisMCParticle, const MCParticle *const pHitMCParticle) const
{
    if (NEUTRON == std::abs(pThisMCParticle->GetParticleId()))
        return false;

    if ((PHOTON == pThisMCParticle->GetParticleId()) && (PHOTON != pOriginalPrimary->GetParticleId()) && (E_MINUS != std::abs(pOriginalPrimary->GetParticleId())))
    {
        if ((pThisMCParticle->GetEndpoint() - pThisMCParticle->GetVertex()).GetMagnitude() > m_maxPhotonPropagation)
            return false;
    }

    if (pThisMCParticle == pHitMCParticle)
        return true;

    for (const MCParticle *const pDaughterMCParticle : pThisMCParticle->GetDaughterList())
    {
        if (this->PassMCParticleChecks(pOriginalPrimary, pDaughterMCParticle, pHitMCParticle))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::SelectTargetMCParticles( LArMonitoringHelper::MCContributionMap candidateTargetToCaloHitListMap, LArMonitoringHelper::MCContributionMap & targetToCaloHitListMap, LArMonitoringHelper::CaloHitToMCMap & caloHitToTargetMap) const
{
  LArMonitoringHelper::MCContributionMap::iterator it;
  for (it=candidateTargetToCaloHitListMap.begin(); it != candidateTargetToCaloHitListMap.end(); ++it) {
    int nHits_U = 0;
    int nHits_V = 0;
    int nHits_W = 0;
    for ( const CaloHit * const pCaloHit : it->second ) {
      switch ( pCaloHit->GetHitType() ) {
        case TPC_VIEW_U:
          nHits_U++;
          break;
        case TPC_VIEW_V:
          nHits_V++;
          break;
        case TPC_VIEW_W:
          nHits_W++;
          break;
        default:
          break;
      }
    }

    if ( ( nHits_U + nHits_V + nHits_W ) < m_matchingMinPrimaryHits )  continue;
    
    int nGoodViews = 0;
    if ( nHits_U >= m_matchingMinHitsForGoodView ) ++nGoodViews; 
    if ( nHits_V >= m_matchingMinHitsForGoodView ) ++nGoodViews; 
    if ( nHits_W >= m_matchingMinHitsForGoodView ) ++nGoodViews; 

    if (nGoodViews < m_matchingMinPrimaryGoodViews) continue;
 
    targetToCaloHitListMap.insert(std::make_pair(it->first, it->second));
    for ( const CaloHit * const pCaloHit : it->second ) {
      caloHitToTargetMap.insert(std::make_pair(pCaloHit, it->first));
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::GetHitOrigins( const CaloHitList *const pCaloHitList, LArMonitoringHelper::CaloHitToMCMap caloHitToTargetMap, CaloHitToOriginMap & caloHitToOriginMap ) const
{
  for ( const CaloHit * const pCaloHit : *pCaloHitList ) {
    CaloHitOrigin origin = COSMIC;
    caloHitToOriginMap.insert(std::make_pair(pCaloHit, origin));
 
    if ( caloHitToTargetMap.find( pCaloHit ) != caloHitToTargetMap.end() ) {
      caloHitToOriginMap[pCaloHit] = NEUTRINO_TARGET; 
      continue;
    }
    
    const MCParticle * mainMCParticle = this->GetMainMCParticle(pCaloHit);
    if ( ! mainMCParticle ) {
      caloHitToOriginMap[pCaloHit] = GHOST; 
      continue;
    }

    if ( LArMCParticleHelper::IsNeutrinoInduced(mainMCParticle) ) {
      caloHitToOriginMap[pCaloHit] = NEUTRINO_NON_TARGET; 
      continue;
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::GetPurities( const PfoList * const pPfoList, CaloHitToOriginMap caloHitToOriginMap, PfoToDoubleMap & pfoToPurityMap ) const 
{
  for ( const ParticleFlowObject * const pPfo : *pPfoList ) {
    double purity = this->GetPurity( pPfo, caloHitToOriginMap );
    pfoToPurityMap.insert(std::make_pair( pPfo, purity ));
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

double CRTaggingAlgorithm::GetPurity( const ParticleFlowObject * const pPfo, CaloHitToOriginMap caloHitToOriginMap ) const
{
  int nHits       = 0;
  int nTargetHits = 0;

  ClusterList clusters2D;
  LArPfoHelper::GetTwoDClusterList(pPfo, clusters2D);

  for (const Cluster * const pCluster : clusters2D){
    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit * const pCaloHit : caloHitList){
      if (caloHitToOriginMap[pCaloHit] == NEUTRINO_TARGET) nTargetHits++;
      nHits++;
    }
  }

  double purity = double (nTargetHits) / double (nHits);
  
  return purity;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::GetSignificances( const PfoList * const pPfoList, CaloHitToOriginMap caloHitToOriginMap, LArMonitoringHelper::MCContributionMap & targetToCaloHitListMap, LArMonitoringHelper::CaloHitToMCMap caloHitToTargetMap, PfoToDoubleMap & pfoToSignificanceMap ) const
{
  for ( const ParticleFlowObject * const pPfo : *pPfoList ) {
    double significance = this->GetSignificance( pPfo, caloHitToOriginMap, targetToCaloHitListMap, caloHitToTargetMap );
    pfoToSignificanceMap.insert(std::make_pair( pPfo, significance ));
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

double CRTaggingAlgorithm::GetSignificance( const ParticleFlowObject * const pPfo, CaloHitToOriginMap caloHitToOriginMap, LArMonitoringHelper::MCContributionMap & targetToCaloHitListMap, LArMonitoringHelper::CaloHitToMCMap caloHitToTargetMap ) const
{
  double significance = 0;

  for ( LArMonitoringHelper::MCContributionMap::iterator it = targetToCaloHitListMap.begin(); it != targetToCaloHitListMap.end(); ++it ) {
    int nSharedHits = 0;
    int nTargetHits = it->second.size();

    ClusterList clusters2D;
    LArPfoHelper::GetTwoDClusterList(pPfo, clusters2D);
  
    for (const Cluster * const pCluster : clusters2D){
      CaloHitList caloHitList;
      pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
  
      for (const CaloHit * const pCaloHit : caloHitList){
        if (caloHitToOriginMap[pCaloHit] == NEUTRINO_TARGET) {
          if ( caloHitToTargetMap[pCaloHit] == it->first ) nSharedHits++;
        }
      }
    }

    significance += double (nSharedHits) / double (nTargetHits) ;
  }

  return significance;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void CRTaggingAlgorithm::GetPfoIds( const PfoList * const pPfoList, PfoToIntMap & pfoIdMap ) const
{
  int thisId = 0;
  for ( const ParticleFlowObject * const pPfo : *pPfoList ) {
    pfoIdMap.insert(std::make_pair(pPfo, thisId));
    thisId++;
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void CRTaggingAlgorithm::GetTargetIds( LArMonitoringHelper::MCContributionMap targetToCaloHitListMap, MCToIntMap & targetIdMap ) const
{
  int thisId = 0;
  for ( LArMonitoringHelper::MCContributionMap::iterator it = targetToCaloHitListMap.begin(); it != targetToCaloHitListMap.end(); ++it) {
    targetIdMap.insert(std::make_pair(it->first, thisId));
    thisId++;
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::GetPrimaryCosmicMuons( const MCParticleList *const pMCParticleList, MCParticleList & cosmicMuonList ) const {

  for ( const MCParticle * const pMCParticle : *pMCParticleList ) {
    if ( std::abs( pMCParticle->GetParticleId() ) == MU_MINUS && LArMCParticleHelper::IsPrimary( pMCParticle ) ) {
      cosmicMuonList.push_back( pMCParticle );
    }
  }

}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::GetIsCosmicMuon( const PfoList * const pPfoList, CaloHitToOriginMap caloHitToOriginMap, PfoToBoolMap & pfoToIsCosmicMuonMap) const {

  for ( const ParticleFlowObject * const pPfo : *pPfoList ) {
    // Get the list of associated MCParticles
    std::map< const MCParticle * const, int > mcParticleContributionMap;

    // For each hit in the PFO, find the associated cosmic MC primary (if any) 
    // count the number of hits associated with each such MCParticle in the PFO
    ClusterList clusters2D;
    LArPfoHelper::GetTwoDClusterList(pPfo, clusters2D);

    int nHits = 0;
    for (const Cluster * const pCluster : clusters2D){
      CaloHitList caloHitList;
      pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

      for (const CaloHit * const pCaloHit : caloHitList){
        nHits++;
        if ( caloHitToOriginMap[ pCaloHit ] == COSMIC ) { 
          const MCParticle * mainMCParticle = this->GetMainMCParticle(pCaloHit);
          if ( mainMCParticle != nullptr ) {
            
            if ( mcParticleContributionMap.find( mainMCParticle ) != mcParticleContributionMap.end() ) {
              mcParticleContributionMap.insert( std::make_pair( mainMCParticle, 0 ) );
            } 
            mcParticleContributionMap[ mainMCParticle ]++;

          }
        }
      }
    }

    bool isCosmicMuon = false;

    // Now find the MC primary with the most associated hits
    if ( mcParticleContributionMap.size() != 0 ){
      const MCParticle * mainMCParticle = nullptr;
      int maxHits = 0;
      for ( std::map< const MCParticle * const, int >::iterator it = mcParticleContributionMap.begin(); it != mcParticleContributionMap.end(); ++it ){
        if ( it->second > maxHits ) {
          mainMCParticle = it->first;
          maxHits        = it->second;
        }
      } 

      // Determine if is a primary cosmic ray muon
      isCosmicMuon = ( std::abs( mainMCParticle->GetParticleId() ) == MU_MINUS && LArMCParticleHelper::IsPrimary( mainMCParticle ) ) ;

    }
    
    pfoToIsCosmicMuonMap.insert( std::make_pair( pPfo, isCosmicMuon ) );

  }

}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::GetPfoAssociations( const PfoList * const pPfoList, PfoToPfoListMap & pfoAssociationMap ) const 
{
  unsigned int minimum3DHits = 15;

  for ( const ParticleFlowObject * const pPfo1 : *pPfoList ) {
    if ( pfoAssociationMap.find( pPfo1 ) == pfoAssociationMap.end() ) {
      PfoList blankPfoList;
      pfoAssociationMap.insert(std::make_pair( pPfo1, blankPfoList ));
    }

    ClusterList clusters3D1;
    LArPfoHelper::GetThreeDClusterList(pPfo1, clusters3D1);
    if ( clusters3D1.size() == 0 ) continue;
    const Cluster * const pCluster1 = clusters3D1.front();
    CaloHitList caloHitList1;
    pCluster1->GetOrderedCaloHitList().FillCaloHitList(caloHitList1);
    if ( caloHitList1.size() < minimum3DHits ) continue;

 
    ThreeDSlidingFitResult fit1 (pCluster1, 5, LArGeometryHelper::GetWireZPitch(this->GetPandora())); 
    CartesianVector min1    = fit1.GetGlobalMinLayerPosition();
    CartesianVector max1    = fit1.GetGlobalMaxLayerPosition();
    
    ThreeDSlidingFitResult fit1b (pCluster1, 100, LArGeometryHelper::GetWireZPitch(this->GetPandora())); 
    CartesianVector minDir1 = fit1b.GetGlobalMinLayerDirection() * ( -1 );
    CartesianVector maxDir1 = fit1b.GetGlobalMaxLayerDirection();

    for ( const ParticleFlowObject * const pPfo2 : *pPfoList ) {
      if ( pPfo1 != pPfo2 ) {
        if ( pfoAssociationMap.find( pPfo2 ) == pfoAssociationMap.end() ) {
          PfoList blankPfoList;
          pfoAssociationMap.insert(std::make_pair( pPfo2, blankPfoList ));
        }

        ClusterList clusters3D2;
        LArPfoHelper::GetThreeDClusterList(pPfo2, clusters3D2);
        if ( clusters3D2.size() == 0 ) continue;
        const Cluster * const pCluster2 = clusters3D2.front();
        CaloHitList caloHitList2;
        pCluster2->GetOrderedCaloHitList().FillCaloHitList(caloHitList2);
        if ( caloHitList2.size() < minimum3DHits ) continue;

        ThreeDSlidingFitResult fit2 (pCluster2, 5, LArGeometryHelper::GetWireZPitch(this->GetPandora())); 
        CartesianVector min2    = fit2.GetGlobalMinLayerPosition();
        CartesianVector max2    = fit2.GetGlobalMaxLayerPosition();
        
        ThreeDSlidingFitResult fit2b (pCluster2, 100, LArGeometryHelper::GetWireZPitch(this->GetPandora())); 
        CartesianVector minDir2 = fit2b.GetGlobalMinLayerDirection() * ( -1 );
        CartesianVector maxDir2 = fit2b.GetGlobalMaxLayerDirection();

        bool isAssociated = false;

        if ( this->CheckAssociation( min1, minDir1, min2, minDir2 ) ) isAssociated = true;
        if ( this->CheckAssociation( min1, minDir1, max2, maxDir2 ) ) isAssociated = true;
        if ( this->CheckAssociation( max1, maxDir1, min2, minDir2 ) ) isAssociated = true;
        if ( this->CheckAssociation( max1, maxDir1, max2, maxDir2 ) ) isAssociated = true;

        if ( isAssociated ) {
          if ( std::find( pfoAssociationMap[ pPfo1 ].begin(), pfoAssociationMap[ pPfo1 ].end(), pPfo2 ) == pfoAssociationMap[ pPfo1 ].end() ) {
            pfoAssociationMap[ pPfo1 ].push_back( pPfo2 );
          }
          if ( std::find( pfoAssociationMap[ pPfo2 ].begin(), pfoAssociationMap[ pPfo2 ].end(), pPfo1 ) == pfoAssociationMap[ pPfo2 ].end() ) {
            pfoAssociationMap[ pPfo2 ].push_back( pPfo1 );
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CRTaggingAlgorithm::CheckAssociation( CartesianVector endPoint1, CartesianVector endDir1, CartesianVector endPoint2, CartesianVector endDir2 ) const 
{
  
  CartesianVector n = endDir1 * ( 1. / endDir1.GetMagnitude() );
  CartesianVector m = endDir2 * ( 1. / endDir2.GetMagnitude() );
 
  CartesianVector a = endPoint2 - endPoint1;
  
  double b = n.GetDotProduct(m);
 
  if ( b == 1 ) return false;
 
  // Distance from endPoint1 along endDir1 to the point of closest approach
  double lambda = ( n - m*b ).GetDotProduct( a ) / ( 1 - b*b );
  
  // Distance from endPoint2 along endDir2 to the point of closest approach
  double mu     = ( n*b - m ).GetDotProduct( a ) / ( 1 - b*b );


  // Calculate the maximum vertex uncertainty
  double deltaTheta = m_angularUncertainty * std::asin(1) / 180;
  double maxVertexUncertainty = m_maxAssociationDist * std::sin( deltaTheta ) + m_positionalUncertainty; 

  // Ensure that the distances to the point of closest approch are within the limits
  if ( lambda < -maxVertexUncertainty )                       return false;
  if ( mu     < -maxVertexUncertainty )                       return false;
  if ( lambda > m_maxAssociationDist + maxVertexUncertainty ) return false;
  if ( mu     > m_maxAssociationDist + maxVertexUncertainty ) return false;

  // Ensure the point of closest approach is within the detector volume 
  // (Because we only want to associate neutrinos we can use the x-direction)
  double face_Xa =    0.00;  // Anode      X face
  double face_Xc =  256.00;  // Cathode    X face
  double face_Yb = -116.25;  // Bottom     Y face
  double face_Yt =  116.25;  // Top        Y face
  double face_Zu =    0.00;  // Upstream   Z face
  double face_Zd = 1036.80;  // Downstream Z face

  CartesianVector impactPos = (endPoint1 + n*lambda + endPoint2 + m*mu ) * 0.5 ;
  double impactX = impactPos.GetX();
  double impactY = impactPos.GetY();
  double impactZ = impactPos.GetZ();

  if ( impactX < face_Xa - maxVertexUncertainty ) return false;
  if ( impactX > face_Xc + maxVertexUncertainty ) return false;
  if ( impactY < face_Yb - maxVertexUncertainty ) return false;
  if ( impactY > face_Yt + maxVertexUncertainty ) return false;
  if ( impactZ < face_Zu - maxVertexUncertainty ) return false;
  if ( impactZ > face_Zd + maxVertexUncertainty ) return false;

  // Maximum uncertainty in impact parameter
  double maxImpactDist = std::sin( deltaTheta ) * ( std::abs(mu) + std::abs(lambda) ) + m_positionalUncertainty;

  // Calculate the distance of closest approach
  CartesianVector D = a - n*lambda + m*mu;

  return ( D.GetMagnitude() <= maxImpactDist );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::SliceEvent( PfoToPfoListMap pfoAssociationMap, PfoToIntMap & pfoToSliceIdMap ) const
{
  std::vector< PfoList > slices;

  for ( PfoToPfoListMap::iterator it = pfoAssociationMap.begin(); it != pfoAssociationMap.end(); ++it ) {
    const ParticleFlowObject * const pPfo = it->first;
    
    bool inSlice = false;
    for ( PfoList & slice : slices ) {
      if ( std::find(slice.begin(), slice.end(), pPfo) != slice.end() ) {
        inSlice = true;
        break;
      }
    }

    if ( ! inSlice ) {
      PfoList newSlice;
      this->FillSlice( pfoAssociationMap, pPfo, newSlice );
      slices.push_back( newSlice );
    }
  }

  int sliceId = 0;
  for ( PfoList slice : slices ) {
    for ( const ParticleFlowObject * const pPfo : slice ) {
      pfoToSliceIdMap.insert( std::make_pair( pPfo, sliceId ) );
    }
    sliceId++;
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::FillSlice( PfoToPfoListMap pfoAssociationMap, const ParticleFlowObject * const pPfo, PfoList & slice ) const
{

  if ( std::find( slice.begin(), slice.end(), pPfo ) == slice.end() ) {
    slice.push_back( pPfo );
    for ( const ParticleFlowObject * const pPfo2 : pfoAssociationMap[ pPfo ] ) {
      this->FillSlice( pfoAssociationMap, pPfo2, slice );
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::GetCRCandidates( const PfoList * const pPfoList, PfoToDoubleMap pfoToPurityMap, PfoToDoubleMap pfoToSignificanceMap, PfoToIntMap pfoIdMap, PfoToIntMap pfoToSliceIdMap, PfoToBoolMap pfoToIsCosmicMuonMap, CRCandidateList & candidates ) const
{
  for ( const ParticleFlowObject * const pPfo : *pPfoList ) {
    CRCandidate candidate(this, pPfo, pfoIdMap[pPfo], pfoToSliceIdMap[pPfo], pfoToPurityMap[pPfo], pfoToSignificanceMap[pPfo], pfoToIsCosmicMuonMap[pPfo]);
    candidates.push_back(candidate);
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::GetCaloHitToPfoMap( const PfoList * const pPfoList, CaloHitToPfoMap & caloHit2DToPfoMap, CaloHitToPfoMap & caloHit3DToPfoMap ) const 
{
  for ( const ParticleFlowObject * const pPfo : *pPfoList ) {
    ClusterList clusters2D;
    LArPfoHelper::GetTwoDClusterList(pPfo, clusters2D);
      
    for (const Cluster * const pCluster : clusters2D){
      CaloHitList caloHitList;
      pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
      for (const CaloHit * const pCaloHit : caloHitList){
        caloHit2DToPfoMap.insert(std::make_pair(pCaloHit, pPfo));
      }
    }

    ClusterList clusters3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);
      
    for (const Cluster * const pCluster : clusters3D){
      CaloHitList caloHitList;
      pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
      for (const CaloHit * const pCaloHit : caloHitList){
        caloHit3DToPfoMap.insert(std::make_pair(pCaloHit, pPfo));
      }
    }

  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::Write2DHits( CaloHitToOriginMap caloHitToOriginMap, LArMonitoringHelper::CaloHitToMCMap caloHitToTargetMap, CaloHitToPfoMap caloHitToPfoMap, PfoToIntMap pfoIdMap, MCToIntMap targetIdMap ) const 
{
  for ( CaloHitToOriginMap::iterator it = caloHitToOriginMap.begin(); it != caloHitToOriginMap.end(); ++it ) {
      const CaloHit * const pCaloHit = it->first;

      int viewInt;
      switch ( pCaloHit->GetHitType() ) { 
        case TPC_VIEW_U:
          viewInt = 0;
          break;   
        case TPC_VIEW_V:
          viewInt = 1;
          break;   
        case TPC_VIEW_W:
          viewInt = 2;
          break;   
        case TPC_3D:
          viewInt = 3;
          break;
        default:
          viewInt = -1;
      }

      if ( viewInt == 0 || viewInt == 1 || viewInt == 2 ){
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "FileId"      , m_fileId       ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "EventId"     , m_eventNumber  ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "View"        , viewInt ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "X"           , pCaloHit->GetPositionVector().GetX() ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "Y"           , pCaloHit->GetPositionVector().GetY() ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "Z"           , pCaloHit->GetPositionVector().GetZ() ));

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "InputEnergy" , pCaloHit->GetInputEnergy()           ));
  
        if ( caloHitToPfoMap.find(pCaloHit) == caloHitToPfoMap.end() ){
          PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "PfoId", -1));  
        } 
        else{ 
          PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "PfoId", pfoIdMap[caloHitToPfoMap[pCaloHit]]));
        }
  
        if ( caloHitToTargetMap.find(pCaloHit) == caloHitToTargetMap.end() ) {
          PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "TargetId" , -1));  
          PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "TargetPdg", std::numeric_limits<int>::max()));  
        }
        else{
          PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "TargetId" , targetIdMap[caloHitToTargetMap[pCaloHit]] ));  
          PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "TargetPdg", caloHitToTargetMap[pCaloHit]->GetParticleId() ));  
        }
  
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "Origin" , int (it->second) ));
      
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), "Hits"));
      }
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::Write3DHits( CaloHitToPfoMap caloHitToPfoMap, PfoToIntMap pfoIdMap ) const 
{
  for ( CaloHitToPfoMap::iterator it = caloHitToPfoMap.begin(); it != caloHitToPfoMap.end(); ++it ) {
      const CaloHit * const pCaloHit = it->first;

      int viewInt;
      switch ( pCaloHit->GetHitType() ) { 
        case TPC_VIEW_U:
          viewInt = 0;
          break;   
        case TPC_VIEW_V:
          viewInt = 1;
          break;   
        case TPC_VIEW_W:
          viewInt = 2;
          break;   
        case TPC_3D:
          viewInt = 3;
          break;
        default:
          viewInt = -1;
      }

      if ( viewInt == 3 ){
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "FileId"      , m_fileId       ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "EventId"     , m_eventNumber  ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "View"        , viewInt ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "X"           , pCaloHit->GetPositionVector().GetX() ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "Y"           , pCaloHit->GetPositionVector().GetY() ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "Z"           , pCaloHit->GetPositionVector().GetZ() ));

        if ( caloHitToPfoMap.find(pCaloHit) == caloHitToPfoMap.end() ){
          PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "PfoId", -1));  
        } 
        else{ 
          PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "PfoId", pfoIdMap[caloHitToPfoMap[pCaloHit]]));
        }
          
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "TargetId" , -1));  
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "TargetPdg", std::numeric_limits<int>::max()));  
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "Origin"   , -1));
  
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), "Hits"));
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::WriteTargets( MCToIntMap targetIdMap ) const 
{
  for ( MCToIntMap::iterator it=targetIdMap.begin(); it != targetIdMap.end(); ++it ){
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Targets", "FileId"      , m_fileId       ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Targets", "EventId"     , m_eventNumber  ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Targets", "Id"          , it->second ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Targets", "Pdg"         , it->first->GetParticleId()  ));
    
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "Targets"));
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::WriteCosmics( MCParticleList cosmicMuonList ) const
{
  for ( const MCParticle * const pMCParticle : cosmicMuonList ) {
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "FileId"      , m_fileId                          ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "EventId"     , m_eventNumber                     ));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "Energy"      , pMCParticle->GetEnergy()          ));
    
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "StartX"      , pMCParticle->GetVertex().GetX()  ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "StartY"      , pMCParticle->GetVertex().GetY()  ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "StartZ"      , pMCParticle->GetVertex().GetZ()  ));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "EndX"        , pMCParticle->GetEndpoint().GetX()  ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "EndY"        , pMCParticle->GetEndpoint().GetY()  ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "EndZ"        , pMCParticle->GetEndpoint().GetZ()  ));

    double P  = pMCParticle->GetMomentum().GetMagnitude();
    double Px = pMCParticle->GetMomentum().GetX();
    double Py = pMCParticle->GetMomentum().GetY();
    double Pz = pMCParticle->GetMomentum().GetZ();

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "Momentum"    , P  ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "MomentumX"   , Px ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "MomentumY"   , Py ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "MomentumZ"   , Pz ));

    double thetaX = std::acos( Px / P );
    double thetaY = std::acos( Py / P );
    double thetaZ = std::acos( Pz / P );

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "ThetaX"      , thetaX ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "ThetaY"      , thetaY ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CosmicMuons", "ThetaZ"      , thetaZ ));
    
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "CosmicMuons"));

  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::WritePfos( CRCandidateList candidates ) const
{
  for ( CRCandidate candidate : candidates ){
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "FileId"      , m_fileId                   ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "EventId"     , m_eventNumber              ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Id"          , candidate.m_id             ));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "N2DHits"     , candidate.m_n2DHits        ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "N3DHits"     , candidate.m_n3DHits        ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "NDaughters"  , candidate.m_nDaughters     ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Slice"       , candidate.m_slice          ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "TotalEnergy" , candidate.m_totalEnergy    ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "MeanEnergy"  , candidate.m_meanEnergy     ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "CanFit"      , candidate.m_canFit ? 1 : 0 ));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "X1", candidate.m_X1 ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Y1", candidate.m_Y1 ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Z1", candidate.m_Z1 ));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "X2", candidate.m_X2 ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Y2", candidate.m_Y2 ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Z2", candidate.m_Z2 ));
    
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Length", candidate.m_length       ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Curvature", candidate.m_curvature ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "FitRMS", candidate.m_fitRMS       ));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Purity"       , candidate.m_purity       ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Significance" , candidate.m_significance ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Class"        , int(candidate.m_class)   ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "IsCosmicMuon" , ( candidate.m_isCosmicMuon ? 1 : 0 ) ));
    
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "PFOs"));
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CRTaggingAlgorithm::CRCandidate::CRCandidate(const CRTaggingAlgorithm * const algorithm, const pandora::ParticleFlowObject * const pPfo, int id, int slice, double purity, double significance, bool isCosmicMuon ) :
  m_algorithm ( algorithm ), 
  m_pPfo( pPfo ),
  m_id( id ),
  m_slice( slice ),
  m_n2DHits( -1 ),
  m_n3DHits( -1 ),
  m_nDaughters( -1 ),
  m_totalEnergy( -1 ),
  m_meanEnergy( -1 ),
  m_canFit( false ),
  m_fitRMS ( std::numeric_limits<double>::max() ),
  m_X1 ( std::numeric_limits<double>::max() ),
  m_Y1 ( std::numeric_limits<double>::max() ),
  m_Z1 ( std::numeric_limits<double>::max() ),
  m_X2 ( std::numeric_limits<double>::max() ),
  m_Y2 ( std::numeric_limits<double>::max() ),
  m_Z2 ( std::numeric_limits<double>::max() ),
  m_length ( std::numeric_limits<double>::max() ),
  m_curvature ( std::numeric_limits<double>::max() ),
  m_purity( purity ),
  m_significance( significance ),
  m_isCosmicMuon ( isCosmicMuon )
{
  // Get the hit lists
  ClusterList clusters2D;
  LArPfoHelper::GetTwoDClusterList(pPfo, clusters2D);
    
  for (const Cluster * const pCluster : clusters2D){
    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
    for (const CaloHit * const pCaloHit : caloHitList){
      m_hitList2D.push_back(pCaloHit);
    }
  }

  ClusterList clusters3D;
  LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);

  for (const Cluster * const pCluster : clusters3D){
    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
    for (const CaloHit * const pCaloHit : caloHitList){
      m_hitList3D.push_back(pCaloHit);
    }
  }


  // NB. The order of these calculations matters!
  //     Some calculations depend on the results of the previous

  this->CalculateNHits();
  this->CalculateNDaughters();
  this->CalculateTotalEnergy();
  this->CalculateMeanEnergy();
  this->CalculateFitVariables();  

  this->DetermineClass();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::CRCandidate::CalculateNHits()
{
  m_n2DHits = m_hitList2D.size();
  m_n3DHits = m_hitList3D.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::CRCandidate::CalculateNDaughters()
{
  m_nDaughters = m_pPfo->GetNDaughterPfos();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::CRCandidate::CalculateTotalEnergy()
{
  m_totalEnergy = 0;
  for ( const CaloHit * const pCaloHit : m_hitList2D ) {
    if ( pCaloHit->GetHitType() == TPC_VIEW_W ){
      m_totalEnergy += pCaloHit->GetInputEnergy();
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::CRCandidate::CalculateMeanEnergy()
{
  m_meanEnergy = m_totalEnergy / ( double (m_n2DHits) );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::CRCandidate::CalculateFitVariables()
{
  int minimum3DHits = 15;

  m_canFit = ( m_n3DHits >= minimum3DHits );
  if ( ! m_canFit ) return;

  ClusterList clusters3D;
  LArPfoHelper::GetThreeDClusterList(m_pPfo, clusters3D);
  const Cluster * const pCluster = clusters3D.front();

  ThreeDSlidingFitResult fit (pCluster, 5, LArGeometryHelper::GetWireZPitch(m_algorithm->GetPandora())); 

  CartesianVector minPos = fit.GetGlobalMinLayerPosition();
  CartesianVector maxPos = fit.GetGlobalMaxLayerPosition();
  
  // End Points 
  m_X1 = minPos.GetX();
  m_Y1 = minPos.GetY();
  m_Z1 = minPos.GetZ();
 
  m_X2 = maxPos.GetX();
  m_Y2 = maxPos.GetY();
  m_Z2 = maxPos.GetZ();

  // Straight line length between end points
  m_length = (maxPos - minPos).GetMagnitude();

  //RMS
  m_fitRMS = fit.GetMinLayerRms();

  // Curvature
  int minLayer = std::max(0, fit.GetMinLayer());
  int maxLayer = fit.GetMaxLayer();

  unsigned int nSamples = maxLayer - minLayer;
  
  m_curvature = 0;

  // Calculate the mean direction
  CartesianVector meanDir(0, 0, 0);
  for ( unsigned int i=0; i<nSamples; i++ ) {
    CartesianVector dir (0, 0, 0);
    fit.GetGlobalFitDirection( ( minLayer + i ) * LArGeometryHelper::GetWireZPitch(m_algorithm->GetPandora()), dir);
    meanDir += dir * ( 1/ double ( nSamples ));
  }

  // Now calculate the average deviation from the mean direction
  for ( unsigned int i=0; i<nSamples; i++ ) {
    CartesianVector dir (0, 0, 0);
    fit.GetGlobalFitDirection( ( minLayer + i ) * LArGeometryHelper::GetWireZPitch(m_algorithm->GetPandora()), dir);
    m_curvature += (dir - meanDir).GetMagnitude() / double ( nSamples );
  }

}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTaggingAlgorithm::CRCandidate::DetermineClass()
{
  int    minimum2DHits        = 15;
  int    minimum3DHits        = 15;
  double pureThreshold        = 0.95;
  double impureThreshold      = 0.05;
  double significantThreshold = 0.1;

  if ( m_n2DHits < minimum2DHits || m_n3DHits < minimum3DHits ) {
    m_class = UNCLASSIFIABLE; 
  }
  else{
    if ( m_purity > pureThreshold ) {
      if ( m_significance >= significantThreshold ) {
        m_class = CLEAR_NEUTRINO;
      }
      else {
        m_class = FRAGMENTED;
      }
    }
    else if ( m_purity < impureThreshold ) {
      if ( m_significance >= significantThreshold ) {
        m_class = ABSORBED;
      }
      else {
        m_class = CLEAR_COSMIC;
      }
    }
    else {
      m_class = MIXED;
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CRTaggingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName"       , m_inputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName"   , m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileId"            , m_fileId));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxPhotonPropagation", m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MatchingMinPrimaryHits", m_matchingMinPrimaryHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MatchingMinHitsForGoodView", m_matchingMinHitsForGoodView));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MatchingMinPrimaryGoodViews", m_matchingMinPrimaryGoodViews));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AngularUncertainty"   , m_angularUncertainty ));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PositionalUncertainty", m_positionalUncertainty ));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxAssociationDist"   , m_maxAssociationDist));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
