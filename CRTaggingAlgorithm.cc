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

using namespace pandora;

namespace lar_content
{

CRTaggingAlgorithm::CRTaggingAlgorithm() :
  m_eventNumber(0),
  m_maxPhotonPropagation(2.5f),
  m_matchingMinPrimaryHits(15),
  m_matchingMinHitsForGoodView(5),
  m_matchingMinPrimaryGoodViews(2)
{
}

CRTaggingAlgorithm::~CRTaggingAlgorithm() {
  PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "PFOs"   , "TaggedPFOs.root", "UPDATE"));
  PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "Hits"   , "TaggedPFOs.root", "UPDATE"));
  PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "Targets", "TaggedPFOs.root", "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

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
 
  // ------------------------------------------------------------------------------------------------------------------------------------------
  // Here begins what might actually form the algorithm - everything else is for development only
  // ------------------------------------------------------------------------------------------------------------------------------------------

  // ===================================================
  // Calculate interesting things
  // ===================================================

  PfoToIntMap pfoIdMap;
  this->GetPfoIds( pPfoList, pfoIdMap );

  CRCandidateList candidates;
  this->GetCRCandidates( pPfoList, pfoToPurityMap, pfoToSignificanceMap, pfoIdMap, candidates );

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

void CRTaggingAlgorithm::GetCRCandidates( const PfoList * const pPfoList, PfoToDoubleMap pfoToPurityMap, PfoToDoubleMap pfoToSignificanceMap, PfoToIntMap pfoIdMap, CRCandidateList & candidates ) const
{
  for ( const ParticleFlowObject * const pPfo : *pPfoList ) {
    CRCandidate candidate(this, pPfo, pfoIdMap[pPfo], pfoToPurityMap[pPfo], pfoToSignificanceMap[pPfo] );
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

void CRTaggingAlgorithm::WritePfos( CRCandidateList candidates ) const
{
  for ( CRCandidate candidate : candidates ){
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "FileId"      , m_fileId                   ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "EventId"     , m_eventNumber              ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Id"          , candidate.m_id             ));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "N2DHits"     , candidate.m_n2DHits        ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "N3DHits"     , candidate.m_n3DHits        ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "TotalEnergy" , candidate.m_totalEnergy    ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "MeanEnergy"  , candidate.m_meanEnergy     ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "CanFit"      , candidate.m_canFit ? 1 : 0 ));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "X1", candidate.m_X1 ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Y1", candidate.m_Y1 ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Z1", candidate.m_Z1 ));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "X2", candidate.m_X2 ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Y2", candidate.m_Y2 ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Z2", candidate.m_Z2 ));
    
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Length", candidate.m_length ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "FitRMS", candidate.m_fitRMS ));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Purity"      , candidate.m_purity       ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Significance", candidate.m_significance ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "PFOs", "Class"       , int(candidate.m_class)   ));
    
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "PFOs"));
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

CRTaggingAlgorithm::CRCandidate::CRCandidate(const CRTaggingAlgorithm * const algorithm, const pandora::ParticleFlowObject * const pPfo, int id, double purity, double significance ) :
  m_algorithm ( algorithm ), 
  m_pPfo( pPfo ),
  m_id( id ),
  m_n2DHits( -1 ),
  m_n3DHits( -1 ),
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
  m_purity( purity ),
  m_significance( significance ) 
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

void CRTaggingAlgorithm::CRCandidate::CalculateTotalEnergy()
{
  m_totalEnergy = 0;
  for ( const CaloHit * const pCaloHit : m_hitList2D ) {
    m_totalEnergy += pCaloHit->GetInputEnergy();
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

  ThreeDSlidingFitResult fit (pCluster, 10000, LArGeometryHelper::GetWireZPitch(m_algorithm->GetPandora())); 

  CartesianVector minPos = fit.GetGlobalMinLayerPosition();
  CartesianVector maxPos = fit.GetGlobalMaxLayerPosition();
  
  // End Points 
  m_X1 = minPos.GetX();
  m_Y1 = minPos.GetY();
  m_Z1 = minPos.GetZ();
 
  m_X2 = maxPos.GetX();
  m_Y2 = maxPos.GetY();
  m_Z2 = maxPos.GetZ();

  // Straight line length
  m_length = (maxPos - minPos).GetMagnitude();

  //RMS
  m_fitRMS = fit.GetMinLayerRms();
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

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
