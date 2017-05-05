/**
 *  @file   larpandoracontent/LArCRTagging/HitDumpingAlgorithm.cc
 * 
 *  @brief  Implementation of the cosmic ray tagging algorithm class.
 * 
 *  $Log: $
 */

#include "larpandoracontent/LArCRTagging/HitDumpingAlgorithm.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

using namespace pandora;

namespace lar_content
{

HitDumpingAlgorithm::HitDumpingAlgorithm() :
  m_eventNumber(0)
{}

HitDumpingAlgorithm::~HitDumpingAlgorithm() {
  PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "Hits"   , "PostHitData.root", "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitDumpingAlgorithm::Run()
{

  // ===================================================
  //  the input collections
  // ===================================================

  //  the list of all the CaloHits (2D)
  const CaloHitList *pCaloHitList = nullptr;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList)); 

  for ( const CaloHit * const pCaloHit : *pCaloHitList ) {
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
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "Hits", "Z"           , pCaloHit->GetPositionVector().GetZ() ));
      PANDORA_MONITORING_API(FillTree(this->GetPandora(), "Hits"));
    }
  }

  m_eventNumber++;
  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitDumpingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName"   , m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileId"            , m_fileId));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
