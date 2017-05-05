/**
 *  @file   larpandoracontent/LArCRTagging/HitDumpingAlgorithm.h
 * 
 *  @brief  Header file for the hit matching algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_HIT_DUMPING_ALGORITHM_H
#define LAR_HIT_DUMPING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "PandoraMonitoringApi.h"

namespace lar_content
{

/**
 *  @brief  HitDumpingAlgorithm class
 */
class HitDumpingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief  Default constructor
     */
    HitDumpingAlgorithm();

    ~HitDumpingAlgorithm();

private:

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_caloHitListName;              ///< The name of the input CaloHit list

    int                     m_eventNumber;
    int                     m_fileId;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *HitDumpingAlgorithm::Factory::CreateAlgorithm() const
{
    return new HitDumpingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CR_TAGGING_ALGORITHM_H
