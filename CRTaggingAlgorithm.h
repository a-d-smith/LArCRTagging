/**
 *  @file   larpandoracontent/LArCRTagging/CRTaggingAlgorithm.h
 * 
 *  @brief  Header file for the pfo characterisation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CR_TAGGING_ALGORITHM_H
#define LAR_CR_TAGGING_ALGORITHM_H 1

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

enum CaloHitOrigin {
  NEUTRINO_TARGET,
  NEUTRINO_NON_TARGET,
  COSMIC, 
  GHOST
};

enum CRTagClass {
  CLEAR_COSMIC,
  CLEAR_NEUTRINO,
  FRAGMENTED,
  ABSORBED,
  MIXED,
  UNCLASSIFIABLE
};


typedef std::map<const pandora::CaloHit * const, CaloHitOrigin >      CaloHitToOriginMap;

typedef std::map<const pandora::MCParticle * const, int >             MCToIntMap;

typedef std::map<const pandora::ParticleFlowObject * const, double >           PfoToDoubleMap;
typedef std::map<const pandora::ParticleFlowObject * const, int >              PfoToIntMap;
typedef std::map<const pandora::ParticleFlowObject * const, bool >             PfoToBoolMap;
typedef std::map<const pandora::ParticleFlowObject * const, pandora::PfoList > PfoToPfoListMap;

typedef std::map<const pandora::ParticleFlowObject * const, pandora::CaloHitList >            PfoToCaloHitListMap;
typedef std::map<const pandora::CaloHit * const, const pandora::ParticleFlowObject * const >  CaloHitToPfoMap;

/**
 *  @brief  CRTaggingAlgorithm class
 */
class CRTaggingAlgorithm : public pandora::Algorithm
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
    CRTaggingAlgorithm();

    ~CRTaggingAlgorithm();

private:

    /**
     *  @brief  Class to encapsulate the logic required determine if a PFO should or shouldn't be tagged as a cosmic ray
     */
    class CRCandidate
    {
    public:
      /**
       *  @brief  Parametrised constructor
       */ 
      CRCandidate(const CRTaggingAlgorithm * const algorithm, const pandora::ParticleFlowObject * const pPfo, int id, int slice, double purity, double significance, bool isCosmicMuon );

      const CRTaggingAlgorithm * const          m_algorithm;    ///< The algorithm using this candidate

      // Data on the candidate PFO
      const pandora::ParticleFlowObject * const m_pPfo;         ///< Address of the candidate PFO
      pandora::CaloHitList                      m_hitList2D;    ///< List of all the 2D hits associated with the PFO
      pandora::CaloHitList                      m_hitList3D;    ///< List of all the 3D hits associated with the PFO
      int                                       m_id;           ///< Unique idendifier for the PFO

      int                                       m_slice;        ///< Slice ID
 
      // Variables used to identify cosmic rays
      int                                       m_n2DHits;      ///< Number of 2D hits over all views
      int                                       m_n3DHits;      ///< Number of 3D hits
    
      int                                       m_nDaughters;   ///< Number of daughter PFOs

      double                                    m_totalEnergy;  ///< Sum of all input hit energies
      double                                    m_meanEnergy;   ///< Mean of all input hit energies 

      bool                                      m_canFit;       ///< If there are a sufficient number of 3D hits to perform a fitting

      double                                    m_fitRMS;       ///< RMS of the 3D linear fit result 
      
      double                                    m_X1;           ///< x-coordinate of 1st fitted hit end point in 3D
      double                                    m_Y1;           ///< y-coordinate of 1st fitted hit end point in 3D
      double                                    m_Z1;           ///< z-coordinate of 1st fitted hit end point in 3D

      double                                    m_X2;           ///< x-coordinate of 2nd fitted hit end point in 3D
      double                                    m_Y2;           ///< y-coordinate of 2nd fitted hit end point in 3D
      double                                    m_Z2;           ///< z-coordinate of 2nd fitted hit end point in 3D

      double                                    m_length;       ///< Straight line length of the linear fit
      double                                    m_curvature;    ///< Measure of the curvature of the track
      /*
      double                                    m_theta1;       ///< Direction of the fit 1
      double                                    m_theta2;       ///< Direction of the fit 2
      */

      // Data used for classification
      double                                    m_purity;       ///< Neutrino purity of the PFO 
      double                                    m_significance; ///< Neutrino significance of the PFO 
      bool                                      m_isCosmicMuon; ///< If the main MCParticle associated with the PFO is a primary cosmic ray muon
      CRTagClass                                m_class;        ///< Classification of the PFO

      // Functions to calcluate the required variables for classification
      void CalculateNHits();                                    ///< Calculate m_nHits
      void CalculateNDaughters();                               ///< Calculate m_nDaughters
      void CalculateTotalEnergy();                              ///< Calculate m_totalEnergy
      void CalculateMeanEnergy();                               ///< Calculate m_meanEnergy 
      void CalculateFitVariables();                             ///< Calculate all variables which require the fit

      // Function to classify PFOs
      void DetermineClass();                                    ///< Determine m_class    

    };

    typedef std::list<CRCandidate> CRCandidateList;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Select a subset of true neutrinos representing those that should be used in performance metrics
     * 
     *  @param  pAllMCParticleList        address of the input mc particle list
     *  @param  selectedMCNeutrinoVector  to receive the populated selected true neutrino vector
     */
    void SelectTrueNeutrinos(const pandora::MCParticleList *const pAllMCParticleList, pandora::MCParticleVector &selectedMCNeutrinoVector) const;


    /**
     *  @brief  Select the MCParticles which are primary, visible and originating from a neutrino in selectedMCNeutrinoVector
     *
     *  @param  pAllMCParticleList             address of the input mc particle list
     *  @param  selectedMCNeutrinoVector       selected neutrinos to be used in performance metrics
     *  @param  candidateTargetMCParticleList  to receive the populated list of candidate target MCParticles
     */
    void GetCandidateTargetMCParticles( const pandora::MCParticleList *const pAllMCParticleList, pandora::MCParticleVector selectedMCNeutrinoVector, pandora::MCParticleList &candidateTargetMCParticleList) const;

    /**
     *  @brief  Relates all MCParticles back a candidate target, and associates their hits with the target if the MCParticle passes the checks for performance metrics (downstream of neutron)
     *
     *  @param  pCaloHitList                    list of all CaloHits in the event
     *  @param  pAllMCParticleList              list of all MCParticles in the event
     *  @param  candidateTargetMCParticleList   list of all candidate target MCParticles
     *  @param  candidateTargetToCaloHitListMap to reveive the populated list of hits associated with each candidate target MCParticle
     */
    void GetTargetToCaloHitMatches( const pandora::CaloHitList *const pCaloHitList, const pandora::MCParticleList *const pAllMCParticleList, pandora::MCParticleList candidateTargetMCParticleList, LArMonitoringHelper::MCContributionMap & candidateTargetToCaloHitListMap ) const;


    /**
     *  @brief  Get the main MCParticle associated with a given CaloHit
     *
     *  @param  pCaloHit input CaloHit
     *
     *  @return address of the mc particle associated to the input CaloHit (or nullptr if no such particle exists)
     */
    const pandora::MCParticle * GetMainMCParticle( const pandora::CaloHit * const pCaloHit ) const;

    /**
     *  @brief  Whether it is possible to navigate from a primary mc particle to a downstream mc particle without "passing through" a neutron
     * 
     *  @param  pOriginalPrimary the address of the original primary mc particle
     *  @param  pThisMCParticle the address of the current mc particle in the primary decay chain
     *  @param  pHitMCParticle the address of the mc particle associated to a calo hit
     * 
     *  @return boolean
     */
    bool PassMCParticleChecks(const pandora::MCParticle *const pOriginalPrimary, const pandora::MCParticle *const pThisMCParticle,
        const pandora::MCParticle *const pHitMCParticle) const;

    /**
     *  @brief  Selects target MCParticles from a list of candidate target MCParticles, by insisting a minimum number of hits
     *
     *  @param  candidateTargetToCaloHitListMap  input map between candidate target MCParticles and their associated hits
     *  @param  targetToCaloHitListMap           output map between chosen target MCParticles and their hits
     *  @param  caloHitToTargetMap               output map between hits and their assoaciated target MCParticle
     */
    void SelectTargetMCParticles( LArMonitoringHelper::MCContributionMap candidateTargetToCaloHitListMap, LArMonitoringHelper::MCContributionMap & targetToCaloHitListMap, LArMonitoringHelper::CaloHitToMCMap & caloHitToTargetMap) const;

    /**
     *  @brief  Get the origin of the supplied hits
     *
     *  @param  pCaloHitList        List of all CaloHits for which the origin should be found
     *  @param  caloHitToTargetMap  map between hits and their associated target MCParticle
     *  @param  caloHitToOriginMap  map between hits and their origin
     */
    void GetHitOrigins( const pandora::CaloHitList *const pCaloHitList, LArMonitoringHelper::CaloHitToMCMap caloHitToTargetMap, CaloHitToOriginMap & caloHitToOriginMap ) const;

    /**
     *  @brief  Get the purity of a list of PFOs
     *
     *  @param  pPfoList            input PfoList
     *  @param  caloHitToOriginMap  input map between hits and their origin
     *  @param  pfoToPurityMap      output map between pfos and their purity
     */ 
    void GetPurities( const pandora::PfoList * const pPfoList, CaloHitToOriginMap caloHitToOriginMap, PfoToDoubleMap & pfoToPurityMap ) const;

    /**
     *  @brief  Get the purity of a specific PFO
     *
     *  @param  pPfo                input PFO
     *  @param  caloHitToOriginMap  input map between hits and their origin
     *
     *  @return purity of the input PFO
     */
    double GetPurity( const pandora::ParticleFlowObject * const pPfo, CaloHitToOriginMap caloHitToOriginMap ) const;

    /**
     *  @brief  Get the significance of a list of PFOs
     *
     *  @param  pPfoList                input PfoList
     *  @param  caloHitToOriginMap      input map between hits and their origin
     *  @param  targetToCaloHitListMap  input map between target MCParticles and their associated hits
     *  @param  caloHitToTargetMap      input map between calohits and their associated target MCParticle
     *  @param  pfoToSignificanceMap    output map between pfos and their significance
     */ 
    void GetSignificances( const pandora::PfoList * const pPfoList, CaloHitToOriginMap caloHitToOriginMap, LArMonitoringHelper::MCContributionMap & targetToCaloHitListMap, LArMonitoringHelper::CaloHitToMCMap caloHitToTargetMap, PfoToDoubleMap & pfoToSignificanceMap ) const;

    /**
     *  @brief  Get the significance of a specific PFO
     *
     *  @param  pPfo                    input PFO
     *  @param  caloHitToOriginMap      input map between hits and their origin
     *  @param  targetToCaloHitListMap  input map between target MCParticles and their associated hits
     *  @param  caloHitToTargetMap      input map between calohits and their associated target MCParticle
     *
     *  @return significance of the input PFO
     */
    double GetSignificance( const pandora::ParticleFlowObject * const pPfo, CaloHitToOriginMap caloHitToOriginMap, LArMonitoringHelper::MCContributionMap & targetToCaloHitListMap, LArMonitoringHelper::CaloHitToMCMap caloHitToTargetMap ) const;


    /**
     *  @brief  Get mapping between PFOs and their IDs
     *
     *  @param  pPfoList  input list of PFOs
     *  @param  pfoIdMap  output mapping between PFOs and their ID
     */
    void GetPfoIds( const pandora::PfoList * const pPfoList, PfoToIntMap & pfoIdMap ) const;

    /**
     *  @brief  Get mapping between Target MCParticles and their IDs
     *
     *  @param  targetToCaloHitListMap  input mapping between target MCParticles and their hits
     *  @param  targetIdMap             output mapping between target MCParticles and their IDs
     */
    void GetTargetIds( LArMonitoringHelper::MCContributionMap targetToCaloHitListMap, MCToIntMap & targetIdMap ) const;


    /**
     *  @brief  Get the list of primary cosmic ray muons
     *
     *  @param  pMCParticleList input list of MCParticles
     *  @param  cosmicMuonList  output list of primary cosmic ray muons
     *
     */
    void GetPrimaryCosmicMuons( const pandora::MCParticleList *const pMCParticleList, pandora::MCParticleList & cosmicMuonList ) const;

    /**
     *  @brief  Determines if a pfo is a primary cosmic ray muon
     *
     *  @param  pPfoList             input list of PFOs to classify
     *  @param  caloHitToOriginMap   input mapping between calo hits and their origin
     *  @param  pfoToIsCosmicMuonMap output mapping between PFOs and a bool which indicated if they are a primary cosmic ray muon
     */
    void GetIsCosmicMuon( const pandora::PfoList * const pPfoList, CaloHitToOriginMap caloHitToOriginMap, PfoToBoolMap & pfoToIsCosmicMuonMap) const;

    /**
     *  @brief  Get mapping between PFOs that are associated with it other by pointing
     *
     *  @param  pPfoList            input list of PFOs
     *  @param  pfoAssociationsMap  output mapping between associated PFOs
     */
    void GetPfoAssociations( const pandora::PfoList * const pPfoList, PfoToPfoListMap & pfoAssociationMap ) const;

    /**
     *  @brief  Check if two PFO endpoints are associated by distance of closest approach
     *
     *  @param  endPoint1  position vector of an endpoint of PFO 1
     *  @param  endDir1    direction vector of an endpoint of PFO 1
     *  @param  endPoint2  position vector of an endpoint of PFO 2
     *  @param  endDir2    direction vector of an endpoint of PFOs
     *
     *  @return If the PFOs are assoicated
     */
    bool CheckAssociation( pandora::CartesianVector endPoint1, pandora::CartesianVector endDir1, pandora::CartesianVector endPoint2, pandora::CartesianVector endDir2 ) const;

    /**
     *  @brief  Break the event up into slices of associated PFOs
     *
     *  @param  pfoAssociationMap  mapping between PFOs and other associated PFOs
     *  @param  pfoToSliceIdMap    mapping between PFOs and thier slice ID
     */   
    void SliceEvent( PfoToPfoListMap pfoAssociationMap, PfoToIntMap & pfoToSliceIdMap ) const;


    /**
     *  @brief  Fill a slice iteratively using PFO associations
     *
     *  @param  pfoAssociationMap  mapping between PFOs and other associated PFOs
     *  @param  pPfo               PFO to add to the slice
     *  @param  slice              the slice to add PFOs to
     */
    void FillSlice( PfoToPfoListMap pfoAssociationMap, const pandora::ParticleFlowObject * const pPfo, pandora::PfoList & slice ) const;

    /**
     *  @brief  Make a list of CRCandidates
     *
     *  @param  pPfoList              input list of PFOs
     *  @param  pfoToPurityMap        input mapping between PFOs and their purity
     *  @param  pfoToSignificanceMap  input mapping between PFOs and their significance
     *  @param  pfoIdMap              input mapping between PFOs and their ID
     *  @param  pfoToSliceIdMap       input mapping between PFOs and their slice id
     *  @param  candidates            output list of CRCandidates 
     */
    void GetCRCandidates( const pandora::PfoList * const pPfoList, PfoToDoubleMap pfoToPurityMap, PfoToDoubleMap pfoToSignificanceMap, PfoToIntMap pfoIdMap, PfoToIntMap pfoToSliceIdMap, PfoToBoolMap pfoToIsCosmicMuonMap, CRCandidateList & candidates ) const;
  
    /**
     *  @brief  Make a list of CRCandidates
     *
     *  @param  pPfoList           input list of PFOs
     *  @param  caloHit2DToPfoMap  output mapping between 2D CaloHits and their associated PFO
     *  @param  caloHit3DToPfoMap  output mapping between 3D CaloHits and their associated PFO
     */
    void GetCaloHitToPfoMap( const pandora::PfoList * const pPfoList, CaloHitToPfoMap & caloHit2DToPfoMap, CaloHitToPfoMap & caloHit3DToPfoMap ) const;

    /**
     *  @brief  Write information on all 2D hits to a root file
     *
     *  @param  caloHitsToOriginMap  mapping between CaloHits and their origin
     *  @param  caloHitToTargetMap   mapping between CaloHits and the associated target MCParticle (if any)
     *  @param  caloHitToPfoMap      mapping between CaloHits and their associated PFO
     *  @param  pfoIdMap             mapping between PFOs and their IDs
     *  @param  targetIdMap          mapping between target MCParticles and their IDs
     */
    void Write2DHits( CaloHitToOriginMap caloHitToOriginMap, LArMonitoringHelper::CaloHitToMCMap caloHitToTargetMap, CaloHitToPfoMap caloHitToPfoMap, PfoToIntMap pfoIdMap, MCToIntMap targetIdMap ) const;

    /**
     *  @brief  Write information on all 3D hits to a root file
     *
     *  @param  caloHitToPfoMap      mapping between CaloHits and their associated PFO
     *  @param  pfoIdMap             mapping between PFOs and their IDs
     */
    void Write3DHits( CaloHitToPfoMap caloHitToPfoMap, PfoToIntMap pfoIdMap ) const;

    /**
     *  @brief  Write information on all target MCParticles
     *
     *  @param  targetIdMap  mapping between target MCParticles and their IDs
     */
    void WriteTargets( MCToIntMap targetIdMap ) const;

    /**
     *  @brief  Write information on all primary muon cosmic MCParticles
     *
     *  @param  cosmicMuonList  list of primary cosmic muon MCParticles
     */
    void WriteCosmics( pandora::MCParticleList cosmicMuonList ) const;

    /**
     *  @brief  Write information on all PFOs
     *
     *  @param  candidates  list of all cosmic ray candidates 
     */
    void WritePfos( CRCandidateList candidates ) const;


    std::string             m_caloHitListName;              ///< The name of the input CaloHit list
    std::string             m_inputPfoListName;             ///< The name of the input PFO list
    std::string             m_mcParticleListName;           ///< The name of the input MCParticle list

    int                     m_eventNumber;
    int                     m_fileId;

    float                   m_maxPhotonPropagation;         ///< Maximum distance travelled by photon, downstream of a track, in mc particle hierarchy

    int                     m_matchingMinPrimaryHits;       ///< The minimum number of good mc primary hits required
    int                     m_matchingMinHitsForGoodView;   ///< The minimum number of good mc primary hits in given view to declare view to be good
    int                     m_matchingMinPrimaryGoodViews;  ///< The minimum number of good views for a mc primary
    
    double                  m_angularUncertainty;           ///< The uncertainty in degrees for the angle of a PFO
    double                  m_positionalUncertainty;        ///< The uncertainty in cm for the position of PFO endpoint in 3D
    double                  m_maxAssociationDist;           ///< The maximum distance from endpoint to point of closest approach
    
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CRTaggingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CRTaggingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CR_TAGGING_ALGORITHM_H
