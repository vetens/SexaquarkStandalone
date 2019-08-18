#include "../interface/AnalyzerRECO.h"
#include <typeinfo>

AnalyzerRECO::AnalyzerRECO(edm::ParameterSet const& pset):
  m_lookAtAntiS(pset.getUntrackedParameter<bool>("lookAtAntiS")),
  m_bsTag(pset.getParameter<edm::InputTag>("beamspot")),
  m_offlinePVTag(pset.getParameter<edm::InputTag>("offlinePV")),
  m_generalTracksTag(pset.getParameter<edm::InputTag>("generalTracksCollection")),
  m_sCandsTag(pset.getParameter<edm::InputTag>("sexaqCandidates")),
  m_V0KsTag(pset.getParameter<edm::InputTag>("V0KsCollection")),
  m_V0LTag(pset.getParameter<edm::InputTag>("V0LCollection")),

  m_bsToken    (consumes<reco::BeamSpot>(m_bsTag)),
  m_offlinePVToken    (consumes<vector<reco::Vertex>>(m_offlinePVTag)),
  m_generalTracksToken(consumes<View<reco::Track> >(m_generalTracksTag)),
  m_sCandsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_sCandsTag)),
  m_V0KsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_V0KsTag)),
  m_V0LToken(consumes<vector<reco::VertexCompositeCandidate> >(m_V0LTag))
  


{

}


void AnalyzerRECO::beginJob() {

     TFileDirectory dir_RECO = m_fs->mkdir("RECO"); 

     TFileDirectory dir_PV = m_fs->mkdir("PV");
     histos_th1f["h_PV_n"] = dir_PV.make<TH1F>(b+"h_PV_n", "; #Primary Vertices; #entries ",100,0,100);
     histos_th1f["h_PV_vx"] = dir_PV.make<TH1F>(b+"h_PV_vx", "; vx (cm); #entries ",200,-10,10);
     histos_th1f["h_PV_vy"] = dir_PV.make<TH1F>(b+"h_PV_vy", "; vy (cm); #entries ",200,-10,10);
     histos_th1f["h_PV_vz"] = dir_PV.make<TH1F>(b+"h_PV_vz", "; vz (cm); #entries ",600,-300,300);

     TFileDirectory dir_RECO_Ks = dir_RECO.mkdir("RECO_Ks");
     histos_th1f["h_RECO_Ks_pt"] = dir_RECO_Ks.make<TH1F>(b+"h_RECO_Ks_pt", "; Ks pT (GeV); #entries ",200,0,20);
     histos_th1f["h_RECO_Ks_eta"] = dir_RECO_Ks.make<TH1F>(b+"h_RECO_Ks_eta", "; Ks #eta; #entries ",200,-10,10);
     histos_th1f["h_RECO_Ks_phi"] = dir_RECO_Ks.make<TH1F>(b+"h_RECO_Ks_phi", "; Ks #phi (rad); #entries ",200,-4,4);
     histos_th2f["h2_RECO_Ks_vx_vy"] = dir_RECO_Ks.make<TH2F>(b+"h2_RECO_Ks_vx_vy", "; Ks vx (cm); Ks vy (cm)",400,-200,200,400,-200,200);
     histos_th2f["h2_RECO_Ks_vx_vz"] = dir_RECO_Ks.make<TH2F>(b+"h2_RECO_Ks_vx_vz", "; Ks vx (cm); Ks vz (cm)",400,-200,200,800,-400,400);
     histos_th1f["h_RECO_Ks_lxy"] = dir_RECO_Ks.make<TH1F>(b+"h_RECO_Ks_lxy", "; Ks lxy(beamspot, Ks vertex) (cm); #entries ",200,0,200);
     histos_th1f["h_RECO_Ks_vz"] = dir_RECO_Ks.make<TH1F>(b+"h_RECO_Ks_vz", "; Ks vz(Ks vertex) (cm); #entries ",600,-300,300);
     histos_th1f["h_RECO_Ks_m"] = dir_RECO_Ks.make<TH1F>(b+"h_RECO_Ks_m", "; Ks m (GeV); #entries ",4000,0,4);

     TFileDirectory dir_RECO_Ks_extra_mass_cut = dir_RECO_Ks.mkdir("RECO_Ks_extra_mass_cut");
     histos_th2f["h2_RECO_Ks_lxy_vz"] = dir_RECO_Ks_extra_mass_cut.make<TH2F>(b+"h2_RECO_Ks_lxy_vz", "; Ks |vz| (cm); Ks l_{0} (cm)",300,0,300,120,0,120);
     histos_th2f["h2_RECO_Ks_pt_pz"] = dir_RECO_Ks_extra_mass_cut.make<TH2F>(b+"h2_RECO_Ks_pt_pz", "; Ks pt (GeV); Ks |pz| (GeV)",100,0,10,100,0,40);
     histos_th2f["h2_RECO_Ks_dxy_dz"] = dir_RECO_Ks_extra_mass_cut.make<TH2F>(b+"h2_RECO_Ks_dxy_dz", "; Ks dxy (cm); Ks dz (cm)",100,-20,20,100,-40,40);

     histos_th2f["h2_RECO_Ks_lxy_dz"] = dir_RECO_Ks_extra_mass_cut.make<TH2F>(b+"h2_RECO_Ks_lxy_dz", "; Ks |dz| (cm); Ks l_{0} (cm)",300,0,300,120,0,120);
     histos_th1f["h_RECO_Ks_dxy"] = dir_RECO_Ks_extra_mass_cut.make<TH1F>(b+"h_RECO_Ks_dxy", "; Ks dxy(beamspot) (cm); #entries ",400,-20,20);
     histos_th1f["h_RECO_Ks_dz"] = dir_RECO_Ks_extra_mass_cut.make<TH1F>(b+"h_RECO_Ks_dz", "; Ks dz(beamspot) (cm); #entries ",400,-200,200);
     histos_th1f["h_RECO_Ks_dz_min"] = dir_RECO_Ks_extra_mass_cut.make<TH1F>(b+"h_RECO_Ks_dz_min", "; Ks dz(beamspot) (cm); #entries ",400,-200,200);
     histos_th2f["h2_RECO_Ks_mass_dz"] = dir_RECO_Ks_extra_mass_cut.make<TH2F>(b+"h2_RECO_Ks_mass_dz", "; Ks |d_{z}| (cm); Ks mass (GeV)",300,0,300,100,0,1);

     TFileDirectory dir_RECO_Ks_extra_mass_cut_kinregA = dir_RECO_Ks_extra_mass_cut.mkdir("RECO_Ks_extra_mass_cut_kinregA");
     histos_th2f["h2_RECO_Ks_lxy_vz_kinregA"] = dir_RECO_Ks_extra_mass_cut_kinregA.make<TH2F>(b+"h2_RECO_Ks_lxy_vz_kinregA", "; Ks |vz| (cm); Ks l_{0} (cm)",300,0,300,120,0,120);
     histos_th2f["h2_RECO_Ks_pt_pz_kinregA"] = dir_RECO_Ks_extra_mass_cut_kinregA.make<TH2F>(b+"h2_RECO_Ks_pt_pz_kinregA", "; Ks pt (GeV); Ks |pz| (GeV)",100,0,10,100,0,40);
     histos_th2f["h2_RECO_Ks_dxy_dz_kinregA"] = dir_RECO_Ks_extra_mass_cut_kinregA.make<TH2F>(b+"h2_RECO_Ks_dxy_dz_kinregA", "; Ks dxy (cm); Ks dz (cm)",100,-2,2,100,-40,40);
     histos_th2f["h2_RECO_Ks_eta_phi_kinregA"] = dir_RECO_Ks_extra_mass_cut_kinregA.make<TH2F>(b+"h2_RECO_Ks_eta_phi_kinregA", "; Ks #eta (cm); Ks #phi (cm)",100,-4,4,100,-4,4);

     TFileDirectory dir_RECO_Lambda = dir_RECO.mkdir("RECO_Lambda");
     histos_th1f["h_RECO_Lambda_pt"] = dir_RECO_Lambda.make<TH1F>(b+"h_RECO_Lambda_pt", "; #bar{#Lambda} or #Lambda pT (GeV); #entries ",200,0,20);
     histos_th1f["h_RECO_Lambda_eta"] = dir_RECO_Lambda.make<TH1F>(b+"h_RECO_Lambda_eta", "; #bar{#Lambda} or #Lambda #eta; #entries ",200,-10,10);
     histos_th1f["h_RECO_Lambda_phi"] = dir_RECO_Lambda.make<TH1F>(b+"h_RECO_Lambda_phi", "; #bar{#Lambda} or #Lambda #phi (rad); #entries ",200,-4,4);
     histos_th2f["h2_RECO_Lambda_vx_vy"] = dir_RECO_Lambda.make<TH2F>(b+"h2_RECO_Lambda_vx_vy", "; #bar{#Lambda} or #Lambda vx (cm); #bar{#Lambda} vy (cm)",400,-200,200,400,-200,200);
     histos_th2f["h2_RECO_Lambda_vx_vz"] = dir_RECO_Lambda.make<TH2F>(b+"h2_RECO_Lambda_vx_vz", "; #bar{#Lambda} or #Lambda vx (cm); #bar{#Lambda} vz (cm)",400,-200,200,800,-400,400);
     histos_th1f["h_RECO_Lambda_lxy"] = dir_RECO_Lambda.make<TH1F>(b+"h_RECO_Lambda_lxy", "; #bar{#Lambda} or #Lambda lxy(beamspot, #bar{#Lambda} vertex) (cm); #entries ",200,0,200);
     histos_th1f["h_RECO_Lambda_vz"] = dir_RECO_Lambda.make<TH1F>(b+"h_RECO_Lambda_vz", "; #bar{#Lambda} or #Lambda vz(Ks vertex) (cm); #entries ",600,-300,300);
     histos_th1f["h_RECO_Lambda_m"] = dir_RECO_Lambda.make<TH1F>(b+"h_RECO_Lambda_m", "; #bar{#Lambda} or #Lambda m (GeV); #entries ",4000,0,4);
     histos_th1f["h_RECO_Lambda_deltaRMin_GEN_Lambda"] = dir_RECO_Lambda.make<TH1F>(b+"h_RECO_Lambda_deltaRMin_GEN_Lambda", "; #DeltaR_{min} (#bar{#Lambda} or or #Lambda GEN, #bar{#Lambda} or #Lambda RECO); #entries ",1000,0,10);

     TFileDirectory dir_RECO_Lambda_extra_mass_cut = dir_RECO_Lambda.mkdir("RECO_Lambda_extra_mass_cut");
     histos_th2f["h2_RECO_Lambda_lxy_vz"] = dir_RECO_Lambda_extra_mass_cut.make<TH2F>(b+"h2_RECO_Lambda_lxy_vz", "; #bar{#Lambda} or #Lambda |v_{z}| (cm); #bar{#Lambda} l_{0} (cm)",300,0,300,120,0,120);
     histos_th2f["h2_RECO_Lambda_lxy_dz"] = dir_RECO_Lambda_extra_mass_cut.make<TH2F>(b+"h2_RECO_Lambda_lxy_dz", "; #bar{#Lambda} or #Lambda |d_{z}| (cm); #bar{#Lambda} l_{0} (cm)",300,0,300,120,0,120);
     histos_th1f["h_RECO_Lambda_dxy"] = dir_RECO_Lambda_extra_mass_cut.make<TH1F>(b+"h_RECO_Lambda_dxy", "; #bar{#Lambda} or #Lambda dxy(beamspot) (cm); #entries ",400,-20,20);
     histos_th1f["h_RECO_Lambda_dz"] = dir_RECO_Lambda_extra_mass_cut.make<TH1F>(b+"h_RECO_Lambda_dz", "; #bar{#Lambda} or #Lambda dz(beamspot) (cm); #entries ",400,-200,200);
     histos_th2f["h2_RECO_Lambda_mass_dz"] = dir_RECO_Lambda_extra_mass_cut.make<TH2F>(b+"h2_RECO_Lambda_mass_dz", "; #bar{#Lambda} or #Lambda |d_{z}| (cm); #bar{#Lambda} mass (GeV)",300,0,300,100,1,2);

     TFileDirectory dir_RECO_Lambda_extra_mass_cut_kinregA = dir_RECO_Lambda_extra_mass_cut.mkdir("RECO_Lambda_extra_mass_cut_kinregA");
     histos_th2f["h2_RECO_Lambda_lxy_vz_kinregA"] = dir_RECO_Lambda_extra_mass_cut_kinregA.make<TH2F>(b+"h2_RECO_Lambda_lxy_vz_kinregA", "; #bar{#Lambda} or #Lambda |vz| (cm); #bar{#Lambda} or #Lambda l_{0} (cm)",300,0,300,120,0,120);
     histos_th2f["h2_RECO_Lambda_pt_pz_kinregA"] = dir_RECO_Lambda_extra_mass_cut_kinregA.make<TH2F>(b+"h2_RECO_Lambda_pt_pz_kinregA", "; #bar{#Lambda} or #Lambda pt (GeV); #bar{#Lambda} or #Lambda |pz| (GeV)",100,0,10,100,0,40);
     histos_th2f["h2_RECO_Lambda_dxy_dz_kinregA"] = dir_RECO_Lambda_extra_mass_cut_kinregA.make<TH2F>(b+"h2_RECO_Lambda_dxy_dz_kinregA", "; #bar{#Lambda} or #Lambda dxy (cm); #bar{#Lambda} or #Lambda dz (cm)",100,-2,2,100,-40,40);
     histos_th2f["h2_RECO_Lambda_eta_phi_kinregA"] = dir_RECO_Lambda_extra_mass_cut_kinregA.make<TH2F>(b+"h2_RECO_Lambda_eta_phi_kinregA", "; #bar{#Lambda} or #Lambda #eta; #bar{#Lambda} or #Lambda #phi)",100,-4,4,100,-4,4);


     //Just in general for tracks:
     TFileDirectory dir_RECO_Tracks = dir_RECO.mkdir("RECO_Tracks");
     histos_th1f["h_RECO_Track_vz"] = dir_RECO_Tracks.make<TH1F>(b+"h_RECO_Track_vz", "; Track |v_{z}| (cm); #entries ",300,0,300);
     histos_th1f["h_RECO_Track_dz"] = dir_RECO_Tracks.make<TH1F>(b+"h_RECO_Track_dz", "; Track |d_{z}| (cm); #entries ",300,0,300);
     histos_th1f["h_RECO_Track_dsz"] = dir_RECO_Tracks.make<TH1F>(b+"h_RECO_Track_dsz", "; Track |ds_{z}| (cm); #entries ",300,0,300);

     TFileDirectory dir_RECO_AntiS = dir_RECO.mkdir("RECO_AntiS");
     TFileDirectory dir_RECO_AntiS_PRE_BACKGROUND_CUTS_Kinematics = dir_RECO_AntiS.mkdir("RECO_AntiS_PRE_BACKGROUND_CUTS_Kinematics");
     histos_th1f["h_RECO_AntiS_pt"] = dir_RECO_AntiS_PRE_BACKGROUND_CUTS_Kinematics.make<TH1F>(b+"h_RECO_AntiS_pt", "; #bar{S} pT (GeV); #entries ",200,0,20);
     histos_th1f["h_RECO_AntiS_phi"] = dir_RECO_AntiS_PRE_BACKGROUND_CUTS_Kinematics.make<TH1F>(b+"h_RECO_AntiS_phi", "; #bar{S} #phi (rad); #entries ",200,-4,4);
     histos_th2f["h2_RECO_AntiS_vx_vy"] = dir_RECO_AntiS_PRE_BACKGROUND_CUTS_Kinematics.make<TH2F>(b+"h2_RECO_AntiS_vx_vy", "; #bar{S} vx (cm); #bar{S} vy (cm)",400,-200,200,400,-200,200);
     histos_th2f["h2_RECO_AntiS_vx_vz"] = dir_RECO_AntiS_PRE_BACKGROUND_CUTS_Kinematics.make<TH2F>(b+"h2_RECO_AntiS_vx_vz", "; #bar{S} vx (cm); #bar{S} vz (cm)",400,-200,200,800,-400,400);
     histos_th1f["h_RECO_AntiS_dxy"] = dir_RECO_AntiS_PRE_BACKGROUND_CUTS_Kinematics.make<TH1F>(b+"h_RECO_AntiS_dxy", "; #bar{S} dxy(beamspot) (cm); #entries ",400,-20,20);
     histos_th1f["h_RECO_AntiS_m"] = dir_RECO_AntiS_PRE_BACKGROUND_CUTS_Kinematics.make<TH1F>(b+"h_RECO_AntiS_m", "; #bar{S} m (GeV); #entries ",200,0,20);
     histos_th1f["h_RECO_AntiS_deltaRDaughters"] = dir_RECO_AntiS_PRE_BACKGROUND_CUTS_Kinematics.make<TH1F>(b+"h_RECO_AntiS_deltaRDaughters", "; #bar{S} #DeltaR (Ks,#bar{#Lambda}); #entries ",100,0,10);
    
     TFileDirectory dir_RECO_AntiS_cut_variables_distributions = dir_RECO_AntiS.mkdir("RECO_AntiS_cut_variables_distributions");
     histos_th1f["h_RECO_AntiS_lxy"] = dir_RECO_AntiS_cut_variables_distributions.make<TH1F>(b+"h_RECO_AntiS_lxy", "; #bar{S} lxy(beamspot, #bar{S} vertex) (cm); #entries ",200,0,200);
     histos_th1f["h_RECO_AntiS_error_lxy"] = dir_RECO_AntiS_cut_variables_distributions.make<TH1F>(b+"h_RECO_AntiS_error_lxy", "; #bar{S} error lxy(beamspot, #bar{S} vertex) (cm); #entries ",200,0,10);
     histos_th1f["h_RECO_AntiS_mass"] = dir_RECO_AntiS_cut_variables_distributions.make<TH1F>(b+"h_RECO_AntiS_mass", "; S or #bar{S} mass (GeV); #entries ",200,0,10);
     histos_th1f["h_RECO_AntiS_vertex_chi2_ndof"] = dir_RECO_AntiS_cut_variables_distributions.make<TH1F>(b+"h_RECO_AntiS_vertex_chi2_ndof", "; S or #bar{S} vertex #chi{^2}/ndof; #entries ",200,0,20);
     histos_th1f["h_RECO_AntiS_deltaPhiDaughters"] = dir_RECO_AntiS_cut_variables_distributions.make<TH1F>(b+"h_RECO_AntiS_deltaPhiDaughters", "; #bar{S} #Delta#Phi (Ks,#bar{#Lambda}); #entries ",100,-5,5);
     histos_th1f["h_RECO_AntiS_openingsAngleDaughters"] = dir_RECO_AntiS_cut_variables_distributions.make<TH1F>(b+"h_RECO_AntiS_openingsAngleDaughters", ";S or #bar{S} openings angle (Ks,#Lambda or #bar{#Lambda}); #entries ",80,-4,4);
     histos_th1f["h_RECO_AntiS_deltaEtaDaughters"] = dir_RECO_AntiS_cut_variables_distributions.make<TH1F>(b+"h_RECO_AntiS_deltaEtaDaughters", ";S or #bar{S} #Delta#Eta (Ks,#Lambda or #bar{#Lambda}); #entries ",80,-4,4);
     histos_th1f["h_RECO_AntiS_dxy_over_lxy"] = dir_RECO_AntiS_cut_variables_distributions.make<TH1F>(b+"h_RECO_AntiS_dxy_over_lxy", ";dxy S or #bar{S}/lxy interaction vertex; #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_dxy_over_lxy_Ks"] = dir_RECO_AntiS_cut_variables_distributions.make<TH1F>(b+"h_RECO_AntiS_dxy_over_lxy_Ks", ";dxy Ks/lxy interaction vertex; #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_dxy_over_lxy_AntiL"] = dir_RECO_AntiS_cut_variables_distributions.make<TH1F>(b+"h_RECO_AntiS_dxy_over_lxy_AntiL", ";dxy #Lambda or #bar{#Lambda}/lxy interaction vertex; #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_pt_AntiL"] = dir_RECO_AntiS_cut_variables_distributions.make<TH1F>(b+"h_RECO_AntiS_pt_AntiL", ";pt #Lambda or #bar{#Lambda} (GeV); #entries ",100,0,10);
     histos_th1f["h_RECO_AntiS_dzPVmin"] = dir_RECO_AntiS_cut_variables_distributions.make<TH1F>(b+"h_RECO_AntiS_dzPVmin", "; S or #bar{S} min dz(PV); #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_vz"] = dir_RECO_AntiS_cut_variables_distributions.make<TH1F>(b+"h_RECO_AntiS_vz", "; vz S or #bar{S} interaction vertex; #entries ",2000,-100,100);
     histos_th1f["h_RECO_AntiS_eta"] = dir_RECO_AntiS_cut_variables_distributions.make<TH1F>(b+"h_RECO_AntiS_eta", "; #bar{S} #eta; #entries ",200,-10,10);

     TFileDirectory dir_RECO_AntiS_cut_variables_distributions_L1 = dir_RECO_AntiS.mkdir("RECO_AntiS_cut_variables_distributions_L1");
     histos_th1f["h_RECO_AntiS_lxy_L1"] = dir_RECO_AntiS_cut_variables_distributions_L1.make<TH1F>(b+"h_RECO_AntiS_lxy_L1", "; #bar{S} lxy(beamspot, #bar{S} vertex) (cm); #entries ",200,0,200);
     histos_th1f["h_RECO_AntiS_error_lxy_L1"] = dir_RECO_AntiS_cut_variables_distributions_L1.make<TH1F>(b+"h_RECO_AntiS_error_lxy_L1", "; #bar{S} error lxy(beamspot, #bar{S} vertex) (cm); #entries ",200,0,10);
     histos_th1f["h_RECO_AntiS_mass_L1"] = dir_RECO_AntiS_cut_variables_distributions_L1.make<TH1F>(b+"h_RECO_AntiS_mass_L1", "; S or #bar{S} mass (GeV); #entries ",200,0,10);
     histos_th1f["h_RECO_AntiS_vertex_chi2_ndof_L1"] = dir_RECO_AntiS_cut_variables_distributions_L1.make<TH1F>(b+"h_RECO_AntiS_vertex_chi2_ndof_L1", "; S or #bar{S} vertex #chi{^2}/ndof; #entries ",200,0,20);
     histos_th1f["h_RECO_AntiS_deltaPhiDaughters_L1"] = dir_RECO_AntiS_cut_variables_distributions_L1.make<TH1F>(b+"h_RECO_AntiS_deltaPhiDaughters_L1", "; #bar{S} #Delta#Phi (Ks,#bar{#Lambda}); #entries ",100,-5,5);
     histos_th1f["h_RECO_AntiS_openingsAngleDaughters_L1"] = dir_RECO_AntiS_cut_variables_distributions_L1.make<TH1F>(b+"h_RECO_AntiS_openingsAngleDaughters_L1", ";S or #bar{S} openings angle (Ks,#Lambda or #bar{#Lambda}); #entries ",80,-4,4);
     histos_th1f["h_RECO_AntiS_deltaEtaDaughters_L1"] = dir_RECO_AntiS_cut_variables_distributions_L1.make<TH1F>(b+"h_RECO_AntiS_deltaEtaDaughters_L1", ";S or #bar{S} #Delta#Eta (Ks,#Lambda or #bar{#Lambda}); #entries ",80,-4,4);
     histos_th1f["h_RECO_AntiS_dxy_over_lxy_L1"] = dir_RECO_AntiS_cut_variables_distributions_L1.make<TH1F>(b+"h_RECO_AntiS_dxy_over_lxy_L1", ";dxy S or #bar{S}/lxy interaction vertex; #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_dxy_over_lxy_Ks_L1"] = dir_RECO_AntiS_cut_variables_distributions_L1.make<TH1F>(b+"h_RECO_AntiS_dxy_over_lxy_Ks_L1", ";dxy Ks/lxy interaction vertex; #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_dxy_over_lxy_AntiL_L1"] = dir_RECO_AntiS_cut_variables_distributions_L1.make<TH1F>(b+"h_RECO_AntiS_dxy_over_lxy_AntiL_L1", ";dxy #Lambda or #bar{#Lambda}/lxy interaction vertex; #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_pt_AntiL_L1"] = dir_RECO_AntiS_cut_variables_distributions_L1.make<TH1F>(b+"h_RECO_AntiS_pt_AntiL_L1", ";pt #Lambda or #bar{#Lambda} (GeV); #entries ",100,0,10);
     histos_th1f["h_RECO_AntiS_dzPVmin_L1"] = dir_RECO_AntiS_cut_variables_distributions_L1.make<TH1F>(b+"h_RECO_AntiS_dzPVmin_L1", "; S or #bar{S} min dz(PV); #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_vz_L1"] = dir_RECO_AntiS_cut_variables_distributions_L1.make<TH1F>(b+"h_RECO_AntiS_vz_L1", "; vz S or #bar{S} interaction vertex; #entries ",2000,-100,100);
     histos_th1f["h_RECO_AntiS_eta_L1"] = dir_RECO_AntiS_cut_variables_distributions_L1.make<TH1F>(b+"h_RECO_AntiS_eta_L1", "; #bar{S} #eta; #entries ",200,-10,10);

     TFileDirectory dir_RECO_AntiS_cut_variables_distributions_L2 = dir_RECO_AntiS.mkdir("RECO_AntiS_cut_variables_distributions_L2");
     histos_th1f["h_RECO_AntiS_lxy_L2"] = dir_RECO_AntiS_cut_variables_distributions_L2.make<TH1F>(b+"h_RECO_AntiS_lxy_L2", "; #bar{S} lxy(beamspot, #bar{S} vertex) (cm); #entries ",200,0,200);
     histos_th1f["h_RECO_AntiS_error_lxy_L2"] = dir_RECO_AntiS_cut_variables_distributions_L2.make<TH1F>(b+"h_RECO_AntiS_error_lxy_L2", "; #bar{S} error lxy(beamspot, #bar{S} vertex) (cm); #entries ",200,0,10);
     histos_th1f["h_RECO_AntiS_mass_L2"] = dir_RECO_AntiS_cut_variables_distributions_L2.make<TH1F>(b+"h_RECO_AntiS_mass_L2", "; S or #bar{S} mass (GeV); #entries ",200,0,10);
     histos_th1f["h_RECO_AntiS_vertex_chi2_ndof_L2"] = dir_RECO_AntiS_cut_variables_distributions_L2.make<TH1F>(b+"h_RECO_AntiS_vertex_chi2_ndof_L2", "; S or #bar{S} vertex #chi{^2}/ndof; #entries ",200,0,20);
     histos_th1f["h_RECO_AntiS_deltaPhiDaughters_L2"] = dir_RECO_AntiS_cut_variables_distributions_L2.make<TH1F>(b+"h_RECO_AntiS_deltaPhiDaughters_L2", "; #bar{S} #Delta#Phi (Ks,#bar{#Lambda}); #entries ",100,-5,5);
     histos_th1f["h_RECO_AntiS_openingsAngleDaughters_L2"] = dir_RECO_AntiS_cut_variables_distributions_L2.make<TH1F>(b+"h_RECO_AntiS_openingsAngleDaughters_L2", ";S or #bar{S} openings angle (Ks,#Lambda or #bar{#Lambda}); #entries ",80,-4,4);
     histos_th1f["h_RECO_AntiS_deltaEtaDaughters_L2"] = dir_RECO_AntiS_cut_variables_distributions_L2.make<TH1F>(b+"h_RECO_AntiS_deltaEtaDaughters_L2", ";S or #bar{S} #Delta#Eta (Ks,#Lambda or #bar{#Lambda}); #entries ",80,-4,4);
     histos_th1f["h_RECO_AntiS_dxy_over_lxy_L2"] = dir_RECO_AntiS_cut_variables_distributions_L2.make<TH1F>(b+"h_RECO_AntiS_dxy_over_lxy_L2", ";dxy S or #bar{S}/lxy interaction vertex; #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_dxy_over_lxy_Ks_L2"] = dir_RECO_AntiS_cut_variables_distributions_L2.make<TH1F>(b+"h_RECO_AntiS_dxy_over_lxy_Ks_L2", ";dxy Ks/lxy interaction vertex; #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_dxy_over_lxy_AntiL_L2"] = dir_RECO_AntiS_cut_variables_distributions_L2.make<TH1F>(b+"h_RECO_AntiS_dxy_over_lxy_AntiL_L2", ";dxy #Lambda or #bar{#Lambda}/lxy interaction vertex; #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_pt_AntiL_L2"] = dir_RECO_AntiS_cut_variables_distributions_L2.make<TH1F>(b+"h_RECO_AntiS_pt_AntiL_L2", ";pt #Lambda or #bar{#Lambda} (GeV); #entries ",100,0,10);
     histos_th1f["h_RECO_AntiS_dzPVmin_L2"] = dir_RECO_AntiS_cut_variables_distributions_L2.make<TH1F>(b+"h_RECO_AntiS_dzPVmin_L2", "; S or #bar{S} min dz(PV); #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_vz_L2"] = dir_RECO_AntiS_cut_variables_distributions_L2.make<TH1F>(b+"h_RECO_AntiS_vz_L2", "; vz S or #bar{S} interaction vertex; #entries ",2000,-100,100);
     histos_th1f["h_RECO_AntiS_eta_L2"] = dir_RECO_AntiS_cut_variables_distributions_L2.make<TH1F>(b+"h_RECO_AntiS_eta_L2", "; #bar{S} #eta; #entries ",200,-10,10);

     TFileDirectory dir_RECO_AntiS_cut_variables_distributions_L3 = dir_RECO_AntiS.mkdir("RECO_AntiS_cut_variables_distributions_L3");
     histos_th1f["h_RECO_AntiS_lxy_L3"] = dir_RECO_AntiS_cut_variables_distributions_L3.make<TH1F>(b+"h_RECO_AntiS_lxy_L3", "; #bar{S} lxy(beamspot, #bar{S} vertex) (cm); #entries ",200,0,200);
     histos_th1f["h_RECO_AntiS_error_lxy_L3"] = dir_RECO_AntiS_cut_variables_distributions_L3.make<TH1F>(b+"h_RECO_AntiS_error_lxy_L3", "; #bar{S} error lxy(beamspot, #bar{S} vertex) (cm); #entries ",200,0,10);
     histos_th1f["h_RECO_AntiS_mass_L3"] = dir_RECO_AntiS_cut_variables_distributions_L3.make<TH1F>(b+"h_RECO_AntiS_mass_L3", "; S or #bar{S} mass (GeV); #entries ",200,0,10);
     histos_th1f["h_RECO_AntiS_vertex_chi2_ndof_L3"] = dir_RECO_AntiS_cut_variables_distributions_L3.make<TH1F>(b+"h_RECO_AntiS_vertex_chi2_ndof_L3", "; S or #bar{S} vertex #chi{^2}/ndof; #entries ",200,0,20);
     histos_th1f["h_RECO_AntiS_deltaPhiDaughters_L3"] = dir_RECO_AntiS_cut_variables_distributions_L3.make<TH1F>(b+"h_RECO_AntiS_deltaPhiDaughters_L3", "; #bar{S} #Delta#Phi (Ks,#bar{#Lambda}); #entries ",100,-5,5);
     histos_th1f["h_RECO_AntiS_openingsAngleDaughters_L3"] = dir_RECO_AntiS_cut_variables_distributions_L3.make<TH1F>(b+"h_RECO_AntiS_openingsAngleDaughters_L3", ";S or #bar{S} openings angle (Ks,#Lambda or #bar{#Lambda}); #entries ",80,-4,4);
     histos_th1f["h_RECO_AntiS_deltaEtaDaughters_L3"] = dir_RECO_AntiS_cut_variables_distributions_L3.make<TH1F>(b+"h_RECO_AntiS_deltaEtaDaughters_L3", ";S or #bar{S} #Delta#Eta (Ks,#Lambda or #bar{#Lambda}); #entries ",80,-4,4);
     histos_th1f["h_RECO_AntiS_dxy_over_lxy_L3"] = dir_RECO_AntiS_cut_variables_distributions_L3.make<TH1F>(b+"h_RECO_AntiS_dxy_over_lxy_L3", ";dxy S or #bar{S}/lxy interaction vertex; #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_dxy_over_lxy_Ks_L3"] = dir_RECO_AntiS_cut_variables_distributions_L3.make<TH1F>(b+"h_RECO_AntiS_dxy_over_lxy_Ks_L3", ";dxy Ks/lxy interaction vertex; #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_dxy_over_lxy_AntiL_L3"] = dir_RECO_AntiS_cut_variables_distributions_L3.make<TH1F>(b+"h_RECO_AntiS_dxy_over_lxy_AntiL_L3", ";dxy #Lambda or #bar{#Lambda}/lxy interaction vertex; #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_pt_AntiL_L3"] = dir_RECO_AntiS_cut_variables_distributions_L3.make<TH1F>(b+"h_RECO_AntiS_pt_AntiL_L3", ";pt #Lambda or #bar{#Lambda} (GeV); #entries ",100,0,10);
     histos_th1f["h_RECO_AntiS_dzPVmin_L3"] = dir_RECO_AntiS_cut_variables_distributions_L3.make<TH1F>(b+"h_RECO_AntiS_dzPVmin_L3", "; S or #bar{S} min dz(PV); #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_vz_L3"] = dir_RECO_AntiS_cut_variables_distributions_L3.make<TH1F>(b+"h_RECO_AntiS_vz_L3", "; vz S or #bar{S} interaction vertex; #entries ",2000,-100,100);
     histos_th1f["h_RECO_AntiS_eta_L3"] = dir_RECO_AntiS_cut_variables_distributions_L3.make<TH1F>(b+"h_RECO_AntiS_eta_L3", "; #bar{S} #eta; #entries ",200,-10,10);


     TFileDirectory dir_RECO_AntiS_cut_variables_distributions_L4 = dir_RECO_AntiS.mkdir("RECO_AntiS_cut_variables_distributions_L4");
     histos_th1f["h_RECO_AntiS_lxy_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4.make<TH1F>(b+"h_RECO_AntiS_lxy_L4", "; #bar{S} lxy(beamspot, #bar{S} vertex) (cm); #entries ",200,0,200);
     histos_th1f["h_RECO_AntiS_error_lxy_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4.make<TH1F>(b+"h_RECO_AntiS_error_lxy_L4", "; #bar{S} error lxy(beamspot, #bar{S} vertex) (cm); #entries ",200,0,10);
     histos_th1f["h_RECO_AntiS_mass_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4.make<TH1F>(b+"h_RECO_AntiS_mass_L4", "; S or #bar{S} mass (GeV); #entries ",200,0,10);
     histos_th1f["h_RECO_AntiS_vertex_chi2_ndof_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4.make<TH1F>(b+"h_RECO_AntiS_vertex_chi2_ndof_L4", "; S or #bar{S} vertex #chi{^2}/ndof; #entries ",200,0,20);
     histos_th1f["h_RECO_AntiS_deltaPhiDaughters_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4.make<TH1F>(b+"h_RECO_AntiS_deltaPhiDaughters_L4", "; #bar{S} #Delta#Phi (Ks,#bar{#Lambda}); #entries ",100,-5,5);
     histos_th1f["h_RECO_AntiS_openingsAngleDaughters_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4.make<TH1F>(b+"h_RECO_AntiS_openingsAngleDaughters_L4", ";S or #bar{S} openings angle (Ks,#Lambda or #bar{#Lambda}); #entries ",80,-4,4);
     histos_th1f["h_RECO_AntiS_deltaEtaDaughters_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4.make<TH1F>(b+"h_RECO_AntiS_deltaEtaDaughters_L4", ";S or #bar{S} #Delta#Eta (Ks,#Lambda or #bar{#Lambda}); #entries ",80,-4,4);
     histos_th1f["h_RECO_AntiS_dxy_over_lxy_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4.make<TH1F>(b+"h_RECO_AntiS_dxy_over_lxy_L4", ";dxy S or #bar{S}/lxy interaction vertex; #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_dxy_over_lxy_Ks_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4.make<TH1F>(b+"h_RECO_AntiS_dxy_over_lxy_Ks_L4", ";dxy Ks/lxy interaction vertex; #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_dxy_over_lxy_AntiL_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4.make<TH1F>(b+"h_RECO_AntiS_dxy_over_lxy_AntiL_L4", ";dxy #Lambda or #bar{#Lambda}/lxy interaction vertex; #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_pt_AntiL_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4.make<TH1F>(b+"h_RECO_AntiS_pt_AntiL_L4", ";pt #Lambda or #bar{#Lambda} (GeV); #entries ",100,0,10);
     histos_th1f["h_RECO_AntiS_dzPVmin_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4.make<TH1F>(b+"h_RECO_AntiS_dzPVmin_L4", "; S or #bar{S} min dz(PV); #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_vz_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4.make<TH1F>(b+"h_RECO_AntiS_vz_L4", "; vz S or #bar{S} interaction vertex; #entries ",2000,-100,100);
     histos_th1f["h_RECO_AntiS_eta_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4.make<TH1F>(b+"h_RECO_AntiS_eta_L4", "; #bar{S} #eta; #entries ",200,-10,10);


     TFileDirectory dir_RECO_AntiS_cut_variables_distributions_L4_extraplots = dir_RECO_AntiS.mkdir("dir_RECO_AntiS_cut_variables_distributions_L4_extraplot");
     histos_th1f["h_RECO_AntiS_corr_dxy_Ks_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_corr_dxy_Ks_extraplots_L4", ";dxy(Ks); #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_corr_dz_Ks_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_corr_dz_Ks_extraplots_L4", ";dz(Ks); #entries ",1000,-100,100);
     histos_th1f["h_RECO_AntiS_corr_dxy_AntiL_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_corr_dxy_AntiL_extraplots_L4", ";dxy(#Lambda or #bar{#Lambda});  #entries ",100,-10,10);
     histos_th1f["h_RECO_AntiS_corr_dz_AntiL_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_corr_dz_AntiL_extraplots_L4", ";dz(#Lambda or #bar{#Lambda}); #entries ",1000,-100,100);
     histos_th2f["h2_RECO_AntiS_corr_dxy_daughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_corr_dxy_daughters_extraplots_L4", ";dxy(Ks)/lxy(Ks); dxy(#Lambda or #bar{#Lambda})/lxy(#Lambda or #bar{#Lambda}); #entries ",100,-10,10,100,-10,10);
     histos_th2f["h2_RECO_AntiS_corr_dxy_sign_dotProdLxydxy_daughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_corr_dxy_sign_dotProdLxydxy_daughters_extraplots_L4",";|dxy(Ks).|sign(#vec{lxy(Ks)}.#vec{dxy(Ks)});  |dxy(#bar{#Lambda} or #Lambda).|sign(#vec{lxy(#bar{#Lambda} or #Lambda)}.#vec{dxy(#bar{#Lambda} or #Lambda)})  ;#Entries",100,-10,10,100,-10,10);
     histos_th2f["h2_RECO_AntiS_corr_dxy_sign_dotProdptdxy_daughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_corr_dxy_sign_dotProdptdxy_daughters_extraplots_L4",";|dxy(Ks).|sign(#vec{pt(Ks)}.#vec{dxy(Ks)});  |dxy(#bar{#Lambda} or #Lambda).|sign(#vec{pt(#bar{#Lambda} or #Lambda)}.#vec{dxy(#bar{#Lambda} or #Lambda)})  ;#Entries",100,-10,10,100,-10,10);

     histos_th2f["h2_RECO_AntiS_corr_dxy_over_lxy_daughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_corr_dxy_over_lxy_daughters_extraplots_L4", ";dxy(Ks)/lxy(Ks); dxy(#Lambda or #bar{#Lambda})/lxy(#Lambda or #bar{#Lambda}); #entries ",100,-10,10,100,-10,10);
     histos_th2f["h2_RECO_AntiS_corr_dz_daughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_corr_dz_daughters_extraplots_L4", ";dz(Ks); dz(#Lambda or #bar{#Lambda}); #entries ",200,-100,100,200,-100,100);
     histos_th2f["h2_RECO_AntiS_corr_dxyz_daughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_corr_dxyz_daughters_extraplots_L4", ";dxyz(Ks); dxyz(#Lambda or #bar{#Lambda}); #entries ",200,-100,100,200,-100,100);

     histos_th2f["h2_RECO_AntiS_dxy_Ks_vs_lxy_Ks_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_dxy_Ks_vs_lxy_Ks_extraplots_L4", ";lxy(Ks); dxy(Ks); #entries ",200,-100,100,200,-10,10);
     histos_th2f["h2_RECO_AntiS_dxy_AntiL_vs_lxy_AntiL_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_dxy_AntiL_vs_lxy_AntiL_extraplots_L4", ";lxy(#Lambda or #bar{#Lambda}); dxy(#Lambda or #bar{#Lambda}); #entries ",200,-100,100,200,-10,10);

     histos_th2f["h2_RECO_AntiS_dxy_Ks_dz_Ks_daughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_dxy_Ks_dz_Ks_daughters_extraplots_L4", ";dxy(Ks); dz(Ks); #entries ",100,-10,10,200,-100,100);
     histos_th2f["h2_RECO_AntiS_dxy_L_dz_L_daughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_dxy_L_dz_L_daughters_extraplots_L4", ";dxy(#Lambda or #bar{#Lambda}); dz(#Lambda or #bar{#Lambda}); #entries ",100,-10,10,200,-100,100);
     histos_th2f["h2_RECO_AntiS_dxy_Ks_dz_L_daughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_dxy_Ks_dz_L_daughters_extraplots_L4", ";dxy(Ks); dz(#Lambda or #bar{#Lambda}); #entries ",100,-10,10,200,-100,100);
     histos_th2f["h2_RECO_AntiS_dxy_L_dz_Ks_daughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_dxy_L_dz_Ks_daughters_extraplots_L4", ";dxy(#Lambda or #bar{#Lambda}); dz(Ks); #entries ",100,-10,10,200,-100,100);


     histos_th1f["h_RECO_AntiS_relDiff_RECO_dxy_daughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_relDiff_RECO_dxy_daughters_extraplots_L4",";#frac{|dxy(Ks)|-|dxy(#bar{#Lambda})|}{|dxy(Ks)|+|dxy(#bar{#Lambda})|};#Entries",500,-10,10);
     histos_th1f["h_RECO_AntiS_relDiff_RECO_dz_daughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_relDiff_RECO_dz_daughters_extraplots_L4",";#frac{|dz(Ks)|-|dz(#bar{#Lambda})|}{|dz(Ks)|+|dz(#bar{#Lambda})|};#Entries",500,-10,10);
     histos_th1f["h_RECO_AntiS_relDiff_RECO_dxyz_daughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_relDiff_RECO_dxyz_daughters_extraplots_L4",";#frac{|dxyz(Ks)|-|dxyz(#bar{#Lambda})|}{|dxyz(Ks)|+|dxyz(#bar{#Lambda})|};#Entries",500,-10,10);
     histos_th1f["h_RECO_AntiS_RECO_vertexNormalizedChi2_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_RECO_vertexNormalizedChi2_extraplots_L4",";#Chi^{2}/ndof S or #bar{S} interaction vertex;#Entries",300,0,15);
     histos_th2f["h2_RECO_AntiS_deltaPhiDeltaEta_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_deltaPhiDeltaEta_extraplots_L4",";#Delta #Phi(Ks,#Lambda or #bar{Lambda}); #Delta #eta(Ks,#Lambda or #bar{Lambda});#Entries",80,-4,4,80,-4,4);
     histos_th1f["h_RECO_AntiS_openingsAngleDaughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_openingsAngleDaughters_extraplots_L4", ";openings angle(Ks, #Lambda or #bar{#Lambda}); #entries ",80,-4,4);
     histos_th2f["h2_RECO_AntiS_deltaEta_openingsAngleDaughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_deltaEta_openingsAngleDaughters_extraplots_L4", ";|#Delta#eta(Ks, #Lambda or #bar{#Lambda})|;openings angle(Ks, #Lambda or #bar{Lambda}); #entries ",40,0,4,40,0,4);
     histos_th1f["h_RECO_AntiS_deltaRDaughters_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_deltaRDaughters_extraplots_L4", ";#DeltaR(Ks, #Lambda or #bar{#Lambda}); #entries ",100,0,10);
     histos_th1f["h_RECO_AntiS_pt_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_pt_extraplots_L4", ";pt S or #bar{S}; #entries ",100,0,10);
     histos_th1f["h_RECO_AntiS_Ks_pt_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_Ks_pt_extraplots_L4", ";pt Ks; #entries ",100,0,10);
     histos_th1f["h_RECO_AntiS_AntiL_pt_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_AntiL_pt_extraplots_L4", ";pt #Lambda or #bar{Lambda}; #entries ",100,0,10);
     histos_th1f["h_RECO_AntiS_mass_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_mass_extraplots_L4", ";mass S or #bar{S}; #entries ",2000,0,20);
     histos_th1f["h_RECO_AntiS_eta_extraplots_L4_dxy_dz_daughters"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_eta_extraplots_L4_dxy_dz_daughters", "; #bar{S} #eta; #entries ",200,-10,10);
     histos_th1f["h_RECO_AntiS_phi_extraplots_L4_dxy_dz_daughters"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_phi_extraplots_L4_dxy_dz_daughters", "; #bar{S} #phi (rad); #entries ",200,-4,4);
     histos_th2f["h2_RECO_AntiS_vx_vy_extraplots_L4_dxy_dz_daughters"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_vx_vy_extraplots_L4_dxy_dz_daughters", "; #bar{S} vx (cm); #bar{S} vy (cm)",400,-200,200,400,-200,200);
     histos_th2f["h2_RECO_AntiS_vx_vz_extraplots_L4_dxy_dz_daughters"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH2F>(b+"h2_RECO_AntiS_vx_vz_extraplots_L4_dxy_dz_daughters", "; #bar{S} vx (cm); #bar{S} vz (cm)",400,-200,200,800,-400,400);
     histos_th1f["h_RECO_AntiS_lxy_extraplots_L4_dxy_dz_daughters"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_lxy_extraplots_L4_dxy_dz_daughters", "; #bar{S} lxy(beamspot, #bar{S} vertex) (cm); #entries ",200,0,200);
     histos_th1f["h_RECO_AntiS_vz_extraplots_L4_dxy_dz_daughters"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_vz_extraplots_L4_dxy_dz_daughters", "; #bar{S} vz(Ks vertex) (cm); #entries ",600,-300,300);
     histos_th1f["h_RECO_AntiS_dxy_extraplots_L4_dxy_dz_daughters"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_dxy_extraplots_L4_dxy_dz_daughters", "; #bar{S} dxy(beamspot) (cm); #entries ",400,-20,20);
     histos_th1f["h_RECO_AntiS_deltaR_extraplots_L4_dxy_dz_daughters"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_deltaR_extraplots_L4_dxy_dz_daughters", "; #DeltaR(Ks, #Lambda or #bar{Lambda}) ; #entries ",100,0,10);
     histos_th1f["h_RECO_AntiS_mass_extraplots_L4_dxy_dz_daughters"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_mass_extraplots_L4_dxy_dz_daughters", "; S or #bar{S} mass (GeV); #entries ",2000,0,20);
     histos_th1f["h_RECO_AntiS_dxyAntiSPVmin_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_dxyAntiSPVmin_extraplots_L4", "; S or #bar{S} dxy(PV with smallest dz)  (GeV); #entries ",200,-10,10);
     histos_th1f["h_RECO_AntiS_dzAntiSPVmin_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_dzAntiSPVmin_extraplots_L4", "; S or #bar{S} dz(PV with smallest dz)  (GeV); #entries ",200,-10,10);
     histos_th1f["h_RECO_AntiS_dzAntiSPVmin_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_dzAntiSPVmin_extraplots_L4", "; S or #bar{S} dz(PV with smallest dz)  (GeV); #entries ",200,-10,10);
     histos_th1f["h_RECO_AntiS_pz_extraplots_L4"] = dir_RECO_AntiS_cut_variables_distributions_L4_extraplots.make<TH1F>(b+"h_RECO_AntiS_pz_extraplots_L4", ";pz S or #bar{S}; #entries ",200,-50,50);
 
     TFileDirectory dir_RECO_AntiS_cutFlows = dir_RECO_AntiS.mkdir("RECO_AntiS_cutFlows");
     histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"] = dir_RECO_AntiS_cutFlows.make<TH1F>(b+"h_RECO_AntiS_BackgroundCuts_cutFlow", "; see code for definitions; #entries ",20,0,20); 

     TFileDirectory dir_RECO_AntiS_deltaRInvestigation = dir_RECO_AntiS.mkdir("RECO_AntiS_deltaRInvestigation");
     histos_th1f["h_RECO_AntiS_deltaR_cut_error_lxy"] = dir_RECO_AntiS_deltaRInvestigation.make<TH1F>(b+"h_RECO_AntiS_deltaR_cut_error_lxy", "; #Delta R(Ks, #Lambda or #bar{#Lambda}); #entries ",100,0,10); 
     histos_th1f["h_RECO_AntiS_deltaR_cut_lxy"] = dir_RECO_AntiS_deltaRInvestigation.make<TH1F>(b+"h_RECO_AntiS_deltaR_cut_lxy", "; #Delta R(Ks, #Lambda or #bar{#Lambda}); #entries ",100,0,10); 
     histos_th1f["h_RECO_AntiS_deltaR_cut_mass"] = dir_RECO_AntiS_deltaRInvestigation.make<TH1F>(b+"h_RECO_AntiS_deltaR_cut_mass", "; #Delta R(Ks, #Lambda or #bar{#Lambda}); #entries ",100,0,10); 

    

}

void AnalyzerRECO::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {

 
  //beamspot
  edm::Handle<reco::BeamSpot> h_bs;

  //primary vertex
  edm::Handle<vector<reco::Vertex>> h_offlinePV;
  iEvent.getByToken(m_offlinePVToken, h_offlinePV);

  //general tracks
  edm::Handle<View<reco::Track>> h_generalTracks;
  iEvent.getByToken(m_generalTracksToken, h_generalTracks);

  //lambdaKshortVertexFilter sexaquark candidates
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_sCands;
  iEvent.getByToken(m_sCandsToken, h_sCands);

  //V0 Kshorts
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks;
  iEvent.getByToken(m_V0KsToken, h_V0Ks);

  //V0 Lambdas
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L;
  iEvent.getByToken(m_V0LToken, h_V0L);



  //beamspot
  TVector3 beamspot(0,0,0);
  TVector3 beamspotVariance(0,0,0);
  if(h_bs.isValid()){  
	beamspot.SetXYZ(h_bs->x0(),h_bs->y0(),h_bs->z0());
	beamspotVariance.SetXYZ(pow(h_bs->x0Error(),2),pow(h_bs->y0Error(),2),pow(h_bs->z0Error(),2));			
  }

  TVector3 FirstOfflinePV(0.,0.,0.);
  if(h_offlinePV.isValid()){ 
	FirstOfflinePV.SetX(h_offlinePV->at(0).x()); FirstOfflinePV.SetY(h_offlinePV->at(0).y()); FirstOfflinePV.SetZ(h_offlinePV->at(0).z());
	histos_th1f["h_PV_n"]->Fill(h_offlinePV->size());
	for(size_t i=0; i<h_offlinePV->size(); ++i) {
		reco::Vertex PrimVertex = h_offlinePV->at(i);		
		FillHistosPV(PrimVertex, beamspot);
	}
   }

  //loop over the RECO Ks to plot the kinematics
  if(h_V0Ks.isValid()){
      for(unsigned int i = 0; i < h_V0Ks->size(); ++i){//loop all Ks candidates
	const reco::VertexCompositeCandidate * Ks = &h_V0Ks->at(i);	
	FillHistosRECOKs(Ks,beamspot,h_offlinePV);
      }
  }

  //loop over the RECO AntiLambda to plot the kinematics
  if(h_V0L.isValid()){
      for(unsigned int i = 0; i < h_V0L->size(); ++i){//loop all Ks candidates
	const reco::VertexCompositeCandidate * L = &h_V0L->at(i);	
	FillHistosRECOLambda(L,beamspot,h_offlinePV);
      }
  }


  //loop over the RECO AntiS to plot the kinematics
  if(h_sCands.isValid()){
      for(unsigned int i = 0; i < h_sCands->size(); ++i){//loop all S candidates
	const reco::VertexCompositeCandidate * antiS = &h_sCands->at(i);	
	//the if below is important, if it is -1 you are looking at signal (antiS). If it is +1 you are looking at background (S).
	int chargeAntiProton = 1; //by default look at the background
	if(m_lookAtAntiS) chargeAntiProton = -1;//only when the m_lookAtAntiS flag is enabled look at the antiS, which has a charge of -1 (representing the charge of the proton)
	if(antiS->charge()==chargeAntiProton)FillHistosRECOAntiS(antiS, beamspot, beamspotVariance, iEvent.id().event(), h_offlinePV);
      }
  }

  if(h_generalTracks.isValid()){
	for(unsigned int i = 0; i < h_generalTracks->size(); ++i){
		const reco::Track *track =  &h_generalTracks->at(i);
		FillHistosRECOTracks(track,beamspot);
	}
  }


 } //end of analyzer

void AnalyzerRECO::FillHistosPV(reco::Vertex PrimVertex, TVector3 beamspot){
	histos_th1f["h_PV_vx"]->Fill(PrimVertex.x());
	histos_th1f["h_PV_vy"]->Fill(PrimVertex.y());
	histos_th1f["h_PV_vz"]->Fill(PrimVertex.z());
}


void AnalyzerRECO::FillHistosRECOKs(const reco::VertexCompositeCandidate * RECOKs, TVector3 beamspot, edm::Handle<vector<reco::Vertex>> h_offlinePV){
	TVector3 KsCreationVertex(RECOKs->vx(),RECOKs->vy(),RECOKs->vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,KsCreationVertex);
	TVector3 KsMomentum(RECOKs->px(),RECOKs->py(),RECOKs->pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(KsCreationVertex,KsMomentum,beamspot);
	double dz = AnalyzerAllSteps::dz_line_point(KsCreationVertex,KsMomentum,beamspot);
	TVector3 PVmin = AnalyzerAllSteps::dz_line_point_min(KsCreationVertex,KsMomentum,h_offlinePV);
	double dz_min = AnalyzerAllSteps::dz_line_point(KsCreationVertex,KsMomentum,PVmin);

	histos_th1f["h_RECO_Ks_pt"]->Fill(RECOKs->pt());	
	histos_th1f["h_RECO_Ks_eta"]->Fill(RECOKs->eta());	
	histos_th1f["h_RECO_Ks_phi"]->Fill(RECOKs->phi());	
	histos_th2f["h2_RECO_Ks_vx_vy"]->Fill(RECOKs->vx(),RECOKs->vy());
	histos_th2f["h2_RECO_Ks_vx_vz"]->Fill(RECOKs->vx(),RECOKs->vz());
	histos_th1f["h_RECO_Ks_lxy"]->Fill(Lxy);	
	histos_th1f["h_RECO_Ks_vz"]->Fill(RECOKs->vz());	
	histos_th1f["h_RECO_Ks_m"]->Fill(RECOKs->mass());
	
	if(RECOKs->mass() > 0.48 && RECOKs->mass() < 0.52){
		histos_th2f["h2_RECO_Ks_lxy_vz"]->Fill(abs(RECOKs->vz()),Lxy);
		histos_th2f["h2_RECO_Ks_pt_pz"]->Fill(RECOKs->pt(),abs(RECOKs->pz()));
		histos_th2f["h2_RECO_Ks_dxy_dz"]->Fill(dxy,dz);

		histos_th2f["h2_RECO_Ks_lxy_dz"]->Fill(abs(dz),Lxy);
		histos_th1f["h_RECO_Ks_dxy"]->Fill(dxy);
		histos_th1f["h_RECO_Ks_dz"]->Fill(dz);
		histos_th1f["h_RECO_Ks_dz_min"]->Fill(dz_min);
		histos_th2f["h2_RECO_Ks_mass_dz"]->Fill(abs(dz),RECOKs->mass());
		//define kinematic region 'A'
		//if(RECOKs->pt() < RECOKs->pz()/5 && RECOKs->pt() > 1.8 && abs(dz_min) < 1 && abs(dxy) > 0.5 && abs(dxy) < 2.4){
		//if(2 < RECOKs->pt() && RECOKs->pt() < 2.5 &&   5 < abs(RECOKs->pz()) && abs(RECOKs->pz()) < 7.5  && abs(dxy)  < 0.1  && abs(dz_min) < 1 ){
		if(abs(dz_min) < 1 && abs(dxy) < 0.06 && abs(RECOKs->eta()) < 0.5 ){
			histos_th2f["h2_RECO_Ks_lxy_vz_kinregA"]->Fill(abs(RECOKs->vz()),Lxy);
			histos_th2f["h2_RECO_Ks_pt_pz_kinregA"]->Fill(RECOKs->pt(),abs(RECOKs->pz()));
			histos_th2f["h2_RECO_Ks_dxy_dz_kinregA"]->Fill(dxy,dz_min);
			histos_th2f["h2_RECO_Ks_eta_phi_kinregA"]->Fill(RECOKs->eta(),RECOKs->phi());
		}

	}

}

void AnalyzerRECO::FillHistosRECOLambda(const reco::VertexCompositeCandidate * RECOLambda, TVector3 beamspot, edm::Handle<vector<reco::Vertex>> h_offlinePV){
	TVector3 LambdaCreationVertex(RECOLambda->vx(),RECOLambda->vy(),RECOLambda->vz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,LambdaCreationVertex);
	TVector3 LambdaMomentum(RECOLambda->px(),RECOLambda->py(),RECOLambda->pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(LambdaCreationVertex,LambdaMomentum,beamspot);
	double dz = AnalyzerAllSteps::dz_line_point(LambdaCreationVertex,LambdaMomentum,beamspot);
	TVector3 PVmin = AnalyzerAllSteps::dz_line_point_min(LambdaCreationVertex,LambdaMomentum,h_offlinePV);
	double dz_min = AnalyzerAllSteps::dz_line_point(LambdaCreationVertex,LambdaMomentum,PVmin);

	histos_th1f["h_RECO_Lambda_pt"]->Fill(RECOLambda->pt());	
	histos_th1f["h_RECO_Lambda_eta"]->Fill(RECOLambda->eta());	
	histos_th1f["h_RECO_Lambda_phi"]->Fill(RECOLambda->phi());	
	histos_th2f["h2_RECO_Lambda_vx_vy"]->Fill(RECOLambda->vx(),RECOLambda->vy());
	histos_th2f["h2_RECO_Lambda_vx_vz"]->Fill(RECOLambda->vx(),RECOLambda->vz());
	histos_th1f["h_RECO_Lambda_lxy"]->Fill(Lxy);	
	histos_th1f["h_RECO_Lambda_vz"]->Fill(RECOLambda->vz());	
	histos_th1f["h_RECO_Lambda_m"]->Fill(RECOLambda->mass());	
	
	if(RECOLambda->mass() > 1.105 && RECOLambda->mass() < 1.125){
		histos_th2f["h2_RECO_Lambda_lxy_vz"]->Fill(abs(RECOLambda->vz()),Lxy);
		histos_th2f["h2_RECO_Lambda_lxy_dz"]->Fill(abs(dz),Lxy);
		histos_th1f["h_RECO_Lambda_dxy"]->Fill(dxy);
		histos_th1f["h_RECO_Lambda_dz"]->Fill(dz);
		histos_th2f["h2_RECO_Lambda_mass_dz"]->Fill(abs(dz),RECOLambda->mass());
		//if(2 < RECOLambda->pt() && RECOLambda->pt() < 2.5 &&   5 < abs(RECOLambda->pz()) && abs(RECOLambda->pz()) < 7.5  && abs(dxy)  < 0.1  && abs(dz_min) < 1 ){
		if(abs(dz_min) < 1 && abs(dxy) < 0.06 && abs(RECOLambda->eta()) < 0.5){
			histos_th2f["h2_RECO_Lambda_lxy_vz_kinregA"]->Fill(abs(RECOLambda->vz()),Lxy);
			histos_th2f["h2_RECO_Lambda_pt_pz_kinregA"]->Fill(RECOLambda->pt(),abs(RECOLambda->pz()));
			histos_th2f["h2_RECO_Lambda_dxy_dz_kinregA"]->Fill(dxy,dz_min);
			histos_th2f["h2_RECO_Lambda_eta_phi_kinregA"]->Fill(RECOLambda->eta(),RECOLambda->phi());
		}
	}
}

void AnalyzerRECO::FillHistosRECOTracks(const reco::Track *track, TVector3 beamspot){
	histos_th1f["h_RECO_Track_vz"]->Fill(abs(track->vz()));
	histos_th1f["h_RECO_Track_dz"]->Fill(abs(track->dz()));
	histos_th1f["h_RECO_Track_dsz"]->Fill(abs(track->dsz()));
}

void AnalyzerRECO::FillHistosRECOAntiS(const reco::VertexCompositeCandidate * RECOAntiS, TVector3 beamspot, TVector3 beamspotVariance, int eventId, edm::Handle<vector<reco::Vertex>> h_offlinePV){
	//calculate some kinematic variables of the for the S:
	TVector3 AntiSInteractionVertex(RECOAntiS->vx(),RECOAntiS->vy(),RECOAntiS->vz());
	TVector3 AntiSMomentumVertex(RECOAntiS->px(),RECOAntiS->py(),RECOAntiS->pz());
	double Lxy = AnalyzerAllSteps::lxy(beamspot,AntiSInteractionVertex);
	double error_Lxy = AnalyzerAllSteps::std_dev_lxy(RECOAntiS->vx(),RECOAntiS->vy(), RECOAntiS->vertexCovariance(0,0), RECOAntiS->vertexCovariance(1,1), beamspot.X(), beamspot.Y(), beamspotVariance.X(), beamspotVariance.Y());
	TVector3 AntiSMomentum(RECOAntiS->px(),RECOAntiS->py(),RECOAntiS->pz());
	double dxy = AnalyzerAllSteps::dxy_signed_line_point(AntiSInteractionVertex,AntiSMomentum,beamspot);
	
	double deltaPhiDaughters = reco::deltaPhi(RECOAntiS->daughter(0)->phi(),RECOAntiS->daughter(1)->phi());
	double deltaEtaDaughters = RECOAntiS->daughter(0)->eta()-RECOAntiS->daughter(1)->eta();
	double deltaRDaughters = sqrt(deltaPhiDaughters*deltaPhiDaughters+deltaEtaDaughters*deltaEtaDaughters);
	reco::LeafCandidate::LorentzVector n_(0,0,0,0.939565);
	double massAntiS = (RECOAntiS->p4()-n_).mass();

	//calculate some kinematic variables for the daughters of the antiS
	TVector3 AntiSDaug0Momentum(RECOAntiS->daughter(0)->px(),RECOAntiS->daughter(0)->py(),RECOAntiS->daughter(0)->pz());
	TVector3 AntiSDaug1Momentum(RECOAntiS->daughter(1)->px(),RECOAntiS->daughter(1)->py(),RECOAntiS->daughter(1)->pz());
	double dxy_daughter0 = AnalyzerAllSteps::dxy_signed_line_point(AntiSInteractionVertex, AntiSDaug0Momentum,beamspot);
	double dxy_daughter1 = AnalyzerAllSteps::dxy_signed_line_point(AntiSInteractionVertex, AntiSDaug1Momentum,beamspot);
	
	double dz_daughter0 = AnalyzerAllSteps::dz_line_point(AntiSInteractionVertex, AntiSDaug0Momentum,beamspot);
        double dz_daughter1 = AnalyzerAllSteps::dz_line_point(AntiSInteractionVertex, AntiSDaug1Momentum,beamspot);

	double dxyz_daughter0 = AnalyzerAllSteps::dxyz_signed_line_point(AntiSInteractionVertex, AntiSDaug0Momentum,beamspot);
	double dxyz_daughter1 = AnalyzerAllSteps::dxyz_signed_line_point(AntiSInteractionVertex, AntiSDaug1Momentum,beamspot);

	//collapse the RECO_dxy, RECO_dz, RECO_dxyz variables on one variable
	double relDiff_dxy_daughters = (abs(dxy_daughter0)-abs(dxy_daughter1))/(abs(dxy_daughter0)+abs(dxy_daughter1));
	double relDiff_dz_daughters = (abs(dz_daughter0)-abs(dz_daughter1))/(abs(dz_daughter0)+abs(dz_daughter1));
	double relDiff_dxyz_daughters = (abs(dxyz_daughter0)-abs(dxyz_daughter1))/(abs(dxyz_daughter0)+abs(dxyz_daughter1));

	//calculate dxy and dz for the antiS itself
	double dxy_AntiS = AnalyzerAllSteps::dxy_signed_line_point(AntiSInteractionVertex,AntiSMomentumVertex,beamspot);
	//double dz_AntiS = AnalyzerAllSteps::dz_line_point(AntiSInteractionVertex,AntiSMomentumVertex,beamspot);

	//calculate the dxy of the antiS with respect to the PV which gives the best match:
	double dzAntiSPVmin = 999.;
        TVector3  bestPVdzAntiS;
        for(unsigned int i = 0; i < h_offlinePV->size(); ++i){
                TVector3 PV(h_offlinePV->at(i).x(),h_offlinePV->at(i).y(),h_offlinePV->at(i).z());
		double dzAntiSPV = AnalyzerAllSteps::dz_line_point(AntiSInteractionVertex,AntiSMomentumVertex,PV);
                if(abs(dzAntiSPV) < abs(dzAntiSPVmin)) {dzAntiSPVmin = dzAntiSPV; bestPVdzAntiS = PV;}
        }

        double dxyAntiSPVmin = AnalyzerAllSteps::dxy_signed_line_point(AntiSInteractionVertex,AntiSMomentumVertex,bestPVdzAntiS);

	reco::Candidate::Vector AntiSDaug0MomentumVec(RECOAntiS->daughter(0)->px(),RECOAntiS->daughter(0)->py(),RECOAntiS->daughter(0)->pz());
	reco::Candidate::Vector AntiSDaug1MomentumVec(RECOAntiS->daughter(1)->px(),RECOAntiS->daughter(1)->py(),RECOAntiS->daughter(1)->pz());
	double openingsAngleDaughters = AnalyzerAllSteps::openings_angle(AntiSDaug0MomentumVec,AntiSDaug1MomentumVec);

	//calculate the sign of the dot product between the displacement vector and the dxy vector for both the Ks and the Lambda
	double signLxyDotdxy_daughter0 = AnalyzerAllSteps::sgn(AnalyzerAllSteps::vec_dxy_line_point(AntiSInteractionVertex, AntiSDaug0Momentum,beamspot)*AntiSInteractionVertex);
	double signLxyDotdxy_daughter1 = AnalyzerAllSteps::sgn(AnalyzerAllSteps::vec_dxy_line_point(AntiSInteractionVertex, AntiSDaug1Momentum,beamspot)*AntiSInteractionVertex);
	//calculate the sign of the dot product between the displacement vector and the pt of both the Ks and the Lambda
	double signPtDotdxy_daughter0 = AnalyzerAllSteps::sgn(AnalyzerAllSteps::vec_dxy_line_point(AntiSInteractionVertex, AntiSDaug0Momentum,beamspot)*AntiSDaug0Momentum);
	double signPtDotdxy_daughter1 = AnalyzerAllSteps::sgn(AnalyzerAllSteps::vec_dxy_line_point(AntiSInteractionVertex, AntiSDaug1Momentum,beamspot)*AntiSDaug1Momentum);

	//plot kinematics of *all* the RECO S candidates, this is before the background cuts
	histos_th1f["h_RECO_AntiS_pt"]->Fill(RECOAntiS->pt());	
	histos_th1f["h_RECO_AntiS_phi"]->Fill(RECOAntiS->phi());	
	histos_th2f["h2_RECO_AntiS_vx_vy"]->Fill(RECOAntiS->vx(),RECOAntiS->vy());
	histos_th2f["h2_RECO_AntiS_vx_vz"]->Fill(RECOAntiS->vx(),RECOAntiS->vz());
	histos_th1f["h_RECO_AntiS_dxy"]->Fill(dxy);	
	histos_th1f["h_RECO_AntiS_m"]->Fill(massAntiS);
	histos_th1f["h_RECO_AntiS_deltaRDaughters"]->Fill(deltaRDaughters);

	//try to find out which of the fundamental background cuts get rid of the deltaR=0 peak--> it appears to be the M>0 cut
	if(error_Lxy < 0.1) histos_th1f["h_RECO_AntiS_deltaR_cut_error_lxy"]->Fill(deltaRDaughters);	
	if(Lxy < 1.9) histos_th1f["h_RECO_AntiS_deltaR_cut_lxy"]->Fill(deltaRDaughters);	
	if(massAntiS > 0.) histos_th1f["h_RECO_AntiS_deltaR_cut_mass"]->Fill(deltaRDaughters);


	//Fill a cut flow histogram:
    	histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(0);
	if(Lxy > AnalyzerAllSteps::MinLxyCut)                               histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(1);
	if(error_Lxy < AnalyzerAllSteps::MaxErrorLxyCut)                    histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(2);
	if(massAntiS > AnalyzerAllSteps::MinErrorMassCut)                               histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(3);
	if(RECOAntiS->vertexNormalizedChi2() < AnalyzerAllSteps::MaxNormChi2Cut)        histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(4);

	if(abs(deltaPhiDaughters) > AnalyzerAllSteps::MinAbsDeltaPhiDaughtersCut)           histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(5);
	if(abs(deltaPhiDaughters) < AnalyzerAllSteps::MaxAbsDeltaPhiDaughtersCut)           histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(6);
	if(openingsAngleDaughters < AnalyzerAllSteps::MaxOpeningsAngleCut)                  histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(7);
	if(abs(deltaEtaDaughters)<AnalyzerAllSteps::MaxAbsDeltaEtaDaughCut)                 histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(8);
	if(dxy_AntiS/Lxy > AnalyzerAllSteps::MinDxyOverLxyCut)       histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(9);
	if(dxy_AntiS/Lxy < AnalyzerAllSteps::MaxDxyOverLxyCut)       histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(10);

	if(dxy_daughter1/Lxy>AnalyzerAllSteps::DxyKsExclusionRangeMaxCut || dxy_daughter1/Lxy<AnalyzerAllSteps::DxyKsExclusionRangeMinCut)histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(11);
	if(dxy_daughter0/Lxy>AnalyzerAllSteps::DxyAntiLExclusionRangeMaxCut || dxy_daughter0/Lxy<AnalyzerAllSteps::DxyAntiLExclusionRangeMinCut) histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(12);

	if(RECOAntiS->daughter(0)->pt() > AnalyzerAllSteps::MinLambdaPtCut)                 histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(13);
	if(abs(dzAntiSPVmin)<AnalyzerAllSteps::dzAntiSPVminCut)                         histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(14);
	if(abs(RECOAntiS->vz())>AnalyzerAllSteps::vzAntiSInteractionVertexCut)              histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(15);
	if(abs(RECOAntiS->eta())>AnalyzerAllSteps::antiSEtaCut)                             histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(16);

	if(     Lxy>AnalyzerAllSteps::MinLxyCut &&
		error_Lxy < AnalyzerAllSteps::MaxErrorLxyCut &&
		massAntiS > AnalyzerAllSteps::MinErrorMassCut &&
		RECOAntiS->vertexNormalizedChi2() < AnalyzerAllSteps::MaxNormChi2Cut &&

		abs(deltaPhiDaughters) > AnalyzerAllSteps::MinAbsDeltaPhiDaughtersCut &&
		abs(deltaPhiDaughters) < AnalyzerAllSteps::MaxAbsDeltaPhiDaughtersCut &&
		openingsAngleDaughters < AnalyzerAllSteps::MaxOpeningsAngleCut &&
		abs(deltaEtaDaughters) < AnalyzerAllSteps::MaxAbsDeltaEtaDaughCut &&
		dxy_AntiS/Lxy > AnalyzerAllSteps::MinDxyOverLxyCut &&
		dxy_AntiS/Lxy < AnalyzerAllSteps::MaxDxyOverLxyCut &&

		(dxy_daughter1/Lxy > AnalyzerAllSteps::DxyKsExclusionRangeMaxCut || dxy_daughter1/Lxy < AnalyzerAllSteps::DxyKsExclusionRangeMinCut) &&
		(dxy_daughter0/Lxy>AnalyzerAllSteps::DxyAntiLExclusionRangeMaxCut || dxy_daughter0/Lxy<AnalyzerAllSteps::DxyAntiLExclusionRangeMinCut) &&

		RECOAntiS->daughter(0)->pt() > AnalyzerAllSteps::MinLambdaPtCut &&
		abs(dzAntiSPVmin) < AnalyzerAllSteps::dzAntiSPVminCut &&
		abs(RECOAntiS->vz())>AnalyzerAllSteps::vzAntiSInteractionVertexCut &&
		abs(RECOAntiS->eta())>AnalyzerAllSteps::antiSEtaCut)histos_th1f["h_RECO_AntiS_BackgroundCuts_cutFlow"]->Fill(17);	

	
	//start investigating background discriminating variables
	//the stuff we cut on for *all* the RECO S candidates:
	histos_th1f["h_RECO_AntiS_lxy"]->Fill(Lxy);	
	histos_th1f["h_RECO_AntiS_error_lxy"]->Fill(error_Lxy);
	histos_th1f["h_RECO_AntiS_mass"]->Fill(massAntiS);
	histos_th1f["h_RECO_AntiS_vertex_chi2_ndof"]->Fill(RECOAntiS->vertexNormalizedChi2());
	histos_th1f["h_RECO_AntiS_deltaPhiDaughters"]->Fill(deltaPhiDaughters);
	histos_th1f["h_RECO_AntiS_openingsAngleDaughters"]->Fill(openingsAngleDaughters);
	histos_th1f["h_RECO_AntiS_deltaEtaDaughters"]->Fill(deltaEtaDaughters);
	histos_th1f["h_RECO_AntiS_dxy_over_lxy"]->Fill(dxy_AntiS/Lxy);
	histos_th1f["h_RECO_AntiS_dxy_over_lxy_Ks"]->Fill(dxy_daughter1/Lxy);
        histos_th1f["h_RECO_AntiS_dxy_over_lxy_AntiL"]->Fill(dxy_daughter0/Lxy);
	histos_th1f["h_RECO_AntiS_pt_AntiL"]->Fill(RECOAntiS->daughter(0)->pt());
	histos_th1f["h_RECO_AntiS_dzPVmin"]->Fill(dzAntiSPVmin);
	histos_th1f["h_RECO_AntiS_vz"]->Fill(RECOAntiS->vz());
	histos_th1f["h_RECO_AntiS_eta"]->Fill(RECOAntiS->eta());
		

	
}

void AnalyzerRECO::endJob()
{
}

AnalyzerRECO::~AnalyzerRECO()
{
}


DEFINE_FWK_MODULE(AnalyzerRECO);
