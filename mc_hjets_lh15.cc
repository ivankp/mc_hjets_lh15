// -*- C++ -*-

// ChangeLog
// -----------------------------------------------------------------------------
// 02/06
// * filling underflow and overflow bins
// -----------------------------------------------------------------------------
// 13/06
// * bugfix
// -----------------------------------------------------------------------------
// 15/05 
// * added some WBF and WBF2 observables
// -----------------------------------------------------------------------------
// 19/05 MS
// * changed name to MC_HJETS_LH15
// * switched all histograms to NLO histograms to properly do fixed-order NLO 
//   statistical uncertainties
// * make analysis work with h->yy decays and stable higgs final state alike
// * remove all fiducial cuts
// * fix RapJets definition
// * fix sum_tau_jet histogramming
// * added exclusive jet multiplicities where missing
// * commented needless filling of underflow and overflow bins
// * change deltaphi_jj_bins to be in multiples of pi
// * added 3j observables (pT(H),pT(j) incl. and excl.)
// * grouped histogram initialisation to be readable/maintainable
// * used multiple bins to book-keep loose and tight cross sections
//   (0..dijet, 1..VBF, 2..VBF2)
// * some code cosmetics
// -----------------------------------------------------------------------------
// 21/05
// * bugfix in tight cross section histogramming
// * bugfix in histogram synchronisation
// -----------------------------------------------------------------------------
// 03/06
// * renamed some observables sensibly
// -----------------------------------------------------------------------------
// 05/06
// * add dR(y,j1) and dR(y,j2)
// * add jet veto cross sections
// -----------------------------------------------------------------------------
// 06/06
// * remove coarsly binned histograms
// * adjusted binnings of many others to half-way match statistical population
// -----------------------------------------------------------------------------
// 09/07
// * make FastJets object a member to fix memleak
// -----------------------------------------------------------------------------


#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/Logging.hh"
#include <map>
#include <sstream>

namespace Rivet {


  class MC_HJETS_LH15 : public Analysis {
  private:

    class NLOHisto1D : public YODA::Histo1D {
    private:

      YODA::Histo1D* _tmphist;
      int _current_event_number;

      void _syncHists() {
	for (size_t i=0; i<_tmphist->bins().size(); ++i) {
	  if (_tmphist->bin(i).area()) YODA::Histo1D::fillBin(i, _tmphist->bin(i).area());
	}
	if (_tmphist->overflow().sumW())  YODA::Histo1D::overflow()+=_tmphist->overflow();
	if (_tmphist->underflow().sumW()) YODA::Histo1D::underflow()+=_tmphist->underflow();
	_tmphist->reset();
      }

    public:

      NLOHisto1D(size_t nbins, double lower, double upper, const string& path) :
	YODA::Histo1D(nbins, lower, upper, path),
	_current_event_number(-1)
      {
	_tmphist = new Histo1D(nbins, lower, upper, path+"_tmp");
      }
      
      NLOHisto1D(const vector<double>& binedges, const string& path) :
	YODA::Histo1D(binedges, path),
	_current_event_number(-1)
      {
	_tmphist = new Histo1D(binedges, path+"_tmp");
      }

      ~NLOHisto1D()
      {
	delete _tmphist;
      }
	
      void fill(double x, const Event& event)
      {
	if (_current_event_number==-1)
	  _current_event_number = event.genEvent()->event_number();

	if (event.genEvent()->event_number()!=_current_event_number) {
	  _syncHists();
	  _current_event_number = event.genEvent()->event_number();
	}

	_tmphist->fill(x, event.weight());
      }

      void fillBin(size_t i, const Event& event, const double& fac)
      {
	if (_current_event_number==-1)
	  _current_event_number = event.genEvent()->event_number();

	if (event.genEvent()->event_number()!=_current_event_number) {
	  _syncHists();
	  _current_event_number = event.genEvent()->event_number();
	}

	_tmphist->fillBin(i, event.weight()*fac);
      }

      void finalize()
      {
	_syncHists();
      }
      
    };

    typedef shared_ptr<NLOHisto1D> NLOHisto1DPtr;


    NLOHisto1DPtr bookNLOHisto1D(const string& hname,
				      size_t nbins, double lower, double upper)
    {
      NLOHisto1DPtr hist(new NLOHisto1D(nbins, lower, upper, histoPath(hname)));
      addAnalysisObject(hist);
      return hist;
    }

    NLOHisto1DPtr bookNLOHisto1D(const string& hname,
				 const vector<double>& binedges)
    {
      NLOHisto1DPtr hist(new NLOHisto1D(binedges, histoPath(hname)));
      addAnalysisObject(hist);
      return hist;
    }

  private:

    double _jrap, _jR, _jpT;
    double _mH, _mHdev;
    double _wbfdyjj, _wbfmjj;
    double _dphiHjjtight;
    double _taujcut;
    FastJets _jetalgo;
    std::map<std::string,NLOHisto1DPtr> histos;

  public:
    MC_HJETS_LH15() : 
      Analysis("MC_HJETS_LH15"),
      _jrap(4.4), _jR(0.4), _jpT(30.*GeV), _mH(125.*GeV), _mHdev(1.*GeV),
      _wbfdyjj(2.8), _wbfmjj(400.*GeV), _dphiHjjtight(2.6), _taujcut(8.*GeV),
      _jetalgo(FastJets::ANTIKT, _jR)
    {}
  
    double sqr(const double& x) { return x*x; }
  
    // Will be used with p1=(1,0,0,1),p3=(1,0,0,-1). Only non-zero terms kept.
    std::complex<double> EPSTENSOR(const FourMomentum &p1,
				   const FourMomentum &p2,
				   const FourMomentum &p3,
				   const FourMomentum &p4)
    {
      return -std::complex<double>(0.,1.)
             *(-p1.z()*p2.x()*p3.E()*p4.y()+p1.z()*p2.y()*p3.E()*p4.x()
	       +p1.E()*p2.x()*p3.z()*p4.y()-p1.E()*p2.y()*p3.z()*p4.x());
    }
  
    void fillVetoCrossSection(double pT, string id, const Event & e) {
      NLOHisto1DPtr jvh = histos[id];
      int index = jvh->binIndexAt(pT);
      for (unsigned int i=index;i<jvh->numBins();++i) {
        jvh->fillBin(i,e,jvh->bin(i).xWidth());
      }
    }

    void init() {
      FinalState fs;
      IdentifiedFinalState higgses(PID::HIGGS);
      IdentifiedFinalState photons(PID::PHOTON);
      VetoedFinalState rest(fs);
      rest.addVetoOnThisFinalState(higgses);
      rest.addVetoOnThisFinalState(photons);
      addProjection(fs, "FS");
      addProjection(higgses, "Higgses");
      addProjection(photons, "Photons");
      addProjection(rest, "Rest");

      inithistos();
    }


    void inithistos() {
    
      histos["XS"] =
	bookNLOHisto1D("XS",1,0.,1.);
      histos["m_gammagamma"] =
	bookNLOHisto1D("m_gammagamma",bwspace(21,124.99,125.01,125.,0.00407));
    
      std::vector<double> H_pT_bins,H_pT_jj_bins,Hj_pT_bins,
			  H_pT_0j_excl_bins,H_pT_1j_excl_bins,
			  H_pT_2j_excl_bins,H_pT_3j_excl_bins,
			  H_pT_0j_incl_bins,H_pT_1j_incl_bins,
			  H_pT_2j_incl_bins,H_pT_3j_incl_bins,
			  jet1_pT_bins,jet2_pT_bins,jet3_pT_bins,
			  H_y_bins,jet_y_bins,
			  deltaphi_jj_bins,deltay_jj_bins,
			  dijet_mass_bins,
			  deltay_yy_bins,dR_y_j_bins,
			  tau_jet_bins,
			  HT_bins,
			  pTt_bins,
			  deltay_H_jj_bins,
			  delta_phi2_bins,
			  H_dijet_mass_bins,
			  deltaphi_Hjj_bins;

      H_pT_bins += 0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,
		   65.,70.,75.,80.,85.,90.,95.,100.,110.,120.,130.,140.,150.,
		   160.,170.,180.,190.,200.;

      H_pT_jj_bins += 0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,
		      130.,140.,150.,160.,170.,180.,190.,200.,220.,240.,260.,
		      280.,300.,320.,340.,360.,380.,400.;

      Hj_pT_bins += 0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,
		    65.,70.,75.,80.,85.,90.,95.,100.,110.,120.,130.,140.,150.,
		    160.,170.,180.,190.,200.,220.,240.,260.,280.,300.,
		    320.,340.,360.,380.,400.;
      
      H_pT_0j_excl_bins += 0.,20.,30.,45.,200.;
      H_pT_1j_excl_bins += 0.,40.,60.,95.,200.;
      H_pT_2j_excl_bins += 0.,90.,140.,200.;
      H_pT_3j_excl_bins += 0.,90.,140.,200.;
      H_pT_0j_incl_bins += 0.,20.,30.,45.,200.;
      H_pT_1j_incl_bins += 0.,40.,60.,95.,200.;
      H_pT_2j_incl_bins += 0.,90.,140.,200.;
      H_pT_3j_incl_bins += 0.,90.,140.,200.;

      jet1_pT_bins += 0.,30.,50.,70.,100.,140.,500.;
      jet2_pT_bins += 0,30,40,50,140,500.;
      jet3_pT_bins += 0,30,50,150,500.;

      H_y_bins += 0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,
		  1.4,1.5,1.6,1.7,1.8,1.9,2.,2.1,2.2,2.3,2.4,2.5,2.7,2.9,
		  3.1,3.3,3.5,4.0,4.5,5;
      
      jet_y_bins += 0.,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25,2.5,
		    3.,3.5,4.,4.4;

      deltaphi_jj_bins += 0.,PI/16.,2.*PI/16.,3.*PI/16.,4.*PI/16.,
			  5.*PI/16.,6.*PI/16.,7.*PI/16.,8.*PI/16.,
			  9.*PI/16.,10.*PI/16.,11.*PI/16.,12.*PI/16.,
			  13.*PI/16.,14.*PI/16.,15.*PI/16.,PI;

      deltay_jj_bins += 0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,
			6.5,7.,7.5,8.,8.5,9.,9.5,10.;

      dijet_mass_bins += 0.,50.,100.,150.,200.,250.,300.,350.,400.,450.,
			 500.,550.,600.,650.,700.,750.,800.,850.,900.,
			 950.,1000.;
      deltay_yy_bins += 0,0.3,0.6,0.9,1.2,1.5,2.0,2.55;
      dR_y_j_bins += 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,
		     2.0,2.4,2.8,3.2,3.6,4.0,4.5,5.0,5.5,6.0,7.0,8.0,10.0;

      tau_jet_bins += 0.,4.,8.,12.,16.,20.,25.,30.,40.,60.,85.;
      HT_bins += 0.,30.,40.,50.,60.,70.,90.,110.,130.,150.,200.,
		 250.,300.,400.,500.,600.,800.,1000.;
      pTt_bins += 0.,10.,20.,30.,40.,60.,80.,150.,500.;
      deltaphi_Hjj_bins += 0.0,1.0,2.0,2.3,2.6,2.8,2.9,3.,3.05,3.1,PI;

      deltay_H_jj_bins += 0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,
			  2.8,3.0,3.5,4.0,5.0,6.0,8.0;
			   
      delta_phi2_bins += -PI,-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,
			 -1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.4,
			  1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,PI;

      H_dijet_mass_bins += 0.,100.,150.,200.,250.,300.,350.,400.,450.,500.,550.,
			   600.,650.,700.,750.,800.,850.,900.,950.,1000.,1100.,
			   1200.,1300.,1400.,1500.,1700.,2000.;    

      // fine binnings always end on ""
      // incl. and excl. jet multis
      histos["NJet_incl_30"] = bookNLOHisto1D("NJet_incl_30",4,-0.5,3.5);
      histos["NJet_excl_30"] = bookNLOHisto1D("NJet_excl_30",4,-0.5,3.5);
      histos["NJet_incl_50"] = bookNLOHisto1D("NJet_incl_50",4,-0.5,3.5);
      histos["NJet_excl_50"] = bookNLOHisto1D("NJet_excl_50",4,-0.5,3.5);
      histos["NJet_incl_30_jhj"] = bookNLOHisto1D("NJet_incl_30_jhj",4,-0.5,3.5);
      histos["NJet_excl_30_jhj"] = bookNLOHisto1D("NJet_excl_30_jhj",4,-0.5,3.5);		
      histos["NJet_excl_30_VBF"] = bookNLOHisto1D("NJet_excl_30_VBF",4,-0.5,3.5);
      histos["NJet_incl_30_VBF"] = bookNLOHisto1D("NJet_incl_30_VBF",4,-0.5,3.5);
      histos["NJet_excl_30_VBF2"] = bookNLOHisto1D("NJet_excl_30_VBF2",4,-0.5,3.5);
      histos["NJet_incl_30_VBF2"] = bookNLOHisto1D("NJet_incl_30_VBF2",4,-0.5,3.5);

      // pT(H) in incl. and excl. jet bins
      histos["H_pT_incl"] = bookNLOHisto1D("H_pT_incl",H_pT_bins);
      histos["H_pT_excl"]= bookNLOHisto1D("H_pT_excl",H_pT_bins);
    
      histos["H_j_pT_incl"]= bookNLOHisto1D("H_j_pT_incl",H_pT_bins);
      histos["H_j_pT_excl"]= bookNLOHisto1D("H_j_pT_excl",H_pT_bins);
    
      histos["H_jj_pT_incl"] = bookNLOHisto1D("H_jj_pT_incl",H_pT_jj_bins);
      histos["H_jj_pT_excl"] = bookNLOHisto1D("H_jj_pT_excl",H_pT_jj_bins);

      histos["H_jjj_pT_incl"] = bookNLOHisto1D("H_jjj_pT_incl",H_pT_jj_bins);
      histos["H_jjj_pT_excl"] = bookNLOHisto1D("H_jjj_pT_excl",H_pT_jj_bins);

      histos["H_jj_pT_VBF"] = bookNLOHisto1D("H_jj_pT_VBF",H_pT_jj_bins);
      histos["H_jj_pT_VBF2"] = bookNLOHisto1D("H_jj_pT_VBF2",H_pT_jj_bins);

      // pT(H+nj) in incl. and excl. jet bins
      histos["Hj_pT_incl"]= bookNLOHisto1D("Hj_pT_incl",Hj_pT_bins);
      histos["Hj_pT_excl"]= bookNLOHisto1D("Hj_pT_excl",Hj_pT_bins);

      histos["Hjj_pT_incl"] = bookNLOHisto1D("Hjj_pT_incl",Hj_pT_bins);
      histos["Hjj_pT_excl"] = bookNLOHisto1D("Hjj_pT_excl",Hj_pT_bins);

      // pT(j) in incl. and excl. jet bins
      histos["jet1_pT_incl"] = bookNLOHisto1D("jet1_pT_incl",H_pT_bins);
      histos["jet1_pT_excl"] = bookNLOHisto1D("jet1_pT_excl",H_pT_bins);

      histos["jet2_pT_incl"] = bookNLOHisto1D("jet2_pT_incl",H_pT_bins);
      histos["jet2_pT_excl"] = bookNLOHisto1D("jet2_pT_excl",H_pT_bins);

      histos["jet3_pT_incl"] = bookNLOHisto1D("jet3_pT_incl",H_pT_bins);
      histos["jet3_pT_excl"] = bookNLOHisto1D("jet3_pT_excl",H_pT_bins);
    
      // inclusive Higgs and jet rapidities
      histos["H_y"] = bookNLOHisto1D("H_y",H_y_bins);

      histos["jet1_y"] = bookNLOHisto1D("jet1_y",jet_y_bins);

      histos["jet2_y"] = bookNLOHisto1D("jet2_y",jet_y_bins);
    
      histos["jet3_y"] = bookNLOHisto1D("jet3_y",jet_y_bins);

      // photon observables
      histos["cos_theta_star"] = bookNLOHisto1D("cos_theta_star",10,0.,1.);
      histos["cos_theta_star_80"] = bookNLOHisto1D("cos_theta_star_80",4,0.,1.);
      histos["cos_theta_star_200"] = bookNLOHisto1D("cos_theta_star_200",4,0.,1.);
      histos["cos_theta_star_gt200"] = bookNLOHisto1D("cos_theta_star_gt200",4,0.,1.);
      histos["deltay_yy"] = bookNLOHisto1D("deltay_yy",deltay_yy_bins);
      histos["dR_y_j1"] = bookNLOHisto1D("dR_y_j1",dR_y_j_bins);
      histos["dR_y_j2"] = bookNLOHisto1D("dR_y_j2",dR_y_j_bins);

      // book-keep loose and tight cross sections
      histos["loose"] = bookNLOHisto1D("loose",3,-0.5,2.5);
      histos["tight"] = bookNLOHisto1D("tight",3,-0.5,2.5);

      // \Delta\phi(jj) in incl. and excl. jet bins
      histos["deltaphi_jj_incl"] = bookNLOHisto1D("deltaphi_jj_incl",deltaphi_jj_bins);
      histos["deltaphi_jj_excl"] = bookNLOHisto1D("deltaphi_jj_excl",deltaphi_jj_bins);
      histos["deltaphi_jj_VBF"] = bookNLOHisto1D("deltaphi_jj_VBF",deltaphi_jj_bins);
      histos["deltaphi_jj_VBF2"] = bookNLOHisto1D("deltaphi_jj_VBF2",deltaphi_jj_bins);

      // \Delta y(jj) 
      histos["deltay_jj"] = bookNLOHisto1D("deltay_jj",deltay_jj_bins);

      // m(jj)
      histos["dijet_mass"] = bookNLOHisto1D("dijet_mass",dijet_mass_bins);

      // m(Hjj)
      histos["H_dijet_mass"]= bookNLOHisto1D("H_dijet_mass",H_dijet_mass_bins);

      // \Delta\phi(H,jj) incl. and excl.
      histos["deltaphi_Hjj_incl"] = bookNLOHisto1D("deltaphi_Hjj_incl",deltaphi_Hjj_bins);
      histos["deltaphi_Hjj_excl"] = bookNLOHisto1D("deltaphi_Hjj_excl",deltaphi_Hjj_bins);
      histos["deltaphi_Hjj_VBF"] = bookNLOHisto1D("deltaphi_Hjj_VBF",deltaphi_Hjj_bins);
      histos["deltaphi_Hjj_VBF2"] = bookNLOHisto1D("deltaphi_Hjj_VBF2",deltaphi_Hjj_bins);

      // \Delta y(H,jj)
      histos["deltay_H_jj"] = bookNLOHisto1D("deltay_H_jj",deltay_H_jj_bins);

      // HT
      histos["HT_all"] = bookNLOHisto1D("HT_all",HT_bins);
      histos["HT_jets"] = bookNLOHisto1D("HT_jets",HT_bins);

      // tau(j) observables
      histos["tau_jet1"] = bookNLOHisto1D("tau_jet1",tau_jet_bins);
      histos["tau_jet2"] = bookNLOHisto1D("tau_jet2",tau_jet_bins);
      histos["tau_jet3"] = bookNLOHisto1D("tau_jet3",tau_jet_bins);
      histos["tau_jet_max"]  = bookNLOHisto1D("tau_jet_max",tau_jet_bins);
      histos["sum_tau_jet"]  = bookNLOHisto1D("sum_tau_jet",tau_jet_bins);

      // pTt
      histos["pTt"] = bookNLOHisto1D("pTt",pTt_bins);

      // \phi_2
      histos["deltaphi2"] = bookNLOHisto1D("deltaphi2",delta_phi2_bins);
      histos["deltaphi2_VBF"] = bookNLOHisto1D("deltaphi2_VBF",delta_phi2_bins);
      histos["deltaphi2_VBF2"] = bookNLOHisto1D("deltaphi2_VBF2",delta_phi2_bins);
    

      // IVAN *****************************************************
				
      std::string jj[] = { "pT", "dy" };
				
      for (int i=0; i<2; ++i) {
	histos["jj"+jj[i]+"_dy"] = bookNLOHisto1D("jj"+jj[i]+"_dy",18,0,9);
	histos["jj"+jj[i]+"_dy_2j_excl"] = bookNLOHisto1D("jj"+jj[i]+"_dy_2j_excl",18,0,9);
	histos["jj"+jj[i]+"_dy_3j_excl"] = bookNLOHisto1D("jj"+jj[i]+"_dy_3j_excl",18,0,9);
      }
				
      // rapidity distance between forward and backward jets
      histos["jjfb_dy"] = bookNLOHisto1D("jjfb_dy",18,0,9);
				
      for (int i=0; i<2; ++i) {
	for (int j=1; j<=3; ++j) {
	  for (int y=1; y<=6; ++y) {
	    std::stringstream ss;
	    ss << "jet" << j << "_pT_jj" << jj[i] << "_mindy" << y;
	    histos[ss.str()] = bookNLOHisto1D(ss.str(),30,0,300);
	  }
	}
      }
    		
      histos["xs_central_jet_veto"] = bookNLOHisto1D("xs_central_jet_veto",25,0.,5.);
      histos["xs_central_jet_veto_VBF"] = bookNLOHisto1D("xs_central_jet_veto_VBF",25,0.,5.);
      histos["xs_central_jet_veto_VBF2"] = bookNLOHisto1D("xs_central_jet_veto_VBF2",25,0.,5.);
      

      // END IVAN *************************************************
      histos["xs_jet_veto_j0"]
	= bookNLOHisto1D("xs_jet_veto_j0",logspace(300,1.,1000.));
      histos["xs_jet_veto_j1_30"]
	= bookNLOHisto1D("xs_jet_veto_j1_30",logspace(300,1.,1000.));
      histos["xs_jet_veto_j1_50"]
	= bookNLOHisto1D("xs_jet_veto_j1_50",logspace(300,1.,1000.));
      histos["xs_jet_veto_j1_100"]
	= bookNLOHisto1D("xs_jet_veto_j1_100",logspace(300,1.,1000.));
      histos["xs_jet_veto_j1_200"]
	= bookNLOHisto1D("xs_jet_veto_j1_200",logspace(300,1.,1000.));
      histos["xs_jet_veto_j1_500"]
	= bookNLOHisto1D("xs_jet_veto_j1_500",logspace(300,1.,1000.));
      histos["xs_jet_veto_j1_1000"]
	= bookNLOHisto1D("xs_jet_veto_j1_1000",logspace(300,1.,1000.));
      histos["xs_jet_veto_h_50"]
	= bookNLOHisto1D("xs_jet_veto_h_50",logspace(300,1.,1000.));
      histos["xs_jet_veto_h_100"]
	= bookNLOHisto1D("xs_jet_veto_h_100",logspace(300,1.,1000.));
      histos["xs_jet_veto_h_200"]
	= bookNLOHisto1D("xs_jet_veto_h_200",logspace(300,1.,1000.));
      histos["xs_jet_veto_h_500"]
	= bookNLOHisto1D("xs_jet_veto_h_500",logspace(300,1.,1000.));
    }


    /// Do the analysis
    void analyze(const Event & e) {

      ParticleVector higgses =
	applyProjection<IdentifiedFinalState>(e, "Higgses").particles();
      ParticleVector photons =
	applyProjection<IdentifiedFinalState>(e, "Photons").particles();
      ParticleVector rest =
	applyProjection<VetoedFinalState>(e, "Rest").particles();

      // require either one stable Higgs or at least two photons
      if (higgses.size()>1) vetoEvent;
      FourMomentum hmom;
      size_t idph1(0),idph2(0);
      std::vector<FourMomentum> phs;
      if (higgses.size()==1) {
        hmom = higgses[0].momentum();
      }
      else if (photons.size()>1) {
	// reconstruct Higgs from photon pair with correct inv. mass
	// only take first one
	bool foundone(false);
        for (size_t i(0);i<photons.size();++i) {
          for (size_t j(i+1);j<photons.size();++j) {
            if (!foundone &&
	        fabs((photons[i].momentum()
                      +photons[j].momentum()).mass()-_mH)<_mHdev) {
              idph1=i; idph2=j;
              hmom = photons[i].momentum()+photons[j].momentum();
	      phs.push_back(photons[i].momentum());
	      phs.push_back(photons[j].momentum());
	      foundone=true;
	      break;
            }
          }
          if (foundone) break;
        }
      }
      else vetoEvent;

      // check that found one Higgs
      if (higgses.size()==0 && phs.size()!=2) vetoEvent;

      // add remaining photons to the remaining final state
      for (size_t i(0); i<photons.size(); ++i) {
	if (idph1==idph2 || (i!=idph1 && i!=idph2)) rest.push_back(photons[i]);
      }

      // calculate jets
      _jetalgo.calc(rest);

      // Create y ordered and pt ordered jet vectors
      Jets PTJets,RapJets,alljets;
      foreach (const Jet& jetcand, _jetalgo.pseudoJetsByPt(0.*GeV)) {
	if (fabs(jetcand.momentum().rapidity()) < _jrap) {
	  alljets.push_back(jetcand);
	}
      }
      foreach (const Jet& jetcand, _jetalgo.pseudoJetsByPt(_jpT)) {
	if (fabs(jetcand.momentum().rapidity()) < _jrap) {
	  PTJets.push_back(jetcand);
	}
      }
      foreach (const Jet& jetcand, _jetalgo.pseudojetsByRapidity(_jpT)) {
	if (fabs(jetcand.momentum().rapidity()) < _jrap) {
	  RapJets.push_back(jetcand);
	}
      }

      // cross check
      if (PTJets.size()!=RapJets.size()) abort();

      histos["XS"]->fill(0.5,e);

      // photon observables
      if (phs.size()>1) {
	histos["m_gammagamma"]->fill(hmom.mass()/GeV,e);
	// |costheta*| from 1307.1432
	double cts(abs(sinh(phs[0].eta()-phs[1].eta()))/
		  sqrt(1.+sqr(hmom.pT()/hmom.mass()))
		  * 2.*phs[0].pT()*phs[1].pT()/sqr(hmom.mass()));
	histos["cos_theta_star"]->fill(cts,e);

	if (hmom.pT()<80.*GeV) {
	  histos["cos_theta_star_80"]->fill(cts,e);
	}

	if (hmom.pT()>80.*GeV && hmom.pT()<200.*GeV) {
	  histos["cos_theta_star_200"]->fill(cts,e);
	}

	if (hmom.pT()>200.*GeV) {
	  histos["cos_theta_star_gt200"]->fill(cts,e);
	}

	double pTt = fabs(phs[0].px()*phs[1].py()-phs[1].px()*phs[0].py())/
		     ((phs[0]-phs[1]).pT()*2);
	histos["pTt"]->fill(pTt,e);
	double deltay_yy = fabs(phs[0].rapidity()-phs[1].rapidity());
	histos["deltay_yy"]->fill(deltay_yy,e);
      }

      // inclusive histograms

      histos["H_pT_incl"]->fill(hmom.pT()/GeV,e);

      histos["H_y"]->fill(fabs(hmom.rapidity()),e);

      histos["NJet_excl_30"]->fill(PTJets.size(),e);
      for (size_t i(0);i<4;++i) {
	if (PTJets.size()>=i) histos["NJet_incl_30"]->fill(i,e);
      }

      // 0j jet veto cross section
      double jv0pT1 = alljets.size()>0?alljets[0].momentum().pT():0.;
      fillVetoCrossSection(jv0pT1,"xs_jet_veto_j0",e);

      // njets == 0;
      if (PTJets.size()==0) {
	histos["H_pT_excl"]->fill(hmom.pT()/GeV,e);
	// why????
	// 6/2 added fill for jet1_pT for 0-30 GeV bin, i.e. no jets
	// histos["jet1_pT_incl"]->fill(10,e);
	// 6/2 added fill for overflow bin for 0 jets in event
	// histos["deltay_jj"]->fill(9,e);
	// 6/2 added fill for overflow bin for 0 jets in event
	// histos["Hjj_pT_incl"]->fill(160,e);
	// 6/2 added fill for overflow bin for 0 jets in event
	// histos["jet2_y"]->fill(4.6,e);
	// 6/2 added fill for overflow bin for 0 jets in event
	// histos["jet2_pT_incl"]->fill(400,e);
      }

      // njets > 0;
      if (PTJets.size()>0) {
	const FourMomentum& j1(PTJets[0].momentum());

	histos["jet1_pT_incl"]->fill(j1.pT()/GeV,e);
	histos["jet1_y"]->fill(j1.rapidity(),e);
	histos["Hj_pT_incl"]->fill((hmom+j1).pT()/GeV,e);
	histos["H_j_pT_incl"]->fill(hmom.pT()/GeV,e);

	// Calculate tau
	double tauJet1 = sqrt(sqr(j1.pT()) + sqr(j1.mass()))/
			 (2.*cosh(j1.rapidity() - hmom.rapidity()));
	
	histos["tau_jet1"]->fill(tauJet1/GeV,e);

	// njets == 1;
	if (PTJets.size()==1) {
	  histos["Hj_pT_excl"]->fill((hmom+j1).pT()/GeV,e);
	  histos["H_j_pT_excl"]->fill(hmom.pT()/GeV,e);
	  histos["jet1_pT_excl"]->fill(j1.pT()/GeV,e);
	  // again, why???
	  // 6/2 added fill for j2_pT for 0-30 GeV bins, i.e. no 2nd jet
	  // histos["jet2_pT_incl"]->fill(10,e);
	  // 6/2 added fill for overflow bin for 1 jet in event
	  // histos["deltay_jj"]->fill(9,e);
	  // 6/2 added fill for overflow bin for 1 jet in event
	  // histos["Hjj_pT_incl"]->fill(160,e);
	  // 6/2 added fill for overflow bin for 1 jet in event
	  // histos["jet2_y"]->fill(4.6,e);
	}

	// 1j jet veto cross section with minimal pT(j1)
	double jv1pT2 = alljets.size()>1?alljets[1].momentum().pT():0.;
	if (alljets[0].momentum().pT()>30.*GeV)
	  fillVetoCrossSection(jv1pT2,"xs_jet_veto_j1_30",e);
	if (alljets[0].momentum().pT()>50.*GeV)
	  fillVetoCrossSection(jv1pT2,"xs_jet_veto_j1_50",e);
	if (alljets[0].momentum().pT()>100.*GeV)
	  fillVetoCrossSection(jv1pT2,"xs_jet_veto_j1_100",e);
	if (alljets[0].momentum().pT()>200.*GeV)
	  fillVetoCrossSection(jv1pT2,"xs_jet_veto_j1_200",e);
	if (alljets[0].momentum().pT()>500.*GeV)
	  fillVetoCrossSection(jv1pT2,"xs_jet_veto_j1_500",e);
	if (alljets[0].momentum().pT()>1000.*GeV)
	  fillVetoCrossSection(jv1pT2,"xs_jet_veto_j1_1000",e);
      }

      // 1j jet veto cross section with minimal pT(j1)
      double jv1pT2 = alljets.size()>1?alljets[1].momentum().pT():0.;
      if (hmom.pT()>50.*GeV)
	fillVetoCrossSection(jv1pT2,"xs_jet_veto_h_50",e);
      if (hmom.pT()>100.*GeV)
	fillVetoCrossSection(jv1pT2,"xs_jet_veto_h_100",e);
      if (hmom.pT()>200.*GeV)
	fillVetoCrossSection(jv1pT2,"xs_jet_veto_h_200",e);
      if (hmom.pT()>500.*GeV)
	fillVetoCrossSection(jv1pT2,"xs_jet_veto_h_500",e);

      // njets > 1;
      if (PTJets.size()>1) {
	const FourMomentum& j1(PTJets[0].momentum());
	const FourMomentum& j2(PTJets[1].momentum());

	if (phs.size()>1) {
	  // fill dR of j1 and j2 and nearest photon
	  histos["dR_y_j1"]->fill(min(deltaR(j1,phs[0]),deltaR(j1,phs[1])),e);
	  histos["dR_y_j2"]->fill(min(deltaR(j2,phs[0]),deltaR(j2,phs[1])),e);
	}

	const FourMomentum& jb(RapJets.front().momentum());
	const FourMomentum& jf(RapJets.back().momentum());

    	// Calculation of phi_2 from arXiv:1001.3822  
	// phi_2 = azimuthal angle between the vector sum of jets 
	// forward and jets backward of the Higgs boson

	FourMomentum vsumf(0.,0.,0.,0.);
	FourMomentum vsumb(0.,0.,0.,0.);
	FourMomentum p1(1.,0.,0.,1.);
	FourMomentum p2(1.,0.,0.,-1.);

	bool f_nonzero(false),b_nonzero(false);
	foreach (const Jet& jj, RapJets) {
	  if (jj.momentum().rapidity()>hmom.rapidity()) {
	    vsumf += jj.momentum();
	    f_nonzero = true;
	  }
	  else {
	    vsumb += jj.momentum();
	    b_nonzero = true;
	  }
	}

	double phi2(-10.);
	// Calculate phi_2
	if (f_nonzero && b_nonzero) {
	  phi2 = acos((vsumb.x()*vsumf.x()+vsumb.y()*vsumf.y())/
		      (sqrt(vsumb.x()*vsumb.x()+vsumb.y()*vsumb.y())*
		       sqrt(vsumf.x()*vsumf.x()+vsumf.y()*vsumf.y()))); 
	  if (imag(EPSTENSOR(p1,vsumb,p2,vsumf))<0.) phi2 *= -1.;

	}
	else if (!f_nonzero) {
	  vsumb -= jf;
	  phi2 = acos((vsumb.x()*jf.x()+vsumb.y()*jf.y())/
		      (sqrt(vsumb.x()*vsumb.x()+vsumb.y()*vsumb.y())*
		       sqrt(jf.x()*jf.x()+jf.y()*jf.y())));
	  if (imag(EPSTENSOR(p1,vsumb,p2,jf))<0.)    phi2 *= -1.;
	}
	else { 
	  vsumf -= jb;
	  phi2 = acos((jb.x()*vsumf.x()+jb.y()*vsumf.y())/
		      (sqrt(jb.x()*jb.x()+jb.y()*jb.y())*
		       sqrt(vsumf.x()*vsumf.x()+vsumf.y()*vsumf.y())));  
	  if (imag(EPSTENSOR(p1,jb,p2,vsumf))<0.)    phi2 *= -1.;
	}
	histos["deltaphi2"]->fill(phi2,e);

	if (f_nonzero && b_nonzero) {
	  histos["NJet_excl_30_jhj"]->fill(PTJets.size(),e); 
	  for (size_t i(0);i<4;++i) {
	    if (PTJets.size()>=i) histos["NJet_incl_30_jhj"]->fill(i,e);   
	  }
	}

	histos["deltaphi_jj_incl"]->fill(deltaPhi(j1,j2),e);
	histos["deltaphi_Hjj_incl"]->fill(deltaPhi(hmom,j1+j2),e);
	histos["Hjj_pT_incl"]->fill((hmom+j1+j2).pT()/GeV,e);
	histos["H_jj_pT_incl"]->fill(hmom.pT()/GeV,e);
	histos["jet2_pT_incl"]->fill(j2.pT()/GeV,e);
	histos["jet2_y"]->fill(fabs(j2.rapidity()),e);
	histos["dijet_mass"]->fill((j1+j2).mass(),e);
	histos["H_dijet_mass"]->fill((hmom+j1+j2).mass(),e);
	histos["deltay_jj"]->fill(fabs(j1.rapidity()-j2.rapidity()),e);
	histos["deltay_H_jj"]->fill(fabs((hmom.rapidity()-(j1+j2).rapidity())),e);

	histos["loose"]->fill(0.,e);
	if (fabs(deltaPhi(hmom,j1+j2))>_dphiHjjtight) histos["tight"]->fill(0.,e);

	// introduce boolean whether to fill WBF and WBF2 observables
	bool wbf(fabs(j1.rapidity()-j2.rapidity())>_wbfdyjj &&
	         (j1+j2).mass()>_wbfmjj);

	bool wbf2(false);
	std::vector<Jet>::const_iterator itr1,itr2;
	for (itr1=RapJets.begin(); itr1!=RapJets.end()-1; ++itr1) {
	  for (itr2=itr1+1; itr2!=RapJets.end(); ++itr2) {
	    double dy_pair = fabs(itr1->momentum().rapidity()-itr2->momentum().rapidity());
	    double m_pair = (itr1->momentum()+itr2->momentum()).mass();
	    if (dy_pair>_wbfdyjj && m_pair>_wbfmjj) wbf2=true;
  	  }
	}

	// fill WBF histograms
	if (wbf) {
	  histos["deltaphi2_VBF"]->fill(phi2,e);
	  histos["H_jj_pT_VBF"]->fill(hmom.pT()/GeV,e);
	  histos["deltaphi_jj_VBF"]->fill(deltaPhi(j1,j2),e);
	  histos["deltaphi_Hjj_VBF"]->fill(deltaPhi(hmom,j1+j2),e);

	  histos["NJet_excl_30_VBF"]->fill(PTJets.size(),e);
	  for (size_t i(0);i<4;++i) {
	    if (PTJets.size()>=i) histos["NJet_incl_30_VBF"]->fill(i,e);
	  }

	  histos["loose"]->fill(1.,e);
	  if (fabs(deltaPhi(hmom,j1+j2))>_dphiHjjtight) histos["tight"]->fill(1.,e);
	}

	// ** 15/5/2015 **
	
	// WBF2: look for any pair that fulfils the requirements

	// fill WBF2 histograms
	if (wbf2) {
	  histos["deltaphi2_VBF2"]->fill(phi2,e);
	  histos["H_jj_pT_VBF2"]->fill(hmom.pT()/GeV,e);
	  histos["deltaphi_jj_VBF2"]->fill(deltaPhi(j1,j2),e);
	  histos["deltaphi_Hjj_VBF2"]->fill(deltaPhi(hmom,j1+j2),e);

	  histos["NJet_excl_30_VBF2"]->fill(PTJets.size(),e);
	  for (size_t i(0);i<4;++i) {
	    if (PTJets.size()>=i) histos["NJet_incl_30_VBF2"]->fill(i,e);
	  }

	  histos["loose"]->fill(2.,e);
	  if (fabs(deltaPhi(hmom,j1+j2))>_dphiHjjtight) histos["tight"]->fill(2.,e);
	}

	// END ** 15/5/2015 **

	double tauJet2 = sqrt(sqr(j2.pT()) + sqr(j2.mass()))/
	                 (2.0*cosh(j2.rapidity() - hmom.rapidity()));
	histos["tau_jet2"]->fill(tauJet2/GeV,e);

	// njets == 2;
	if (PTJets.size()==2) {
	  histos["deltaphi_jj_excl"]->fill(deltaPhi(j1,j2),e);
	  histos["deltaphi_Hjj_excl"]->fill(deltaPhi(hmom,j1+j2),e);
	  histos["Hjj_pT_excl"]->fill((hmom+j1+j2).pT()/GeV,e);
	  histos["H_jj_pT_excl"]->fill(hmom.pT()/GeV,e);
	  histos["jet2_pT_excl"]->fill(j2.pT()/GeV,e);
	  // again, why?
	  // 6/2 added fill for j3_pT 0-30 GeV, i.e. no jet3
	  // histos["jet3_pT_incl"]->fill(10,e);
	}


	// IVAN *******************************************************

	const size_t nj(PTJets.size());
	double jjpT_dy(0.), jjdy_dy(0.);

	jjpT_dy = fabs(j1.rapidity() - j2.rapidity());
	jjdy_dy = fabs(jf.rapidity() - jb.rapidity());
        
	// Fill histograms ----------------------
        
	// ** Jet pT for different dijet rapidity separations

	histos["jjpT_dy"]->fill(jjpT_dy, e);
	if      (nj==2) histos["jjpT_dy_2j_excl"]->fill(jjpT_dy, e);
	else if (nj==3) histos["jjpT_dy_3j_excl"]->fill(jjpT_dy, e);

	histos["jjdy_dy"]->fill(jjdy_dy, e);
	if      (nj==2) histos["jjdy_dy_2j_excl"]->fill(jjdy_dy, e);
	else if (nj==3) histos["jjdy_dy_3j_excl"]->fill(jjdy_dy, e);

	for (size_t j(0); j<std::min(nj,3lu); ++j) {
	  const double pT = PTJets[j].momentum().pT();

	  for (int y=1; y<=6; ++y) {
	    std::stringstream ss;
	    ss << "jet" << j+1 << "_pT_jjpT_mindy" << y;
	    if ( jjpT_dy > (j+1) ) histos[ss.str()]->fill(pT, e);
	    else break;
	  }

	  for (int y=1; y<=6; ++y) {
	    std::stringstream ss;
	    ss << "jet" << j+1 << "_pT_jjdy_mindy" << y;
	    if ( jjdy_dy > (j+1) ) histos[ss.str()]->fill(pT, e);
	    else break;
	  }
	}
  			

	// ** jet veto

	double ydists(100.);
	double ycenter=(jb.rapidity()+jf.rapidity())/2.;
	double ydistt;

	std::vector<Jet>::const_iterator jitr;
	for (jitr=RapJets.begin()+1; jitr!=RapJets.end()-1; ++jitr) {
	  ydistt=fabs(jitr->momentum().rapidity()-ycenter);
	  if (ydistt<ydists) ydists=ydistt;
	}

	// ydists is now the smallest distance between the centre of the
	// tagging jets and any possible further jet
	// (100 in case of no further jets)
	// Now fill the jet veto histogram
	NLOHisto1DPtr const hcentraljveto = histos["xs_central_jet_veto"];
	NLOHisto1DPtr const hcentraljveto_VBF = histos["xs_central_jet_veto_VBF"];
	NLOHisto1DPtr const hcentraljveto_VBF2 = histos["xs_central_jet_veto_VBF2"];

	for (int i=0, n=hcentraljveto->numBins(); i<n; ++i) {
	  if (hcentraljveto->bin(i).xMin() < ydists) {
	    double xwidth=hcentraljveto->bin(i).xWidth();
	    hcentraljveto->fillBin(i, e, xwidth);
	    
	    if (wbf)  hcentraljveto_VBF->fillBin(i, e, xwidth);
	    // ** 15/5/2015 **
	    if (wbf2) hcentraljveto_VBF2->fillBin(i, e, xwidth);
	    // END ** 15/5/2015 **
	  }
	}

	// rapidity distance between forward and backward jets
	histos["jjfb_dy"]->fill(fabs(jf.rapidity()-jb.rapidity()), e);
      
	// END IVAN ***************************************************
      }

      // njets > 2;
      if (PTJets.size()>2) {
	const FourMomentum& j3(PTJets[2].momentum());

	histos["H_jjj_pT_incl"]->fill(hmom.pT()/GeV,e);
	histos["jet3_pT_incl"]->fill(j3.pT()/GeV,e);
	histos["jet3_y"]->fill(j3.rapidity(),e);

	double tauJet3 = sqrt(sqr(j3.pT()) + sqr(j3.mass()))/
			 (2.0*cosh(j3.rapidity() - hmom.rapidity()));
	histos["tau_jet3"]->fill(tauJet3/GeV,e);

	if (PTJets.size()==3) {
	  histos["H_jjj_pT_excl"]->fill(hmom.pT()/GeV,e);
	  histos["jet3_pT_excl"]->fill(j3.pT()/GeV,e);
	}
      }

      double HT_jets(0.),HT_all(0.); 
      for(size_t i(0);i<PTJets.size();i++) HT_jets += PTJets[i].momentum().pT();
      HT_all=HT_jets+hmom.Et();
      histos["HT_jets"]->fill(HT_jets,e);
      histos["HT_all"]->fill(HT_all,e);

      double max_tj=0;
      double sum_tj=0;
      foreach (const Jet& j,PTJets) {
	double tauJet = sqrt(sqr(j.momentum().pT()) + sqr(j.momentum().mass()))/
			(2.0*cosh(j.momentum().rapidity() - hmom.rapidity()));
	if (tauJet > _taujcut) {
	  sum_tj+=tauJet;
	  if (tauJet > max_tj) max_tj=tauJet;
	}
      }
      if (max_tj>0.) histos["tau_jet_max"]->fill(max_tj,e);
      if (sum_tj>0.) histos["sum_tau_jet"]->fill(sum_tj,e);

      if (PTJets.size()>=0) {
	if (PTJets.size()>=1) {
	  const FourMomentum& j1(PTJets[0].momentum());
	  if (j1.pT()>50.*GeV) {
	    histos["NJet_incl_50"]->fill(1,e);
	  }
	  else {
	    histos["NJet_incl_50"]->fill(0,e);
	    histos["NJet_excl_50"]->fill(0,e);
	  }
	}
	if (PTJets.size()==1) {
	  const FourMomentum& j1(PTJets[0].momentum());
	  if (j1.pT()>50.*GeV) {
	    histos["NJet_excl_50"]->fill(1,e);
	  }
	}
	if (PTJets.size()>=2) {
	  const FourMomentum& j2(PTJets[1].momentum());
	  if (j2.pT()>50.*GeV) {
	    histos["NJet_incl_50"]->fill(2,e);
	  }
	}
	if (PTJets.size()==2) {
	  const FourMomentum& j2(PTJets[1].momentum());
	  if(j2.pT()>50.*GeV) {
	    histos["NJet_excl_50"]->fill(2,e);
	  }
	}
	if (PTJets.size()>=3) {
	  const FourMomentum& j3(PTJets[2].momentum());
	  if (j3.pT()>50.*GeV) {
	    histos["NJet_incl_50"]->fill(3,e);
	  }
	}
	if (PTJets.size()==3) {
	  const FourMomentum& j3(PTJets[2].momentum());
	  if (j3.pT()>50.*GeV) {
	    histos["NJet_excl_50"]->fill(3,e);
	  }
	}
      }
      if (PTJets.size()==0) {
	histos["NJet_incl_50"]->fill(0,e);
	histos["NJet_excl_50"]->fill(0,e);
      }
    }

    /// Finalize
    void finalize() {
      for (std::map<std::string,NLOHisto1DPtr>::iterator 
	     hit=histos.begin(); hit!=histos.end();hit++)
	hit->second->finalize();

      double scalefactor(crossSection()/sumOfWeights());
    
      for (std::map<std::string,NLOHisto1DPtr>::iterator 
	     hit=histos.begin(); hit!=histos.end();hit++)
	scale(hit->second,scalefactor);
    }
  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_HJETS_LH15);

}
