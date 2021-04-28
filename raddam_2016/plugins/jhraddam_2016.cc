// -*- C++ -*-
// 
// Package:    raddam_2016/raddam_2016
// Class:      raddam_2016
//
/**\class raddam_2016 raddam_2016.cc raddam_2016/raddam_2016/plugins/raddam_2016.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Mauricio Thiel
//         Created:  Mon, 03 Aug 2020 21:38:55 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"
#include "DataFormats/CTPPSReco/interface/TotemRPUVPattern.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLiteFwd.h"

#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"
#include "Geometry/VeryForwardRPTopology/interface/RPTopology.h"

#include "RecoCTPPS/TotemRPLocal/interface/FastLineRecognition.h"

#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
#include "RecoCTPPS/TotemRPLocal/interface/TotemRPLocalTrackFitterAlgorithm.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TF1.h>


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



class jhraddam2016 : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit jhraddam2016(const edm::ParameterSet&);
		~jhraddam2016();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

                TFile* efffile;
                TH2F *h_num_45F[5];
                TH2F *h_den_45F[5];
                TH2F *h_num_45N[5];
                TH2F *h_den_45N[5];
                TH2F *h_num_56F[5];
                TH2F *h_den_56F[5];
                TH2F *h_num_56N[5];
                TH2F *h_den_56N[5];

	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------
		edm::InputTag tagRecHit;

                edm::EDGetTokenT<edm::DetSetVector<TotemRPLocalTrack>> tokenLocalTrack_;
                edm::EDGetTokenT<edm::DetSetVector<TotemRPRecHit>> tokenRecHit_;

		edm::ESWatcher<VeryForwardRealGeometryRecord> geometryWatcher;

};


using namespace std;
using namespace edm;


jhraddam2016::jhraddam2016(const edm::ParameterSet& iConfig):
	tokenLocalTrack_(consumes<DetSetVector<TotemRPLocalTrack>>(iConfig.getParameter<edm::InputTag>("tagLocalTrack"))),
        tokenRecHit_(consumes< edm::DetSetVector<TotemRPRecHit> >(iConfig.getParameter<edm::InputTag>("tagRecHit")))
{
  efffile = new TFile("joutput.root", "recreate");
  efffile->cd();

  char name[100];

  for(Int_t ihist = 0; ihist < 5; ihist++)
    {
      sprintf(name,"h_num_45F_superplane%i",ihist);
      h_num_45F[ihist] = new TH2F(name,name,500,0,25,500,-25,25);
      sprintf(name,"h_den_45F_superplane%i",ihist);
      h_den_45F[ihist] = new TH2F(name,name,500,0,25,500,-25,25);

      sprintf(name,"h_num_45N_superplane%i",ihist);
      h_num_45N[ihist] = new TH2F(name,name,500,0,25,500,-25,25);
      sprintf(name,"h_den_45N_superplane%i",ihist);
      h_den_45N[ihist] = new TH2F(name,name,500,0,25,500,-25,25);

      sprintf(name,"h_num_56F_superplane%i",ihist);
      h_num_56F[ihist] = new TH2F(name,name,500,0,25,500,-25,25);
      sprintf(name,"h_den_56F_superplane%i",ihist);
      h_den_56F[ihist] = new TH2F(name,name,500,0,25,500,-25,25);

      sprintf(name,"h_num_56N_superplane%i",ihist);
      h_num_56N[ihist] = new TH2F(name,name,500,0,25,500,-25,25);
      sprintf(name,"h_den_56N_superplane%i",ihist);
      h_den_56N[ihist] = new TH2F(name,name,500,0,25,500,-25,25);
    }
}


jhraddam2016::~jhraddam2016(){ 
  for(Int_t i = 0; i < 5; i++)
    {
      h_num_45F[i]->Write();
      h_den_45F[i]->Write();
      h_num_45N[i]->Write();
      h_den_45N[i]->Write();
      h_num_56F[i]->Write();
      h_den_56F[i]->Write();
      h_num_56N[i]->Write();
      h_den_56N[i]->Write();
    }
  efffile->Write();
  efffile->Close();
}


//
// member functions
//

// ------------ method called for each event  ------------
	void
jhraddam2016::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	edm::Handle<DetSetVector<TotemRPLocalTrack>> tracks;
	iEvent.getByToken(tokenLocalTrack_, tracks);

	edm::Handle<DetSetVector<TotemRPRecHit>> hits;
	iEvent.getByToken(tokenRecHit_, hits);

	bool pass = false;
	double x_official = -999.;
	double y_official = -999.;
	double x_mine = -999.;
	double y_mine = -999.;

	
	int planes45N[10] = {0,0,0,0,0,0,0,0,0,0};
	int planes45F[10] = {0,0,0,0,0,0,0,0,0,0};
	int planes56N[10] = {0,0,0,0,0,0,0,0,0,0};
	int planes56F[10] = {0,0,0,0,0,0,0,0,0,0};

	ESHandle<CTPPSGeometry> geometry;
	iSetup.get<VeryForwardRealGeometryRecord>().get(geometry);

	for (auto &ds : *tracks)
	  {
	    CTPPSDetId rpId(ds.detId());
	    for (auto &ft : ds)
	      {
		if (!ft.isValid())
		  continue;

		double rp_z = geometry->getRPTranslation(rpId).z();
		
		//203827
		//-203827
		//212550
		//-212550

		double x = 0;
		double y = 0;

		for (unsigned int plNum = 0; plNum < 10; ++plNum)
		  {
		    TotemRPDetId plId = rpId;
		    plId.setPlane(plNum);


		    double ft_z = ft.getZ0();
		    double ft_x = ft.getX0() + ft.getTx() * (ft_z - rp_z);
		    double ft_y = ft.getY0() + ft.getTy() * (ft_z - rp_z);
		    
		    x = ft_x;
		    y = ft_y;
		    
		    double ft_v = geometry->globalToLocal(plId, CLHEP::Hep3Vector(ft_x, ft_y, ft_z)).y();
		    
		    const auto &hit_ds_it = hits->find(plId);
		    if (hit_ds_it != hits->end())
		      {
			for (const auto &h : *hit_ds_it)
			  {
			    bool match = (fabs(ft_v - h.getPosition()) < 2.*0.066);
			    if (match)
			      {
				if(rp_z == 203827)
				  planes45N[plNum] = 1;
				if(rp_z == -203827)
				  planes56N[plNum] = 1;
				if(rp_z == 212550)
				  planes45F[plNum] = 1;
				if(rp_z == -212550)
				  planes56F[plNum] = 1;
			      }
			  }
		      } 
		  }
		// 45F
		if(planes45F[2]+planes45F[3]+planes45F[4]+planes45F[5]+planes45F[6]+planes45F[7]+planes45F[8]+planes45F[9] >=4)
		  {
		    h_den_45F[0]->Fill(x,y);
		    if(planes45F[0]+planes45F[1]==2)
		      {
			h_num_45F[0]->Fill(x,y);
		      }
		  }
		if(planes45F[0]+planes45F[1]+planes45F[4]+planes45F[5]+planes45F[6]+planes45F[7]+planes45F[8]+planes45F[9] >=4)
                  {
                    h_den_45F[1]->Fill(x,y);
                    if(planes45F[2]+planes45F[3]==2)
                      {
                        h_num_45F[1]->Fill(x,y);
                      }
                  }
                if(planes45F[0]+planes45F[1]+planes45F[2]+planes45F[3]+planes45F[6]+planes45F[7]+planes45F[8]+planes45F[9] >=4)
                  {
                    h_den_45F[2]->Fill(x,y);
                    if(planes45F[4]+planes45F[4]==2)
                      {
                        h_num_45F[2]->Fill(x,y);
                      }
                  }
                if(planes45F[0]+planes45F[1]+planes45F[2]+planes45F[3]+planes45F[4]+planes45F[5]+planes45F[8]+planes45F[9] >=4)
                  {
                    h_den_45F[3]->Fill(x,y);
                    if(planes45F[6]+planes45F[7]==2)
                      {
                        h_num_45F[3]->Fill(x,y);
                      }
                  }
                if(planes45F[0]+planes45F[1]+planes45F[2]+planes45F[3]+planes45F[4]+planes45F[5]+planes45F[6]+planes45F[7] >=4)
                  {
                    h_den_45F[4]->Fill(x,y);
                    if(planes45F[8]+planes45F[9]==2)
                      {
                        h_num_45F[4]->Fill(x,y);
                      }
                  }
		// 45N
                if(planes45N[2]+planes45N[3]+planes45N[4]+planes45N[5]+planes45N[6]+planes45N[7]+planes45N[8]+planes45N[9] >=4)
                  {
                    h_den_45N[0]->Fill(x,y);
                    if(planes45N[0]+planes45N[1]==2)
                      {
                        h_num_45N[0]->Fill(x,y);
                      }
                  }
                if(planes45N[0]+planes45N[1]+planes45N[4]+planes45N[5]+planes45N[6]+planes45N[7]+planes45N[8]+planes45N[9] >=4)
                  {
                    h_den_45N[1]->Fill(x,y);
                    if(planes45N[2]+planes45N[3]==2)
                      {
                        h_num_45N[1]->Fill(x,y);
                      }
                  }
                if(planes45N[0]+planes45N[1]+planes45N[2]+planes45N[3]+planes45N[6]+planes45N[7]+planes45N[8]+planes45N[9] >=4)
                  {
                    h_den_45N[2]->Fill(x,y);
                    if(planes45N[4]+planes45N[4]==2)
                      {
                        h_num_45N[2]->Fill(x,y);
                      }
                  }
                if(planes45N[0]+planes45N[1]+planes45N[2]+planes45N[3]+planes45N[4]+planes45N[5]+planes45N[8]+planes45N[9] >=4)
                  {
                    h_den_45N[3]->Fill(x,y);
                    if(planes45N[6]+planes45N[7]==2)
                      {
                        h_num_45N[3]->Fill(x,y);
                      }
                  }
                if(planes45N[0]+planes45N[1]+planes45N[2]+planes45N[3]+planes45N[4]+planes45N[5]+planes45N[6]+planes45N[7] >=4)
                  {
                    h_den_45N[4]->Fill(x,y);
                    if(planes45N[8]+planes45N[9]==2)
                      {
                        h_num_45N[4]->Fill(x,y);
                      }
                  }
		// 56F
                if(planes56F[2]+planes56F[3]+planes56F[4]+planes56F[5]+planes56F[6]+planes56F[7]+planes56F[8]+planes56F[9] >=4)
                  {
                    h_den_56F[0]->Fill(x,y);
                    if(planes56F[0]+planes56F[1]==2)
                      {
                        h_num_56F[0]->Fill(x,y);
                      }
                  }
                if(planes56F[0]+planes56F[1]+planes56F[4]+planes56F[5]+planes56F[6]+planes56F[7]+planes56F[8]+planes56F[9] >=4)
                  {
                    h_den_56F[1]->Fill(x,y);
                    if(planes56F[2]+planes56F[3]==2)
                      {
                        h_num_56F[1]->Fill(x,y);
                      }
                  }
                if(planes56F[0]+planes56F[1]+planes56F[2]+planes56F[3]+planes56F[6]+planes56F[7]+planes56F[8]+planes56F[9] >=4)
                  {
                    h_den_56F[2]->Fill(x,y);
                    if(planes56F[4]+planes56F[4]==2)
                      {
                        h_num_56F[2]->Fill(x,y);
                      }
                  }
                if(planes56F[0]+planes56F[1]+planes56F[2]+planes56F[3]+planes56F[4]+planes56F[5]+planes56F[8]+planes56F[9] >=4)
                  {
                    h_den_56F[3]->Fill(x,y);
                    if(planes56F[6]+planes56F[7]==2)
                      {
                        h_num_56F[3]->Fill(x,y);
                      }
                  }
                if(planes56F[0]+planes56F[1]+planes56F[2]+planes56F[3]+planes56F[4]+planes56F[5]+planes56F[6]+planes56F[7] >=4)
                  {
                    h_den_56F[4]->Fill(x,y);
                    if(planes56F[8]+planes56F[9]==2)
                      {
                        h_num_56F[4]->Fill(x,y);
                      }
                  }
		// 56N
                if(planes56N[2]+planes56N[3]+planes56N[4]+planes56N[5]+planes56N[6]+planes56N[7]+planes56N[8]+planes56N[9] >=4)
                  {
                    h_den_56N[0]->Fill(x,y);
                    if(planes56N[0]+planes56N[1]==2)
                      {
                        h_num_56N[0]->Fill(x,y);
                      }
                  }
                if(planes56N[0]+planes56N[1]+planes56N[4]+planes56N[5]+planes56N[6]+planes56N[7]+planes56N[8]+planes56N[9] >=4)
                  {
                    h_den_56N[1]->Fill(x,y);
                    if(planes56N[2]+planes56N[3]==2)
                      {
                        h_num_56N[1]->Fill(x,y);
                      }
                  }
                if(planes56N[0]+planes56N[1]+planes56N[2]+planes56N[3]+planes56N[6]+planes56N[7]+planes56N[8]+planes56N[9] >=4)
                  {
                    h_den_56N[2]->Fill(x,y);
                    if(planes56N[4]+planes56N[4]==2)
                      {
                        h_num_56N[2]->Fill(x,y);
                      }
                  }
                if(planes56N[0]+planes56N[1]+planes56N[2]+planes56N[3]+planes56N[4]+planes56N[5]+planes56N[8]+planes56N[9] >=4)
                  {
                    h_den_56N[3]->Fill(x,y);
                    if(planes56N[6]+planes56N[7]==2)
                      {
                        h_num_56N[3]->Fill(x,y);
                      }
                  }
                if(planes56N[0]+planes56N[1]+planes56N[2]+planes56N[3]+planes56N[4]+planes56N[5]+planes56N[6]+planes56N[7] >=4)
                  {
                    h_den_56N[4]->Fill(x,y);
                    if(planes56N[8]+planes56N[9]==2)
                      {
                        h_num_56N[4]->Fill(x,y);
                      }
                  }
	      }
	  }

	
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
	void
jhraddam2016::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void
jhraddam2016::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
jhraddam2016::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);

	//Specify that only 'tracks' is allowed
	//To use, remove the default given above and uncomment below
	//ParameterSetDescription desc;
	//desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
	//descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(jhraddam2016);
