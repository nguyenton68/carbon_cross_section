/* Nguyen Ton Aug 2017 
 * Get solid angle for data from MC (phase space only). The MC is ran with RC=0(OFF), phase space turn ON.
 * Input: 
 * - Root files (need enough statistic)
 * - Set of cut to define the analysis cut correspond to experimental data cut (this is exactly the same as cut applied in data)
 * Output:
 * - Solid angle in usr
 */
#include <iomanip>
#include <vector>
#include "TMath.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"

void rec_phase_space()
{
  /* Define histograms */
  TH2F * hthph     = new TH2F("hthph","",200,-20,20,200,20,60);
  TH2F * hthphCut  = new TH2F("hthphCut","",200,-20,20,200,20,60);  

  int bcsa=0;
  int scsa=0; 
  /* This analysis cut is applied in experimental data */
  const int DIM = 4;
  double c_phiTg[DIM]={0.006894,0.0005,0.005084,0.01295};
  double c_thetaTg[DIM]={0.04637,0.04616,0.03367,0.03641};
  double edge[DIM];
  
  char command[120];
  sprintf(command,"ls *.root |cat > rootlist");
  gSystem->Exec(command);
  char rootfile[100];
  ifstream rootlist("rootlist",ios_base::in);


  while(true){
    rootlist>>rootfile;
    if(rootlist.eof()) break;
  
    Float_t phfoc,yfoc,phitg,thetatg,deltatg,xfoc,thfoc;
    TFile *rootlist1  = new TFile(rootfile);//TFile::Open(rootfile);
    // TFile *fSim  = new TFile("./rootfiles/center_foil_phase_space_wo_cut.root");

    TTree *ntup = (TTree*)gROOT->FindObject("h1");
    //TTree *ntup = (TTree*)fSim->Get("h1");
    Int_t n_entries = (Int_t)ntup->GetEntries();
    ntup->SetBranchAddress("phfoc",&(phfoc));
    ntup->SetBranchAddress("yfoc",&(yfoc));
    ntup->SetBranchAddress("xfoc",&(xfoc));
    ntup->SetBranchAddress("thfoc",&(thfoc));
    ntup->SetBranchAddress("deltat",&(deltatg));
    ntup->SetBranchAddress("thetat",&(thetatg));
    ntup->SetBranchAddress("phit",&(phitg));
    cout<<"Number of events= "<<n_entries<<endl;
    /* Acceptance cut */
    const int ACC = 6;//corner to define the edge cut
    double y_corner[ACC],ph_corner[ACC],acc_edge[ACC];

    for(int kk=0;kk<n_entries-1;kk++){
      ntup->GetEntry(kk);

      /* Acceptance cuts */
      /* caution: convert delta --> percentage */
      double delta = deltatg/100.;
      y_corner[0] = -0.227484*delta + 0.00965004;
      y_corner[1] = -0.162115*delta - 0.017178;
      y_corner[2] =  0.203004*delta - 0.0157547;
      y_corner[3] =  0.524234*delta + 0.0129056;
      y_corner[4] = -0.260981*delta + 0.024819;
      y_corner[5] = -0.394433*delta + 0.0205985;

      ph_corner[0] = -0.120625*delta + 0.0271181;
      ph_corner[1] = -0.0885365*delta + 0.012184;
      ph_corner[2] = 0.209205*delta - 0.0203365;
      ph_corner[3] = 0.677652*delta - 0.0290798;
      ph_corner[4] = 0.105192*delta + 0.0445469;
      ph_corner[5] = -0.0579263*delta + 0.0470789;
      double slp, xpoint;
      for(int ii=0; ii<ACC; ii++)
	{
	  if(ii==1||ii==3)
	    {
	      slp =(y_corner[ii+1]-y_corner[ii])/(ph_corner[ii+1]-ph_corner[ii]);
	      xpoint = y_corner[ii]-slp*ph_corner[ii];
	      acc_edge[ii]=slp*phfoc+xpoint;
	    }
	  else
	    {
	      if(ii==5)
		{
		  slp=(ph_corner[0]-ph_corner[ii])/(y_corner[0]-y_corner[ii]);
		}
	      else
		{
		  slp=(ph_corner[ii+1]-ph_corner[ii])/(y_corner[ii+1]-y_corner[ii]);
		}
	      xpoint = ph_corner[ii] - slp*y_corner[ii];
	      acc_edge[ii] = slp*yfoc + xpoint;
	    }
	}
      /* Acceptance cut and w (radiative) cut */
      if((phfoc<=acc_edge[0]||phfoc<=acc_edge[5])&&phfoc>=acc_edge[2]&&phfoc<=acc_edge[4]&&yfoc<=acc_edge[3]&&yfoc>=acc_edge[1]
	 &&deltatg>-0.5&&deltatg<0.0)
	{
	  double ph_rec =0.001692 -0.012293*xfoc -0.007802*yfoc-0.40254*thfoc +0.627065*phfoc
	    -5.927*yfoc*phfoc -6.172*yfoc*thfoc
	    -11.326*yfoc*yfoc +6.111*phfoc*phfoc;
	  double th_rec = 0.037990 -0.004055*xfoc +0.849883*yfoc +0.07456*thfoc -0.163464*phfoc;
	  double slope, intercept;
	  for(int i=0; i<DIM; i++)
	    {
	      if(i==(DIM-1))
		{
		  slope =(c_thetaTg[0]-c_thetaTg[i])/(c_phiTg[0]-c_phiTg[i]);
		} else 
		{
		  slope= (c_thetaTg[i+1]-c_thetaTg[i])/(c_phiTg[i+1]-c_phiTg[i]);
		}
	      intercept = c_thetaTg[i] - slope*c_phiTg[i];
	      edge[i] = slope*ph_rec + intercept;
	    }
	  if(th_rec<=edge[0]&&th_rec>=edge[1]&&th_rec>=edge[2]&&th_rec<=edge[3])
	    {
	      hthph->Fill(phitg,thetatg);
	      bcsa++;
	      if(thetatg>=35&&thetatg<=45&&phitg>=4&&phitg<=6){
		scsa++;
		hthphCut->Fill(phitg,thetatg);
	      }
	    }// cut at reconstructed variables
	}// acceptance cut
    }// loop events
    delete rootlist1;
  }
  TCanvas * c = new TCanvas("c","",800,600);
  c->Clear();
  hthph->Draw("colz");
  
  cout<<"**************"<<endl;
  cout<<"Real solid angle= "<<bcsa*20./scsa<<endl;
    
}
