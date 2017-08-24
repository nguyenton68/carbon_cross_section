/* Nguyen Ton Aug 2017
 * Get simulation cross section with cuts:
 * acceptance cut, analysis cut, PID cut
 * Input: MC output root file
 * - Corners (theta,phi) which define the edge cut at target angle
 * - Change the cut according to the edge distribution
 * Output:
 * - Cross section ub
 */
#include <iomanip>
#include "TMath.h"
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TROOT.h"
#include <cstdlib>
const int CORNERS = 7;
const int Nbin=200;
void rec_sim_xs()
{

  int runnum = 1899;

  double w_min(-0.005),w_max(0.01);
  /* Focal plane histogram limits*/
  double y_min(-0.06),y_max(0.06),ph_min(-0.06),ph_max(0.06);
  /* Target plane histogram limits */
  double hph_min(-10.0),hph_max(20.0),hth_min(10.0),hth_max(70.0);

  /* Output root file and histograms */
  TFile * fout = new TFile("sim_1899_reconstructed_target.root","RECREATE");
  TH1F * hy     = new TH1F("hy","",Nbin,y_min,y_max);
  TH1F * hxs     = new TH1F("hxs","",Nbin,0,25000);
  TH2F * hyphi     = new TH2F("hyphi","",Nbin,y_min,y_max,Nbin,ph_min,ph_max);
  TH2F * hthph = new TH2F("hthph","",Nbin,hph_min,hph_max,Nbin,hth_min,hth_max);
  TH2F * hthph_cut = new TH2F("hthph_cut","",Nbin,hph_min,hph_max,Nbin,hth_min,hth_max);

  Float_t deltatg,xs,phfoc,yfoc,phitg,thetatg,xfoc,thfoc,dpor,thor,phor,wor,wmm;
  /* Input root file */
  TFile * f = new TFile("1899_elastic_wo_cut.root");
  TTree *ntup = (TTree*)gROOT->FindObject("h1");

  /* Count for cross section calculation */  
  double xsc=0.0;
  int tacc=0;
  Int_t n_entries = (Int_t)ntup->GetEntries();
  ntup->SetBranchAddress("phfoc",&(phfoc));
  ntup->SetBranchAddress("yfoc",&(yfoc));
  ntup->SetBranchAddress("xfoc",&(xfoc));
  ntup->SetBranchAddress("thfoc",&(thfoc));
  ntup->SetBranchAddress("deltat",&(deltatg));
  ntup->SetBranchAddress("thetat",&(thetatg));
  ntup->SetBranchAddress("phit",&(phitg));
  ntup->SetBranchAddress("wmm",&(wmm));
  ntup->SetBranchAddress("wor",&(wor));
  ntup->SetBranchAddress("xs",&(xs));
  ntup->SetBranchAddress("thor",&(thor));
  ntup->SetBranchAddress("phor",&(phor));
  ntup->SetBranchAddress("dpor",&(dpor));

  cout<<"Number of events= "<<n_entries<<endl;

  /* Define edge for analysis cuts:
   * These edges are determined from phase space run with analysis cuts which are applied to data
   */
  double corner_phi[CORNERS],corner_theta[CORNERS];
  double edge[CORNERS];
  corner_phi[0] = 8.764;//8.806;
  corner_phi[1] = 8.216;//8.441;
  corner_phi[2] = 1.995;//1.842;
  corner_phi[3] = 2.973;//2.895;
  corner_phi[4] = 2.973;//2.895;
  corner_phi[5] = 2.309;//2.328;
  corner_phi[6] = 3.404;//2.895;

  corner_theta[0] = 50.233;//50.;
  corner_theta[1] = 30.016;//34.113;
  corner_theta[2] = 32.306;//32.258;
  corner_theta[3] = 39.871;//39.758;
  corner_theta[4] = 45.363;//45.323;
  corner_theta[5] = 50.233;//49.677;
  corner_theta[6] = 51.788;//51.613;

 /* Acceptance cut */
  const int ACC = 6;//corner to define the edge cut
  double y_corner[ACC],ph_corner[ACC],acc_edge[ACC];

  /* Start to loop the root tree */
  for(int kk=0;kk<n_entries-1;kk++){
    ntup->GetEntry(kk);	

    /* Analysis cut */
    double slope, intercept;
    for(int ll=0; ll<CORNERS; ll++)
      {
	if(ll==(CORNERS-1)){
	  slope = (corner_theta[ll]-corner_theta[0])/(corner_phi[ll]-corner_phi[0]);
	} else {
	  slope = (corner_theta[ll+1]-corner_theta[ll])/(corner_phi[ll+1]-corner_phi[ll]);
	}
	intercept = corner_theta[ll] - slope*corner_phi[ll];
	edge[ll] = slope*phor + intercept;
      }  

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


    /* Count tacc: which is number of event wo radiation and acceptance effect
     * Cut is applied to original thrown quantities
     */
      if(thor>=edge[0]&&thor>=edge[1]
       &&(thor<=edge[2]||phor>=corner_phi[3]||thor>=edge[4])
       &&thor<=edge[5]&&thor<=edge[6]
       &&wor<=w_max&&wor>=w_min)
      tacc++;

      hthph->Fill(phor,thor);//2D histogram wo cut

      /* Acceptance cut and w (radiative) cut */
      if((phfoc<=acc_edge[0]||phfoc<=acc_edge[5])&&phfoc>=acc_edge[2]&&phfoc<=acc_edge[4]&&yfoc<=acc_edge[3]&&yfoc>=acc_edge[1]
	 &&wmm<=w_max&&wmm>=w_min)
	{
      //yfoc<=0.03&&yfoc>=-0.02&&phfoc<=0.05&&phfoc>=-0.05&&thfoc<=0.03&&thfoc>=-0.03// fp limit cuts
	  double ph_rec =0.001692 -0.012293*xfoc -0.007802*yfoc-0.40254*thfoc +0.627065*phfoc
	    -5.927*yfoc*phfoc -6.172*yfoc*thfoc
	    -11.326*yfoc*yfoc +6.111*phfoc*phfoc;
	  double th_rec = 0.037990 -0.004055*xfoc +0.849883*yfoc +0.07456*thfoc -0.163464*phfoc;
	  
	  /* Analysis cut and fill histograms */
	  if(thor>=edge[0]&&thor>=edge[1]
	     &&(thor<=edge[2]||phor>=corner_phi[3]||thor>=edge[4])
	     &&thor<=edge[5]&&thor<=edge[6]) 
	    {
	      hy->Fill(yfoc,xs);
	      hxs->Fill(xs,xs);
	      hyphi->Fill(yfoc,phfoc);
	      hthph_cut->Fill(phor,thor);
	    }// end analysis cut
	}//end acceptance and w cut
  }//end loop root tree
  
  double cacc = hy->GetEntries();//number of event survive after analysis and acceptance cut
  double cross_section = (hxs->Integral())/tacc;

  cout<<"Number of event inside sieve cut = "<<tacc<<endl;
  cout<<"Number of event inside 4D cut = "<<cacc<<endl;
  cout<<"Cross section= "<<cross_section<<endl;
  cout<<"**************"<<endl;
  TCanvas * c= new TCanvas("c","",800,600);
  c->Clear();
  hthph->Draw("colz");
  hthph_cut->Draw("same");
  fout->Write();
  
}
