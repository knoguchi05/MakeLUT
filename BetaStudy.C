#define BetaStudy_cxx
#include "BetaStudy.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void BetaStudy::Loop()
{
  //   In a ROOT session, you can do:
  //      root> .L BetaStudy.C
  //      root> BetaStudy t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
  //

  //############################################
  // parameter
  bool IsManyDisplay    = true;
  bool IsChargeSidePlus = true;
  string outputFileName = "test";
  //############################################

  string rootFile;
  string pdfFile;
  if(IsChargeSidePlus){
    rootFile = "./root/" + outputFileName + "_plus.root";
    pdfFile = "./PDF/" + outputFileName + "_plus.pdf";
  } else {
    rootFile = "./root/" + outputFileName + "_minus.root";
    pdfFile = "./PDF/" + outputFileName + "_minus.pdf";
  }

  TString outputRootFile, outputPdfFile;
  outputRootFile.Form("%s",rootFile.c_str());
  outputPdfFile.Form("%s",pdfFile.c_str());

  TH2D *hs_Beta[30][12];
  for(int etaBin=0;etaBin<30;etaBin++){
    for(int phiBin=0;phiBin<12;phiBin++){
      hs_Beta[etaBin][phiBin] = new TH2D(Form("hs_Beta_%d_%d",etaBin,phiBin),Form("Beta_1/Pt_%d_%d;1/Pt;Beta",etaBin,phiBin),250,0.001,0.25,100,0,0.20);
    }
  }

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
//  Long64_t nentries = 500000;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    double Per = ((double)jentry/(double)nentries)*100.; 
    if(jentry%10000==0){
      cout << "#Event" << jentry << "/" << nentries << "||" << Per << "%/" << "100%" << "\r" << flush;
    }

    if(muon_pt->empty()) continue;
    if(SA_pt->empty()) continue;
    for(int off1=0;off1<muon_pt->size();off1++){
      for(int sa1=0;sa1<SA_pt->size();sa1++){

        double delEta = muon_eta->at(off1) - SA_roiEta->at(sa1);
        double delPhi = muon_phi->at(off1) - SA_roiPhi->at(sa1);
        double deltaR_match = sqrt( pow(delEta,2) + pow(delPhi,2) );
        if((muon_pt->at(off1)/1000.) < 10){
          if(deltaR_match >= -0.01*(muon_pt->at(off1)/1000) + 0.18) continue;

        } else {
          if(deltaR_match >= 0.08) continue;
        }

        if(SA_EndcapBeta->at(sa1) == 0) continue;

        // side and charge
        if(IsChargeSidePlus){ 
          if(muon_charge->at(off1) * SA_eta->at(sa1) < 0) continue;
        } else {
          if(muon_charge->at(off1) * SA_eta->at(sa1) > 0) continue;
        }

        double dividePt = 1./(muon_pt->at(off1)/1000.);

        //    if(SA_etaBin->at(sa1)<24) continue;
        if(SA_etaBin->at(sa1)>29) continue;

        int SAetaBin = SA_etaBin->at(sa1);
        int SAphiBin = SA_phiBin->at(sa1);
        hs_Beta[SAetaBin][SAphiBin]->Fill(dividePt,SA_EndcapBeta->at(sa1));  

      }
    }
  }


//########################################################################################
// create LUT Alg
//########################################################################################
  if(IsManyDisplay){
    TFile *fout = new TFile(outputRootFile,"recreate");
    int tableNum = 0;
    TCanvas *t1 = new TCanvas("t1","table1",1400,700);
    t1->Divide(12,6,0,0,0);
    TCanvas *c1 = new TCanvas("c1","canvas1",1400,700);
    c1->Divide(12,6,0,0,0);
    for(int etaBin=0;etaBin<6;etaBin++){
      for(int phiBin=0;phiBin<12;phiBin++){
        tableNum = tableNum + 1;
        t1->cd(tableNum);
        hs_Beta[etaBin][phiBin]->Draw("colz");
        hs_Beta[etaBin][phiBin]->SetStats(0);
        hs_Beta[etaBin][phiBin]->Write();
        t1->SetGridx();
      }
    }
    t1->Print(outputPdfFile+"(");

    tableNum = 0;
    TCanvas *t2 = new TCanvas("t2","table2",1400,700);
    t2->Divide(12,6,0,0,0);
    for(int etaBin=6;etaBin<12;etaBin++){
      for(int phiBin=0;phiBin<12;phiBin++){
        tableNum = tableNum + 1;
        t2->cd(tableNum);
        hs_Beta[etaBin][phiBin]->Draw("colz");
        hs_Beta[etaBin][phiBin]->SetStats(0);
        hs_Beta[etaBin][phiBin]->Write();
        t2->SetGridx();
      }
    }
    t2->Print(outputPdfFile);

    tableNum = 0;
    TCanvas *t3 = new TCanvas("t3","table3",1400,700);
    t3->Divide(12,6,0,0,0);
    for(int etaBin=12;etaBin<18;etaBin++){
      for(int phiBin=0;phiBin<12;phiBin++){
        tableNum = tableNum + 1;
        t3->cd(tableNum);
        hs_Beta[etaBin][phiBin]->Draw("colz");
        hs_Beta[etaBin][phiBin]->SetStats(0);
        hs_Beta[etaBin][phiBin]->Write();
        t3->SetGridx();
      }
    }
    t3->Print(outputPdfFile);

    tableNum = 0;
    TCanvas *t4 = new TCanvas("t4","table4",1400,700);
    t4->Divide(12,6,0,0,0);
    for(int etaBin=18;etaBin<24;etaBin++){
      for(int phiBin=0;phiBin<12;phiBin++){
        tableNum = tableNum + 1;
        t4->cd(tableNum);
        hs_Beta[etaBin][phiBin]->Draw("colz");
        hs_Beta[etaBin][phiBin]->SetStats(0);
        hs_Beta[etaBin][phiBin]->Write();
        t4->SetGridx();
      }
    }
    t4->Print(outputPdfFile);

    tableNum = 0;
    TCanvas *t5 = new TCanvas("t5","table5",1400,700);
    t5->Divide(12,6,0,0,0);
    for(int etaBin=24;etaBin<30;etaBin++){
      for(int phiBin=0;phiBin<12;phiBin++){
        tableNum = tableNum + 1;
        t5->cd(tableNum);
        hs_Beta[etaBin][phiBin]->Draw("colz");
        hs_Beta[etaBin][phiBin]->SetStats(0);
        hs_Beta[etaBin][phiBin]->Write();
        t5->SetGridx();
      }
    }
    fout->Close();
    t5->Print(outputPdfFile+")");
  } else {
    int tableNum = 0;
    TCanvas *table1 = new TCanvas("table1","table",1000,750);
    table1->Divide(4,3,0.001,0.001);
    for(int etaBin=0;etaBin<30;etaBin++){
      for(int phiBin=0;phiBin<12;phiBin++){
        if(phiBin == 0) tableNum = 0;
        tableNum = tableNum + 1;
        table1->cd(tableNum);
        hs_Beta[etaBin][phiBin]->Draw("colz");
        hs_Beta[etaBin][phiBin]->SetStats(0);
        hs_Beta[etaBin][phiBin]->Write();
        table1->SetGridx();
        TLegend *leg = new TLegend(0.10, 0.70, 0.4, 0.90, "");
        float etaMin = 0.;
        float etaMax = 0.;
        etaMin = 1+0.05*etaBin;
        etaMax = 1+0.05*(etaBin+1.);
        float phiMin = 0.;
        float phiMax = 0.;
        phiMin = (M_PI/96)*phiBin;
        phiMax = (M_PI/96)*(phiBin+1.);
        leg->AddEntry((TObject*)0,Form("%.2f<|eta|<%.2f,%.3f<|phi|<%.3f,Q<0",etaMin,etaMax,phiMin,phiMax),"");
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.035);
        leg->Draw();
      }
      table1->Modified();
      table1->Update();
      if(etaBin == 0) table1->Print(outputPdfFile+"(");
      else if(etaBin == 29) table1->Print(outputPdfFile+")");
      else table1->Print(outputPdfFile);
    }
  }

}
