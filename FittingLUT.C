#include <TFile.h>
#include <TROOT.h>
#include <TChain.h>
#include "vector"

void FittingLUT(){
  
  //############################################
  // parameter
  bool IsChargeSidePlus = true;
  string inputRoot      = "./root/test_plus.root";
  string outputFileName = "test";
  //############################################
  
  string inputLUT;
  string outputLUT;
  string outputPdf;
  string sign;
  if(IsChargeSidePlus){
    inputLUT   = "./LUT/plusBetaEndcapRun2.lut";
    outputLUT  = "./LUT/" + outputFileName + "_newLUT_plus.lut";
    outputPdf  = "./PDF/" + outputFileName + "_fit_plus.pdf";
    sign       = ">0";
  } else {
    inputLUT  = "./LUT/minusBetaEndcapRun2.lut";
    outputLUT = "./LUT/" + outputFileName + "_newLUT_minus.lut";
    outputPdf  = "./PDF/" + outputFileName + "_fit_minus.pdf";
    sign       = "<0";
  }
  
  TString inputRootFile, outputPdfFile;
  inputRootFile.Form("%s",inputRoot.c_str());
  outputPdfFile.Form("%s",outputPdf.c_str());
  
  TFile *file = new TFile(inputRootFile);
  std::ifstream inputLutFile(inputLUT);
  std::ofstream fout(outputLUT);
  
  if(IsChargeSidePlus){
    fout << "side*charge=+ betapol2" << std::endl;
  } else {
    fout << "side*charge=- betapol2" << std::endl;
  }
    
  string line;
  float inputParamA[30][12];
  float inputParamB[30][12];
  while(getline(inputLutFile,line)){
    if(line.front() == 's') continue;
    std::istringstream lineStream(line);
    float num1,num2,num3,num4;
    lineStream >> num1 >> num2 >> num3 >> num4;
    int val1 = num1-1;
    int val2 = num2;
    inputParamA[val1][val2] = num3;
    inputParamB[val1][val2] = num4;
  //  cout << num1 << " " << num2 << " " << inputParamA[val1][val2] << " " << inputParamB[val1][val2] << endl;
  }
  
  TH2D *hs_Beta[30][12];

  int tableNum = 0;
  TCanvas *t1 = new TCanvas("t1","table1",600,600);
  t1->Divide(3,4,0.001,0.001);
  
  for(int etaBin=6;etaBin<30;etaBin++){
    for(int phiBin=0;phiBin<12;phiBin++){
      TH1D *hs_Median[25];
      Double_t x[25] = {0};
      Double_t y[25] = {0};
      Double_t xError[25] = {0};
      Double_t yError[25] = {0};
      
      if(phiBin == 0) tableNum = 0;
      tableNum = tableNum + 1;
      hs_Beta[etaBin][phiBin] = (TH2D*)file->Get(Form("hs_Beta_%d_%d",etaBin,phiBin));
      TCanvas *c6 =new TCanvas("c6","canvas6",700,700);
      c6->Divide(3,4);
      for(int iBin=1;iBin<25;iBin++){  
        hs_Median[iBin] = hs_Beta[etaBin][phiBin]->ProjectionY(Form("hs_Median[%d]",iBin),iBin*10+1,iBin*10+10);
        hs_Median[iBin]->SetTitle(Form("hs_Median_%d_%d_%d;#beta;Entry",etaBin,phiBin,iBin+1));
        
        int entries = hs_Median[iBin]->GetEntries();
        if(entries<20) continue;

        c6->cd(iBin+1);
        hs_Median[iBin]->Draw();
        double mean = 0;
        double minFit = 0;
        double maxFit = 0;
        double median = hs_Median[iBin]->GetMean();
        minFit = median - 0.005;
        if(minFit < 0) minFit = 0;
        maxFit = median + 0.005;

        TF1 *func = new TF1("func","gaus",minFit,maxFit);
        hs_Median[iBin]->Fit("func","","",minFit,maxFit);
        mean = func->GetParameter(1);

        if(iBin>1){
          double delta = 0;
          std::vector<double> vecMean;
          std::vector<double> vecDeltaY;
          double deltaY = median - mean;
          while(fabs(deltaY) > 0.01){
            delta += 0.001;
            if(delta > 0.01) break;
            minFit = median - delta;
            if(minFit < 0) minFit = 0;
            maxFit = median + delta;
            hs_Median[iBin]->Fit("func","","",minFit,maxFit);
            mean = func->GetParameter(1);
            vecMean.push_back(mean);
            deltaY = median - vecMean.back();
            vecDeltaY.push_back(deltaY);
          }
          if(delta > 0.01){
            std::vector<double>::iterator iter = std::min_element(vecDeltaY.begin(), vecDeltaY.end());
            size_t index = std::distance(vecDeltaY.begin(), iter);
            mean = vecMean.at(index);
          }

        }
        
        if(mean < 0.) continue;
        if(mean > 0.2) continue;
        x[iBin] = iBin*0.01+0.005;
        y[iBin] = mean;
        xError[iBin] = 0.005;
        yError[iBin] = func->GetParError(2);
      }
      t1->cd(tableNum);
      
      TGraphErrors *gr = new TGraphErrors(25,x,y,xError,yError);
      gr->SetTitle(Form("#beta_1/p_{T}_%d_%d;1/p_{T,offline}[GeV^{-1}];#beta[rad]",etaBin,phiBin));
      gr->SetMaximum(0.2);
      gr->SetMinimum(0.);
      gr->GetXaxis()->SetLimits(0.,0.25);
      gr->SetMarkerStyle(1);
      gr->SetMarkerSize(1.2);
      gr->Draw("AP");
      TF1 *f1 = new TF1("f1","[0]*x+[1]*x*x");
      TF1 *f2 = new TF1("f1","[0]*x+[1]*x*x");
      f1->SetLineWidth(2);
      f2->SetLineColor(kBlue);
      f2->SetLineWidth(2);
      f2->SetLineStyle(2);
      if(etaBin<24){
        if((etaBin==6 && phiBin==0) || (etaBin==6 && phiBin==1)){
          f1->SetParameter(0,inputParamA[etaBin][phiBin]);
          f1->SetParameter(1,1.0);
        } else {
          f1->SetParameter(0,inputParamA[etaBin][phiBin]);
          f1->SetParameter(1,inputParamB[etaBin][phiBin]);
        }
        f2->SetParameter(0,inputParamA[etaBin][phiBin]);
        f2->SetParameter(1,inputParamB[etaBin][phiBin]);
        f2->Draw("same");
      } else {
        f1->SetParameter(0,inputParamA[etaBin-6][phiBin-6]);
        f1->SetParameter(1,inputParamB[etaBin-6][phiBin-6]);
        f2->SetParameter(0,inputParamA[etaBin][phiBin]);
        f2->SetParameter(1,inputParamB[etaBin][phiBin]);
        f2->Draw("same");
      }
      gr->Fit(f1);
  
      float paramA = f1->GetParameter(0);
      float paramB = f1->GetParameter(1);
      TLegend *leg = new TLegend(0.20, 0.70, 0.50, 0.90, "");
      leg->AddEntry(f1,Form("Run-3 : #beta = %f #times (1/p_{T}) + %f #times (1/p_{T})^{2}",paramA,paramB),"l");
      leg->AddEntry(f2,Form("Run-2 : #beta = %f #times (1/p_{T}) + %f #times (1/p_{T})^{2}",inputParamA[etaBin][phiBin],inputParamB[etaBin][phiBin]),"l");
      float etaMin = 0.;
      float etaMax = 0.;
      etaMin = 1+0.05*etaBin;
      etaMax = 1+0.05*(etaBin+1.);
      float phiMin = 0.;
      float phiMax = 0.;
      phiMin = (M_PI/96)*phiBin;
      phiMax = (M_PI/96)*(phiBin+1.);
      leg->AddEntry((TObject*)0,Form("%.2f<|#eta|<%.2f, %.3f<|#phi|<%.3f, Q#times#eta/|#eta|%s",etaMin,etaMax,phiMin,phiMax,sign.c_str()),"");
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.03);
      leg->Draw();

      fout << etaBin+1 << " " << phiBin << " " << paramA << " " << paramB << std::endl;

    }
  
    t1->Modified();
    t1->Update();
    if(etaBin == 6) t1->Print(outputPdfFile+"(");
    else if(etaBin == 29) t1->Print(outputPdfFile+")");
    else t1->Print(outputPdfFile);

//    if(etaBin == 6) t1->Print("./PDF/FittingEtaBin7M0917_signMinus.pdf(");
//    else if(etaBin == 29) t1->Print("./PDF/FittingEtaBin7M0917_signMinus.pdf)");
//    else t1->Print("./PDF/FittingEtaBin7M0917_signMinus.pdf");
    
//    if(etaBin == 6) t1->Print("./PDF/FittingEtaBin7M0917_signPlus.pdf(");
//    else if(etaBin == 29) t1->Print("./PDF/FittingEtaBin7M0917_signPlus.pdf)");
//    else t1->Print("./PDF/FittingEtaBin7M0917_signPlus.pdf");
  
  }
//  fout.close();

}
