#include "MLEfit.hh"

#include <TString.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraph.h>


int MLEfit::m_fitOption = 0 ;

double MLEfit::vec_nuTime1D[100000] = {-10};
double MLEfit::vec_nuTime2D[100000] = {-10};
double MLEfit::vec_nuEnergy[100000] = {-10};

TH1D* MLEfit::TPdf;
TH2D* MLEfit::TEvisPdf;

double MLEfit::Tmin;
double MLEfit::Tmax;
double MLEfit::Emin;
double MLEfit::Emax;

// pre-setting
double MLEfit::fitTmin = 0.00; //-0.02;
double MLEfit::fitTmax = 0.02;  //0.05;
double MLEfit::extTmin = -0.05;
double MLEfit::extTmax = 0.09;  //0.08;

double MLEfit::m_statScale = 1;

int MLEfit::m_nevents = 1000;

double MLEfit::m_deltaT = 0;

double MLEfit::m_nNu = 1;
double MLEfit::m_scale1D = 1;

double MLEfit::m_chi2;
double MLEfit::m_chi2Min;
int MLEfit::m_nParameter;
double MLEfit::m_bestFit[1];
double MLEfit::m_bestFitError[1];
bool MLEfit::m_DoFit;



MLEfit::MLEfit()
{
    for(int i=0; i<100000; i++) {
        vec_nuTime1D[i] = -10;
        vec_nuTime2D[i] = -10;
        vec_nuEnergy[i] = -10;
    }

    //if (m_fitOption == 0)
    //    cout << "In this script we adopt 1D fitting " << endl;
    //if (m_fitOption == 1)
    //    cout << "In this script we adopt 2D fitting " << endl;
}

MLEfit::~MLEfit()
{;}

void MLEfit::LoadPdf()
{
    TFile* pdfFile = TFile::Open(pdfFileName, "read");
    TEvisPdf = (TH2D*)pdfFile->Get("TEvisPdf");
    //TEvisPdf = (TH2D*)pdfFile->Get("hENuT_mod82503_cha1_mh2");
    //cout << "----------> Loading 2D Pdf histogram <----------" << endl;

    TPdf = new TH1D("TPdf", "", TEvisPdf->GetNbinsX(), -0.1, 0.1);

    Emin = 0.1;
    Emax = 4;

    // projection manual
    for(int xbin=1; xbin<TEvisPdf->GetNbinsX(); xbin++) {
        double cont = 0;
        for(int ybin=1; ybin<TEvisPdf->GetNbinsY(); ybin++) {
            double tmpE = TEvisPdf->GetYaxis()->GetBinCenter(ybin);
            if (tmpE < Emin or tmpE > Emax)
                continue;
            cont += TEvisPdf->GetBinContent(xbin, ybin) * TEvisPdf->GetYaxis()->GetBinWidth(1);
        }
        cont *= m_statScale;
        TPdf->SetBinContent(xbin, cont);
    }

    //TPdf->Scale(TEvisPdf->GetYaxis()->GetBinWidth(1) * m_statScale);
    //cout << "----------> ProjectionX 1D Pdf histogram <----------" << endl;
    

    //double xmin = TEvisPdf->GetXaxis()->FindBin(fitTmin);
    //double xmax = TEvisPdf->GetXaxis()->FindBin(fitTmax);
    double xmin = TEvisPdf->GetXaxis()->FindBin(-0.1);
    double xmax = TEvisPdf->GetXaxis()->FindBin(0.1);
    double ymin = TEvisPdf->GetYaxis()->FindBin(Emin);
    double ymax = TEvisPdf->GetYaxis()->FindBin(Emax);

    //m_nNu = TEvisPdf->Integral(xmin, xmax, ymin, ymax, "width");
    m_nNu = TPdf->Integral(xmin, xmax, "width");
    cout << "expected Nu number in this region [ " << xmin << " , " << xmax << "] : " << m_nNu << endl;

    // PDF normalization
    //TPdf->Scale(1./m_nNu);
    //cout << "After scaling totIntegral = " << TPdf->Integral(1, 4001) << endl;

}


void MLEfit::LoadDataset(int evtID)
{
    TFile* dataFile = TFile::Open(dataFileName, "read");
    TTree* tree = (TTree*)dataFile->Get("tFixedStat");
    int m_evtID;
    double m_nuTime1D;
    //double m_nuTime2D;
    //double m_nuEnergy;
    tree->SetBranchAddress("evtID", &m_evtID);
    tree->SetBranchAddress("nuTime1D", &m_nuTime1D);
    //tree->SetBranchAddress("nuTime2D", &m_nuTime2D);
   // tree->SetBranchAddress("nuEnergy", &m_nuEnergy);

    int iNu = 0;
    //cout << "Total SN event in current dataset : " << tree->GetEntries() << endl;
    for(int i=0; i<tree->GetEntries(); i++) {
        tree->GetEntry(i);
        if (m_evtID == evtID) {
            vec_nuTime1D[iNu] = m_nuTime1D;
            //vec_nuTime2D[iNu] = m_nuTime2D;
            //vec_nuEnergy[iNu] = m_nuEnergy;
            iNu++;
        }

        if (m_evtID > evtID) {
            break;
        }
    }
}


// generate dataset from pdf directly ...
void MLEfit::GenerateDataset()
{
    TFile *out = new TFile("nuMass0.0eV_NO_stat1000.root", "recreate");
    TTree* tt  = new TTree("tFixedStat", "tFixedStat");
    int m_evtID;
    double m_nuTime1D;
    tt->Branch("evtID", &m_evtID, "evtID/I");
    tt->Branch("nuTime1D", &m_nuTime1D, "nuTime1D/D");


    double genTmin = -0.015;
    double genTmax = 0.02;
    for (int ievt = 0; ievt<100; ievt++) {
        cout << "Generating Event " << ievt << endl;
        int inu = 0;
        while(inu < m_nevents) {
            double tmpT = TPdf->GetRandom();
            if (tmpT < genTmin or tmpT > genTmax)
                continue;
            else {
                m_evtID = ievt;
                m_nuTime1D = tmpT;
                tt->Fill();
                inu++;
            }
        }
    }
    
    tt->Write();
    out->Close();
}



double MLEfit::GetChi2()
{
    double tmp_nll = 0;

    int nNu = 0;

    for(int iNu=0; iNu<100000; iNu++) {

        if (vec_nuTime1D[iNu] == -10)
            break;
        if (vec_nuTime1D[iNu] < fitTmin or vec_nuTime1D[iNu]>fitTmax)
            continue;
        if (vec_nuTime1D[iNu]+m_deltaT > extTmax or vec_nuTime1D[iNu]+m_deltaT<extTmin) {
            cout << vec_nuTime1D[iNu] << " " << m_deltaT << endl;
            cout << "Warning : data has been shifted outside the fitting range !" << endl;
            continue; 
        }
        double xbin = TPdf->GetXaxis()->FindBin(vec_nuTime1D[iNu]+m_deltaT);

        // Old likelihood
        double prob = TPdf->GetBinContent(xbin); 
        tmp_nll += -TMath::Log(prob * m_scale1D);
        //cout << iNu << " " << m_deltaT << " " << prob << " " << -TMath::Log(prob) << endl;

        nNu++;
    }

    // extended likelihood 
    int xbin1 = TPdf->GetXaxis()->FindBin(fitTmin+m_deltaT);
    int xbin2 = TPdf->GetXaxis()->FindBin(fitTmax+m_deltaT);
    //double InteR = TPdf->Integral(xbin1, xbin2);
    double InteR = TPdf->Integral(xbin1, xbin2, "width");
    m_nNu = InteR * m_scale1D ;
    //cout << "Current region integral = " << InteR << endl;
    //
    tmp_nll += InteR * m_scale1D;

    return tmp_nll;
}


void MLEfit::SetParameters(double *par)
{
    //m_nNu = par[0];
    m_scale1D = par[0];
}



void MLEfit::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}




double MLEfit::GetChiSquare(double maxChi2)
{
    snMinuit = new TMinuit();
    snMinuit->SetFCN(ChisqFCN);
    snMinuit->SetPrintLevel(-1);

    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    snMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    snMinuit->mnparm(iPar, "nsig", 1, 0.01, 0, 1000, ierrflag);   iPar++;

    snMinuit->SetErrorDef(1);
    arglist[0] = 2;
    snMinuit->mnexcm("SET STR", arglist, 1, ierrflag);

    arglist[0] = 5000;
    arglist[1] = 0.1;
    snMinuit->mnexcm("MIGrad", arglist, 1, ierrflag);

    snMinuit->fCstatu.Data();

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    snMinuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

    m_nParameter = snMinuit->GetNumPars();
	for(int i=0; i<m_nParameter; i++)
	{
	    snMinuit->GetParameter(i, m_bestFit[i], m_bestFitError[i]);
	}

    m_DoFit = true;

    //cout << " ====================== " << endl;
    //cout << "minChi2: " << min << " with nSig = " << m_nNu << " and deltaT = " << m_deltaT << endl;
    //cout << " ====================== " << endl;
    m_chi2Min = min;

    delete snMinuit;
    return min;


    return 1;
}




void MLEfit::Plot() {
    TGraph* gg = new TGraph();
    for (int i=0; i<4000; i++) {
        gg->SetPoint(i, TPdf->GetBinCenter(i+1), TPdf->GetBinContent(i+1));
    }


    TCanvas* ca = new TCanvas();
    gg->SetLineColor(kBlue+1);
    gg->Draw("AL");

    ca->SaveAs("TPdf.root");

}









