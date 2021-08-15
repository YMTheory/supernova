#include "MLEfit.hh"

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TCanvas.h>


int MLEfit::m_fitOption = 0;

double MLEfit::vec_nuTime1D[1000] = {-10};
double MLEfit::vec_nuTime2D[1000] = {-10};
double MLEfit::vec_nuEnergy[1000] = {-10};

TH1D* MLEfit::TPdf;
TH2D* MLEfit::TEvisPdf;

double MLEfit::Tmin;
double MLEfit::Tmax;
double MLEfit::Emin;
double MLEfit::Emax;

// pre-setting
double MLEfit::fitTmin = 0.0;
double MLEfit::fitTmax = 0.02;
double MLEfit::extTmin = -0.05;
double MLEfit::extTmax = 0.05;

double MLEfit::m_deltaT = 0;

double MLEfit::m_nNu = 1;

double MLEfit::m_chi2;
double MLEfit::m_chi2Min;
int MLEfit::m_nParameter;
double MLEfit::m_bestFit[1];
double MLEfit::m_bestFitError[1];
bool MLEfit::m_DoFit;



MLEfit::MLEfit()
{
    for(int i=0; i<1000; i++) {
        vec_nuTime1D[i] = -10;
        vec_nuTime2D[i] = -10;
        vec_nuEnergy[i] = -10;
    }

    if (m_fitOption == 0)
        cout << "In this script we adopt 1D fitting " << endl;
    if (m_fitOption == 1)
        cout << "In this script we adopt 2D fitting " << endl;
}

MLEfit::~MLEfit()
{;}

void MLEfit::LoadPdf()
{
    TFile* pdfFile = TFile::Open(pdfFileName, "read");
    TEvisPdf = (TH2D*)pdfFile->Get("hET_mod82503_cha1_mh2");
    cout << "----------> Loading 2D Pdf histogram <----------" << endl;
    TPdf = (TH1D*)TEvisPdf->ProjectionX();
    cout << "----------> ProjectionX 1D Pdf histogram <----------" << endl;

    Tmin = TEvisPdf->GetXaxis()->GetXmin();
    Tmax = TEvisPdf->GetXaxis()->GetXmax();
    Emin = TEvisPdf->GetYaxis()->GetXmin();
    Emax = TEvisPdf->GetYaxis()->GetXmax();

    double xmin = TEvisPdf->GetXaxis()->FindBin(fitTmin);
    double xmax = TEvisPdf->GetXaxis()->FindBin(fitTmax);
    double ymin = TEvisPdf->GetYaxis()->FindBin(Emin);
    double ymax = TEvisPdf->GetYaxis()->FindBin(Emax);

    m_nNu = TEvisPdf->Integral(xmin, xmax, ymin, ymax, "width");
    cout << "PDF Loading Part -> nNu = " << m_nNu << endl;
}


void MLEfit::LoadDataset(int evtID)
{
    TFile* dataFile = TFile::Open(dataFileName, "read");
    TTree* tree = (TTree*)dataFile->Get("tFixedStat");
    int m_evtID;
    double m_nuTime1D;
    double m_nuTime2D;
    double m_nuEnergy;
    tree->SetBranchAddress("evtID", &m_evtID);
    tree->SetBranchAddress("nuTime1D", &m_nuTime1D);
    tree->SetBranchAddress("nuTime2D", &m_nuTime2D);
    tree->SetBranchAddress("nuEnergy", &m_nuEnergy);

    int iNu = 0;
    cout << "Total SN event in current dataset : " << tree->GetEntries() << endl;
    for(int i=0; i<tree->GetEntries(); i++) {
        tree->GetEntry(i);
        if (m_evtID == evtID) {
            vec_nuTime1D[iNu] = m_nuTime1D;
            vec_nuTime2D[iNu] = m_nuTime2D;
            vec_nuEnergy[iNu] = m_nuEnergy;
            iNu++;
        }

        if (m_evtID > evtID) {
            break;
        }
    }
}

void MLEfit::checkPdf()
{
   for (int i=0; i<100; i++) {
    double time = -0.1 + 0.2/100 * i;
    int bin = TPdf->FindBin(time);
    cout << time << " " << bin << " " << TPdf->GetBinContent(bin) << endl;
   } 
}




double MLEfit::GetChi2()
{
    double tmp_nll = 0;

    int nNu = 0;

    // 1D fitting
    if (m_fitOption == 0) {
        for(int iNu=0; iNu<1000; iNu++) {
            
            if (vec_nuTime1D[iNu] == -10)
                continue;
            if (vec_nuTime1D[iNu] < fitTmin or vec_nuTime1D[iNu]>fitTmax)
                continue;
            if (vec_nuTime1D[iNu]+m_deltaT > extTmax or vec_nuTime1D[iNu]+m_deltaT<extTmin) {
                //cout << "Warning : data has been shifted outside the fitting range !" << endl;
                continue; 
            }
            double xbin = TPdf->GetXaxis()->FindBin(vec_nuTime1D[iNu]+m_deltaT);
            double prob = TPdf->GetBinContent(xbin);
            //cout << iNu << " " << vec_nuTime1D[iNu] << " " << vec_nuTime1D[iNu] + m_deltaT << " "
            //     << xbin << " " << prob << endl;

            tmp_nll += -TMath::Log(prob);
            nNu++;
        }
        
    }


    // 2D fitting
    if (m_fitOption == 1) {
        for (int iNu=0; iNu<1000; iNu++) {
            if (vec_nuTime2D[iNu] == -10 and vec_nuEnergy[iNu] == -10)
                continue;
            if (vec_nuTime2D[iNu] < fitTmin or vec_nuTime2D[iNu]>fitTmax)
                continue;
            if (vec_nuTime2D[iNu]+m_deltaT < extTmin or vec_nuTime2D[iNu]+m_deltaT>extTmax) {
                //cout << "Warning : data has been shifted outside the fitting range !" << endl;
                continue; 
            }
            double xbin = TEvisPdf->GetXaxis()->FindBin(vec_nuTime2D[iNu]+m_deltaT);
            double ybin = TEvisPdf->GetYaxis()->FindBin(vec_nuEnergy[iNu]);
            double prob = TEvisPdf->GetBinContent(xbin, ybin);

            tmp_nll += -TMath::Log(prob);
            nNu++;

        }
    }
    
    double non_ext_nll = tmp_nll;

    // extended likelihood 
    tmp_nll += -( nNu * TMath::Log(m_nNu) - m_nNu);
    //cout << "expected nNu =" << m_nNu << ", observed nNu = " << nNu << ", Current nll = " << tmp_nll << " and non-extended nll = " << non_ext_nll << endl;
    
    return tmp_nll;
}


void MLEfit::SetParameters(double *par)
{
    m_nNu = par[0];
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

    snMinuit->mnparm(iPar, "nsig", 30, 0.1, 0, 200, ierrflag);   iPar++;

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


void MLEfit::Plot()
{
    // 1D pdf plots ...
    TFile* ff = new TFile("TPdf.root", "recreate");
    TPdf->Write();
    ff->Close();

}










