#ifndef MLEfit_h
#define MLEfit_h

#include <iostream>

#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TMinuit.h>

using namespace std;

class MLEfit
{
    public:
        MLEfit();
        ~MLEfit();

        void LoadPdf                      ();
        void LoadDataset                  (int evtID);
        double GetChiSquare               ( double maxChi2 = 1000000 );
        static void SetParameters         ( double *par);
        static double GetChi2             ();
        static void Plot                  ();
        static void checkPdf              ();
        //
        double getDeltaT                  ()  { return m_deltaT; }                 
        void setDeltaT                    ( double deltaT ) { m_deltaT = deltaT; }
        double getChiMin                  () { return m_chi2Min;}
        void setPdfFileName               (TString filename)  { pdfFileName = filename;}
        void setDataFileName              (TString filename)  { dataFileName = filename;}
        double getNsig                    () { return m_bestFit[0];}
        

    private:
        TString pdfFileName;
        TString dataFileName;


    private:
        static void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);

        TMinuit* snMinuit;

        static TH1D* TPdf;
        static TH2D* TEvisPdf;

        static double Tmin;
        static double Tmax;
        static double Emin;
        static double Emax;
        static double m_deltaT;

        static double fitTmin;
        static double fitTmax;
        static double extTmin;
        static double extTmax;

        static double m_nNu;

        static int m_fitOption;   // 0 for 1D, 1 for 2D

        static double vec_nuTime1D[1000];
        static double vec_nuTime2D[1000];
        static double vec_nuEnergy[1000];

        static double m_chi2;
        static double m_chi2Min;
        static int m_nParameter;
        static double m_bestFit[1];
        static double m_bestFitError[1];
        static bool m_DoFit;
        

};

#endif
