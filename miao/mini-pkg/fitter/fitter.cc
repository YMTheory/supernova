#include "fitter.hh"
#include <string>
#include <TMath.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <sstream>


int main(int argc, char* argv[])
{
    std::string fitMass   = argv[1];
    std::string fitMH     = argv[2];
    std::string dataMH    = argv[3];
    std::string dist      = argv[4];

    bool doFit = false;
    bool genData = true;

    std::string MO;
    if (fitMH == "1") 
        MO = "NO";
    if (fitMH == "2")
        MO = "IO";

    double scale =  (10./strtof(dist.c_str(), NULL)) * (10./strtof(dist.c_str(), NULL));

    std::string group;
    if (dist == "10.0")
        group = "1";
    if (dist == "5.0")
        group = "999";
    if (dist == "2.0")
        group = "880" ;
    if (dist == "1.0")
        group = "770" ;

    float nuMass = strtof(fitMass.c_str(), NULL);
    float d = strtof(dist.c_str(), NULL);

    scale = 10000.;

    MLEfit* nll = new MLEfit();
    //nll->setDataFileName("/junofs/users/miaoyu/supernova/wenlj/dataset/fineSpec/10.0kpc/TEvisDATA_mod82503_cha1nue_nuType-1_mh"+dataMH+"_mNu0.0eV_10.0kpc_0.1s_Evmax25_Ethr0.1MeV_group"+group+".root");
    //nll->setDataFileName("/junofs/users/miaoyu/supernova/wenlj/dataset/fineSpec/10.0kpc/TEvisDATA_mod82503_cha1nue_nuType-1_mh"+dataMH+"_mNu0.0eV_10.0kpc_0.1s_Evmax25_Ethr0.1MeV_group"+group+".root");
    nll->setDataFileName("nuMass0.0eV_IO_stat10000_wen.root");
    //nll->setDataFileName("nuMass0.0eV_IO_stat10000.root");
    nll->setStatScale(scale);
    cout << "scale factor = " << scale << " , group id " << group << endl;

    nll->setPdfFileName("/junofs/users/miaoyu/supernova/wenlj/etSpec/fineSpec/10kpc/"+MO+"/TEvisPDF_mod82503_cha1nue_nuType-1_mh"+fitMH+"_mNu"+fitMass+"eV_10.0kpc_0.1s_Evismax25.00_Sum.root");
    //nll->setPdfFileName("/junofs/users/miaoyu/supernova/wenlj/etSpec/fineSpec/TEvisPDF_mod82503_cha1nue_mh2_mNu"+fitMass+"eV_10.0kpc_0.1s_Evmax25.root");
    cout << "PdfFileName: " << nll->getPdfFileName() << " with dataFileName: "<< nll->getDataFileName() << endl;
    nll->LoadPdf();

    if (genData)
        nll->GenerateDataset();

    if (doFit) {

        for (int iEvt=1; iEvt<2; iEvt++) {

            nll->LoadDataset(iEvt);

            // coarse scanning 
            int NScan = 35;
            double Tstep = 0.0005;

            double locMinNll = 100000;
            double dT = 5.14e-3 * (nuMass * nuMass) * (100.0/9/9) * (10./10.);
            //double dT = 5.14e-3 * (nuMass*nuMass) * (100.0/9/9) * (d/10.);
            double bestFitT = 0;

            for (int iScan = 0; iScan<NScan; iScan++) {
                //double deltaT = dT + (2*iScan-NScan)*Tstep/2;
                double deltaT = dT - (NScan-1)/2.*Tstep + Tstep*iScan;
                nll->setDeltaT(deltaT);
                nll->GetChiSquare(100000);
                double cur_chi2 = nll->GetChi2();
                //cout << " coarse " << iEvt << " "  <<  deltaT << " " << cur_chi2 << " " << nll->getNsig() << endl;
                if (locMinNll > cur_chi2 )  {
                    locMinNll = cur_chi2;
                    bestFitT = nll->getDeltaT();
                    //    //cout << "CoarseScanning " << bestFitT << locMinNll << endl;
                }
            }   

            //cout << "bestFitT " << bestFitT << endl;

            // fine scanning
            NScan = 25;
            Tstep  = 0.0001;
            locMinNll = 1000;
            double fineBestT;
            for (int iScan = 0; iScan<NScan; iScan++) {
                //double deltaT = bestFitT + (2*iScan-NScan)*Tstep/2;
                double deltaT = bestFitT - (NScan-1)/2.*Tstep + Tstep*iScan;
                nll->setDeltaT(deltaT);
                //    //nll->GetChi2();
                nll->GetChiSquare(100000);
                double tmp_chi2 = nll->GetChi2();
                cout << iEvt  << " " <<  deltaT << " " << 2*tmp_chi2 << " " << nll->getNsig() << " " << nll->getScale1D() << endl;
                if (locMinNll > tmp_chi2) {
                    locMinNll = tmp_chi2;
                    fineBestT = deltaT;
                }
            }   

            cout << "fineBestT " << iEvt << " " << fineBestT << " " << locMinNll << endl;

        }


    }


    return 0;
}
