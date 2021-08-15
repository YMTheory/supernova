#include "fitter.hh"

int main()
{
    MLEfit* nll = new MLEfit();
    nll->setPdfFileName("/junofs/users/miaoyu/supernova/wenlj/etSpec/fineSpec/TEvisPDF_mod82503_cha1nue_mh2_mNu0.0eV_10.0kpc_0.1s_Evmax25.root");
    nll->setDataFileName("/junofs/users/miaoyu/supernova/wenlj/dataset/fineSpec/10.0kpc/TEvisDATA_mod82503_cha1nue_mh2_mNu0.0eV_10.0kpc_0.1s_Evmax25_Ethr0.1MeV_group1.root");
    nll->LoadPdf();


    for (int iEvt=1; iEvt<101; iEvt++) {

        nll->LoadDataset(iEvt);

        // coarse scanning 
        int NScan = 35;
        double Tstep = 0.0005;

        double bestFitT = 0;
        double locMinNll = 100;

        for (int iScan = 0; iScan<NScan; iScan++) {
            double deltaT = Tstep*iScan + (0-(NScan-1)/2.)*Tstep;
            nll->setDeltaT(deltaT);
            //nll->GetChi2();
            nll->GetChiSquare(100000);
            if (locMinNll > nll->getChiMin() )  {
                locMinNll = nll->getChiMin();
                bestFitT = nll->getDeltaT();
            }
        }   


        // fine scanning
        NScan = 25;
        Tstep  = 0.0001;
        locMinNll = 1000;
        for (int iScan = 0; iScan<NScan; iScan++) {
            double deltaT = bestFitT + (2*iScan-NScan)*Tstep/2;
            nll->setDeltaT(deltaT);
            //nll->GetChi2();
            nll->GetChiSquare(100000);
            cout << iEvt << " " << deltaT << " " << nll->getChiMin() << " " << nll->getNsig() << endl;
        }   



    }


    return 0;
}
