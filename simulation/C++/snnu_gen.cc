/*************************************************************************
 @Author: MiaoYu ---> miaoyu@ihep.ac.cn
 @Created Time : Fri Feb 10 14:27:07 2023
 @File Name: snnu_gen.cc
 ************************************************************************/

#include "TH2D.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"

#include "SNdetect.hh"
#include "SNsource.hh"
#include "SNnumGarchingSrc.hh"
#include "SNGarchingIntegFcn.hh"
//#include "SNnumBurrowsSrc.hh"
#include "SNBurrowsIntegFcn.hh"
#include "SNeffectLS.hh"

int main(int argc, char* argv[]) {
    Int_t imod = -1;//SN neutrino source model
    Int_t icha = -1;//nu interaction channel
    Double_t nuMass = 0.0;// nu mass
    Int_t MH = 0 ; // mass ordering
    Double_t Ethr = 0.2 ; // energy threshold MeV
    Double_t tmin = 0.0;
    Double_t tmax = 0.0;
    Double_t scale = 1;
    Int_t    dist  = 10;
    //==========getopt()===========//
    Int_t opt;
    while((opt = getopt(argc, argv, ":e:m:n:c:t:o:x:y:d:s:"))!= -1){
        switch(opt){
            case'm':
                imod = atoi(optarg);
                break;
            case'c':
                icha = atoi(optarg);
                break;
            case 'n':
                nuMass = atof(optarg);
                break;
            case 'o':
                MH = atoi(optarg);
            case 't':
                Ethr = atof(optarg);
            case 'x':
                tmin = atof(optarg);
            case 'y':
                tmax = atof(optarg);
            case 'd':
                dist = atoi(optarg);
            case 's':
                scale = atof(optarg);
            //case'?':
            //    printf("Error: wrong command input option %s\n",optarg);
            //    break;
            //case':':
            //    printf("Error: Forget to add parameter after option %s\n",optarg);
            //    break;
        }
    }

    std::cout << "\n" 
        << "************************************************** \n"
        << "=========== Configuration: ===========>\n"
        << "-m SN burst neutrino model " << imod << "\n"
        << "-c LS channel (0 for pES, 1 for eES, 2 for IBD, 3 for CEvNS) " << icha<< "\n"
        << "-n neutrino mass (unit: eV) " << nuMass<< "\n"
        << "-o mass ordering (1 for NO, 2 for IO) " << MH<< "\n"
        << "-t threshold (unit: MeV) " << Ethr << "\n"
        << "-x tmin (start time for nu: s) " << tmin << "\n"
        << "-y tmax (end time for nu: s) " << tmax << "\n"
        << "-d distance (unit: kpc) " << dist << "\n"
        << "-s statistics scale factor " << scale << "\n"
        << "**************************************************"
        << "\n"
        << std::endl;

    if(icha == -1){
        std::cout << "interaction channel; -c" << std::endl;
        exit(-1);
    }
    if(imod == -1){
        std::cout << "No SN neutrino model specified; -m" << std::endl;
        exit(-1);
    }

    TString MO[3] = {"NONE", "NO", "IO"};


    //=========================analyze===========================//
    //--------initial SNdetect--------//
    //Int_t icha = 2;
    SNdetect::channelName chaname;
    if(icha == 0)chaname = SNdetect::NuP;
    if(icha == 1)chaname = SNdetect::NuE;
    if(icha == 2)chaname = SNdetect::IBD;
    if(icha == 3)chaname = SNdetect::CEvNS;

    SNdetect* pdet = SNdetect::instance();
    pdet->setSrcModel(imod);
    pdet->getPointerSrc()->setSNDistance(dist);
    pdet->setChannel(chaname);
    Double_t erelres = 3e-2;
    pdet->getPointerEffectLS()->setERelativeRes(erelres);
    pdet->getPointerEffectLS()->setThresholdE(0.);
    pdet->initFCN(); 

    SNsource* mc = pdet->getPointerSrc();
    
    SNeffectLS* peffLS = pdet->getPointerEffectLS();
    double Npar = 0;
    if (chaname == SNdetect::NuP) {
        Npar = peffLS->getNumberOfProton();
    } else if (chaname == SNdetect::NuE) {
        Npar = peffLS->getNumberOfCarbon() * 6 + peffLS->getNumberOfProton();
    } else if (chaname == SNdetect::IBD) {
        Npar = peffLS->getNumberOfProton();
    }

    if(0) {
    std::cout << "-----> fluence = " << pdet->getPointerSrc()->oneSNFluenceDetAtTime(0.02, 60, 1, 1) << "\n" 
              << "-----> total particle number = " << peffLS-> getNumberOfCarbon() << "\n"
              << "-----> differential cross section = " << peffLS->differentialXS(60, 0.02, 0) << " cm^2" << "\n"
              << "-----> getEventAboveEthrVisAtTime = " << pdet->getEventAboveEthrVisAtTime(0.02, 0.15, 1, 1) << "\n"
              << std::endl;
    }

    // Total time window configuration
    Float_t itmin = -0.03;
    Float_t itmax = 0.07;
    if (imod < 80000) // Burrows models
    {
        itmin = 0.;
        itmax = 0.10;
    }
    // time binning
    Double_t step_t         = 0.001;  // unit: s
    Double_t step_t_fine    = step_t / 10.;
    Int_t nbin_it           = int((itmax - itmin)/step_t);
    Int_t nbin_t            = int((tmax - tmin)/step_t);
    Int_t nbin_t_fine       = int((tmax - tmin)/step_t_fine);
    std::cout << "\n " << "====== Time binning info ====== " << std::endl;
    std::cout << "Total time binning: from " << itmin << " to " << itmax << " with total " << nbin_it << " bins." << std::endl;
    std::cout << "Current time binning: from " << tmin << " to " << tmax << " with total " << nbin_t << " bins." << std::endl; 
    std::cout << "Current fine time binning: from " << tmin << " to " << tmax << " with total " << nbin_t_fine << " bins." << std::endl; 
    std::cout << "=====================================" << std::endl;
    
    // energy binning 
    Double_t Evmin = 0;
    Double_t Evmax = 90; 
    Double_t step_Ev = 1;
    Int_t nbins_Ev = (Evmax-Evmin)/step_Ev;
    //Double_t Evismin = Ethr;//0;
    Double_t Evismin = 0;
    Double_t Evismax;
    Double_t step_Evis;
    if(chaname == SNdetect::NuP){
        Evismax = 5.;
        step_Evis = 0.05;
    }
    if(chaname == SNdetect::NuE || chaname == SNdetect::IBD){
        Evismax = 80.;
        step_Evis = 0.2;
    }
    if(chaname == SNdetect::CEvNS){
        Evismax = 0.3;
        step_Evis = 0.001;
    }
    Int_t nbins_Evis = (Evismax-Evismin)/step_Evis;

    std::cout << "\n " << "====== Energy binning info ====== " << std::endl;
    std::cout << "Visible energy binning: from " << Evismin << " to " << Evismax << " with total " 
              << nbins_Evis << " bins and step " << step_Evis << " MeV."
              <<std::endl;
    std::cout << "Neutrino energy binning: from " << Evmin << " to " << Evmax << " with total " 
              << nbins_Ev << " bins and step " << step_Ev << " MeV."
              <<std::endl;
    std::cout << "=====================================" << std::endl;
    

    TString chaName[4] = {"pES","eES","IBD", "CEvNS"};

    bool flag1D     = true;
    bool flag2D     = false;
    bool flag2D_new = false;
    bool flag2D_root = false;

    std::cout << "\n"
              << "========= Output info ==========" << "\n"
              << "\n";

    if (flag1D) {

        //// get the visible events time spectra above Ethr and save into histogram and root file.
        //// 1D histogram
        TString modelName;
        modelName = Form("Garchings%d", imod);
        TString fn = Form("%s_PDF_%s_%dkpc_%s_%.2fMeV_%.3fs-%.3fs_scale%.3f_v2.root",  modelName.Data(), MO[MH].Data(), dist, chaName[icha].Data(), Ethr, tmin, tmax, scale);
        std::cout << "output filename : " << fn << std::endl;
        TFile* f = new TFile(fn, "recreate");
        TH1D* h1 = new TH1D("h1", "visible energy spectrum rate", nbin_t, tmin, tmax);

        double factor = 1;
        double tshift = 0;
        if(imod > 6000 and imod < 7000) { 
            factor = 1e50; 
            tshift = 0.03;
        }
        for (int ipt=0; ipt<nbin_t; ipt++) {
            double t = tmin + (0.5 + ipt) * step_t + tshift;
            std::cout << "Running bin " << ipt << " time " << t << std::endl;
            double flu = 0;
            for(int ii=0; ii<6; ii++) {
                //flu += pdet->getEventAboveEthrObsAtTime(t, Ethr, ii, MH) * factor;  // for 1 ms interval
                flu += pdet->getEventAboveEthrVisAtTime(t, Ethr, ii, MH) * factor;  // for 1 ms interval
                //std::cout << "In main: " << t << " " << ii << " " << flu << std::endl;
            }
            //std::cout << "Filling histogram: " << ipt+1 << " " << flu * scale << std::endl;
            h1->SetBinContent(ipt+1, flu*scale);
        }

        h1->Write();
        f->Close();
    }
    
    if (flag2D) {
        TString modelName;
        modelName = Form("Garching%d", imod);
        TString fn = Form("%s_PDF_%s_%s_%dkpc_%.3fs-%.3fs_nuMass%.1f_scale%.3f_test2D.root",  modelName.Data(), chaName[icha].Data(), MO[MH].Data(), dist,  tmin, tmax, nuMass, scale);
        std::cout << "> Output 2-dimensional PDF filename : " << fn << std::endl;
        
        TFile* f = new TFile(fn, "recreate");
        TH2D* h1 = new TH2D("h1", "time-visible energy 2D hist", nbin_it, itmin, itmax, nbins_Evis, Evismin, Evismax);
        Double_t array[nbin_it][nbins_Evis];
        for(int i =0; i<nbin_it; i++) {
            for (int j=0; j<nbins_Evis; j++) {
                array[i][j] = 0;
            }
        }

        for (int ipt=0; ipt<nbin_t_fine; ipt++) {
            std::cout << ">>> Processing timeBin " << ipt << std::endl;
            double tmp_t;
            for (int iEvis=0; iEvis<nbins_Evis; iEvis++) {
                double EvisTmp = Evismin + (iEvis + 0.5) * step_Evis;

                if  (chaname == SNdetect::NuE or chaname == SNdetect::NuP) {
                    for (int iEv=0; iEv<nbins_Ev; iEv++) {
                        double EvTmp = Evmin + (iEv + 0.5) * step_Ev;
                        double tot_fluence = 0;
                        for (int ii=0; ii<6; ii++) {
                            tmp_t           = tmin + (ipt+0.5) * step_t_fine ; // unit : s
                            double fluence  = mc->snFluenceDetAtTime(tmp_t, nuMass, EvTmp, ii, MH) * step_Evis;
                            double T        = pdet->getTFromEvis(EvisTmp);
                            double dxs      = peffLS->differentialXS(EvTmp, T, ii);
                            double dTdEvis  = (pdet->getTFromEvis(EvisTmp+0.005) - pdet->getTFromEvis(EvisTmp)) / 0.005;
                            tot_fluence    += Npar* fluence * dxs * (step_t_fine / step_t) * dTdEvis;
                            //double dTqdTp   = (pdet->getEvisFromT(T+0.005) - pdet->getEvisFromT(T)) / 0.005;
                            //tot_fluence    += Npar* fluence * dxs * (step_t_fine / step_t) / dTqdTp;
                            if(0) {
                                std::cout << "======= Details check ======= \n" 
                                          << "Flavour " << ii << ", shifted time = " << tmp_t << " s, "
                                          << "Evis = " << EvisTmp << " MeV, "
                                          << "T = " << T << " MeV, "
                                          << "Enu = " << EvTmp << " MeV, "
                                          << "fluence = " << fluence << ", "
                                          << "dxs = " << dxs << " cm2, "
                                          << "total fluence = " << tot_fluence 
                                          << std::endl;
                            }
                        }
                        int idx = (tmp_t - itmin) / step_t;
                        //std::cout << tmp_t << " " << itmin << " " << step_t << " " << idx << std::endl;
                        int idy = EvisTmp / step_Evis;
                        if (idx >=0 and idx < nbin_it and tot_fluence > 0)
                            //std::cout << "** Total fluence -> " << idx << " " << idy << " " << tot_fluence << std::endl;
                            array[idx][idy] += tot_fluence;
                    }
                }
            }
        }

        for(int i =0; i<nbin_it; i++) {
            for (int j=0; j<nbins_Evis; j++) {
                //std::cout << i << " " << j << " " << array[i][j] << std::endl;
                h1->SetBinContent(i+1, j+1, array[i][j]);
            }
        }
        h1->Write();
        f->Close();
    }


    if(flag2D_new) {

        TString modelName;
        modelName = Form("Garching%d", imod);
        TString fn = Form("%s_PDF_%s_%s_%dkpc_nuMass%.1f_scale%.3f_test2Dnew.root",  modelName.Data(), chaName[icha].Data(), MO[MH].Data(), dist, nuMass, scale);
        std::cout << "> Output 2-dimensional PDF filename from new alg: " << fn << std::endl;

        TFile* f = new TFile(fn, "recreate");
        TH2D* h1 = new TH2D("h1", "time-visible energy 2D hist", nbin_it, itmin, itmax, nbins_Evis, Evismin, Evismax);
        Double_t array[nbin_it][nbins_Evis];

        for (int ipt=0; ipt<nbin_it; ipt++) {
            std::cout << ">>> Processing time bin " << ipt << " in total " << nbin_it << " bins."<< std::endl;
            double TTmp = itmin + (ipt + 0.5) * step_t;

            for (int iEvis=0; iEvis < nbins_Evis; iEvis++) {
                double EvisTmp = Evismin + (iEvis + 0.5) * step_Evis;

                double tot_fluence = 0;
                if  (chaname == SNdetect::NuE or chaname == SNdetect::NuP) {
                    for (int iEv=0; iEv<nbins_Ev; iEv++) {
                        double EvTmp = Evmin + (iEv + 0.5) * step_Ev;
                        double deltaT = 5.14e-3 * nuMass*nuMass * (100/EvTmp/EvTmp) * (dist/10.); // unit: s

                        for (int f=0; f<6; f++) {
                            double fluence  = mc->oneSNFluenceDetAtTime(TTmp-deltaT, EvTmp, f, MH);
                            double T        = pdet->getTFromEvis(EvisTmp);
                            double dxs      = peffLS->differentialXS(EvTmp, T, f);
                            double dTdEvis  = (pdet->getTFromEvis(EvisTmp+0.005) - pdet->getTFromEvis(EvisTmp)) / 0.005;
                            tot_fluence    += Npar * fluence * dxs * dTdEvis * step_Ev;
                        }
                    }
                }

                if (chaname == SNdetect::IBD) {
                    double EvTmp =  EvisTmp + 1.8;
                    double deltaT = 5.14e-3 * nuMass*nuMass * (100/EvTmp/EvTmp) * (dist / 10.);
                    int f = 1;
                    double fluence = mc->oneSNFluenceDetAtTime(TTmp-deltaT, EvTmp, f, MH);
                    double xs = peffLS->totalXS(EvTmp, f);
                    tot_fluence += Npar * fluence * xs ;
                }

                array[ipt][iEvis] = tot_fluence;
            }
        }

        for(int i =0; i<nbin_it; i++) {
            for (int j=0; j<nbins_Evis; j++) {
                //std::cout << i << " " << j << " " << array[i][j] << std::endl;
                h1->SetBinContent(i+1, j+1, array[i][j]);
            }
        }
        h1->Write();
        f->Close();

    }

    if (flag2D_root) {
        
        TString modelName;
        modelName = Form("Garching%d", imod);
        TString fn = Form("%s_PDF_%s_%s_%dkpc_nuMass%.1feV_scale%.3f_tmin%.3fstmax%.3fs_response_2Droot.root",  modelName.Data(), chaName[icha].Data(), MO[MH].Data(), dist, nuMass, scale, tmin, tmax);
        std::cout << "> Output 2-dimensional PDF filename from new alg: " << fn << std::endl;

        TFile* f = new TFile(fn, "recreate");
        TH2D* h1 = new TH2D("h1", "time-visible energy 2D hist", nbin_it, itmin, itmax, nbins_Evis, Evismin, Evismax);
        Double_t array[nbin_it][nbins_Evis];

        // Initialisation
        for (int ipt=0; ipt<nbin_it; ipt++) {
            for (int iEvis=0; iEvis<nbins_Evis; iEvis++) {
                array[ipt][iEvis] = 0;
            }
        }
    
        for (int ipt=0; ipt<nbin_it; ipt++) {
            double TTmp = itmin + (ipt + 0.5) * step_t;
            if (TTmp < tmin or TTmp >= tmax)
                continue;
            std::cout << ">>> Processing time bin " << ipt << " at " << TTmp << " s." << std::endl;
            for (int iEvis=0; iEvis<nbins_Evis; iEvis++) {
                double EvisTmp = Evismin + (iEvis + 0.5) * step_Evis;
                //std::cout << ">>>>>> Processing Eobs bin " << iEvis << " at " << EvisTmp << " MeV." << std::endl;
                array[ipt][iEvis] = pdet->getEobsSpectrumAtTimeWithMass(TTmp, EvisTmp, -1, MH, nuMass);
            }
        }
        for(int i =0; i<nbin_it; i++) {
            for (int j=0; j<nbins_Evis; j++) {
                h1->SetBinContent(i+1, j+1, array[i][j]);
                if(0) {
                    std::cout << "Bin content of (" << i+1 << ", " <<  j+1 << ") bin = " << array[i][j] << std::endl;
                }
            }
        }
        h1->Write();
        f->Close();
    }

    
    // test line
    // double Evis = 10;
    // double Eobs = 10;
    // double tt   = 0.02;
    // double Nvis0 = pdet->getEvisSpectrumAtTime(tt, Evis, -1, 1);
    // double Nvis1 = pdet->getEvisSpectrumAtTimeWithMass(tt, Evis, -1, 1, 0.0);
    // double Nvis2 = pdet->getEvisSpectrumAtTimeWithMass(tt, Evis, -1, 1, 1.0);
    // double Nobs0 = pdet->getEobsSpectrumAtTime(tt, Eobs, -1, 1);
    // double Nobs1 = pdet->getEobsSpectrumAtTimeWithMass(tt, Eobs, -1, 1, 0.0);
    // double Nobs2 = pdet->getEobsSpectrumAtTimeWithMass(tt, Eobs, -1, 1, 1.0);
    // 
    // std::cout << "In normal ordering case: \n" 
    //           << "====> w/o energy response: " << "\n"
    //           << "noMass Vis -> " << Nvis0 << "\n"
    //           << "nuMass=0.0 -> " << Nvis1 << "\n"
    //           << "noMass=1.0 -> " << Nvis2 << "\n"
    //           << "====> w/ energy response: " << "\n"
    //           << "noMass Vis -> " << Nobs0 << "\n"
    //           << "nuMass=0.0 -> " << Nobs1 << "\n"
    //           << "noMass=1.0 -> " << Nobs2 << "\n"
    //           << std::endl;



}
