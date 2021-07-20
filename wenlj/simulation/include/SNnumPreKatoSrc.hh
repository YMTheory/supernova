#ifndef _SNnumPreKatoSrc_hh
#define _SNnumPreKatoSrc_hh


class SNsource;
class TGraph;
class TGraph2D;
class SNpreKatoIntegFcn;
class TF1;
class SNnumPreKatoSrc: public SNsource{
    public:
        SNnumPreKatoSrc();
        SNnumPreKatoSrc(int imode);//models in [http://asphwww.ph.noda.tus.ac.jp/snn/] used
        ~SNnumPreKatoSrc();

        virtual double oneSNFluenceDet(double E, int type);
        virtual double oneSNFluenceDet(double E, int type, int MH);
        virtual double totalSNFluenceDet(double E);
        virtual double totalSNFluenceDet(double E, int MH);
        virtual double totalSNFluenceDetTimeIntegE(double time, int MH);
        virtual double oneSNFluenceDetTimeIntegE(double time, int type, int MH); 
        

        virtual double totalSNFluenceDetAtTime(double time, double E, int MH);
        virtual double oneSNFluenceDetAtTime(double time, double E, int type, int MH);
        
        virtual double oneSNFluenceDetTimeInterval(double E, double tfirst,double tlast, int type, int MH);
        virtual double totalSNFluenceDetTimeInterval(double E, double tfirst, double tlast, int MH);
        virtual double oneSNFluenceDetTimeInterval(double E, double tfirst,double tlast, int type);
        virtual double totalSNFluenceDetTimeInterval(double E, double tfirst, double tlast);
        
        virtual double oneSNLuminosityTime(double time, int type);
        virtual double oneSNAverageETime(double time, int type);
        
        virtual void getTimeRange(double& tmin, double& tmax, int type);
    public:
        SNpreKatoIntegFcn* ppreInteg;
        TF1* f1TE;
        ClassDef(SNnumPreKatoSrc, 1);
};

#endif
