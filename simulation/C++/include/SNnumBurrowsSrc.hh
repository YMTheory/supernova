#ifndef _SNnumBurrowsSrc_hh
#define _SNnumBurrowsSrc_hh



class SNsource;
class SNBurrowsIntegFcn;
class TGraph;
class TF1;

class SNnumBurrowsSrc : public SNsource{
    public:
        SNnumBurrowsSrc();
        SNnumBurrowsSrc(int imode);
        ~SNnumBurrowsSrc();

        virtual double oneSNFluenceDet(double E, int type);
        virtual double oneSNFluenceDet(double E, int type, int MH);
        virtual double totalSNFluenceDet(double E);
        virtual double totalSNFluenceDet(double E, int MH);
        virtual double totalSNFluenceDetTimeIntegE(double time, int MH);
        virtual double oneSNFluenceDetTimeIntegE(double time, int type, int MH); 
        
        virtual double totalSNFluenceDetAtTime(double time, double E, int MH);
        virtual double oneSNFluenceDetAtTime(double time, double E, int type, int MH);
        double snFluenceDetAtTime(double &time, double nuMass, double E, int type, int MH);

        virtual double oneSNFluenceDetTimeInterval(double E, double tfirst,double tlast, int type, int MH);
        virtual double totalSNFluenceDetTimeInterval(double E, double tfirst, double tlast, int MH);
        virtual double oneSNFluenceDetTimeInterval(double E, double tfirst,double tlast, int type);
        virtual double totalSNFluenceDetTimeInterval(double E, double tfirst, double tlast);

        virtual double oneSNLuminosityTime(double time, int type);
        virtual double oneSNAverageETime(double time, int type);

        virtual void getTimeRange(double& tmin, double& tmax, int type);

    private:
        void readSpectrum(int imode);
        TGraph* grspectrum[3];
        SNBurrowsIntegFcn* pgarfcn;
        TF1* fcnTimeInterval;

    //public:
    //    ClassDef(SNnumBurrowsSrc, 1);
};



#endif
