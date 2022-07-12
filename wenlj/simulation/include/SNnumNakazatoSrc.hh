#ifndef _SNnumNakazatoSrc_hh
#define _SNnumNakazatoSrc_hh



class SNsource;
class SNNakazatoIntegFcn;
class TGraph;
class TF1;

class SNnumNakazatoSrc : public SNsource{
    public:
        SNnumNakazatoSrc();
        SNnumNakazatoSrc(int imode);
        ~SNnumNakazatoSrc();

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
        SNNakazatoIntegFcn* pgarfcn;
        TF1* fcnTimeInterval;
        //TGraph* grSpecTimeInterv[3]; 
        //double timeInterval[2];
    public:
        ClassDef(SNnumNakazatoSrc, 1);
};



#endif

