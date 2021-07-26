# supernova
supernova study in JUNO

## packages from wen

    超新星模拟中的重要参数：模型名称，探测中微子类型，探测反应道，中微子绝对质量，超新星距离。

    -simulation: packages developed by Huiling Li;
        - data: all numerical models collected from different groups: SN, preSN, long cooling...
        - src: source codes for SN detection simulation;   
            - SNsource: distance
                - SNnumGarchingSrc(imodel) source派生类
                - SNGarchingIntegFcn(imodel)
            - SNeffectLS: energy_threshold
                - SNnueLS: 从SNeffectLS派生，描述探测器响应->微分截面和总截面。

    - generatePDF.py: 产生SNmode的时间谱和能谱 "etSpec/fineSpec"
        - configurations: 需要配置给脚本的inputs
            - detection: channel
            - sourece: distance
            - effectLS: energy threshold
        -Histograms definitions: 时间谱和能谱的分bin方法
            - TimeBinning: -0.1->0.1, 0.00005 / bin
            - EnergyBinning: NuP: Ethr - 5MeV; Nue: Ethre - 25; IBD: 1.022 - 25+1.022
        - Modify PDFs with different neutrino masses:
            - Loop in time bins:
                - Loop in Evis bins:
                1. IBD channel: invoke snFluenceDetAtTime -> neutrino time changes due to nuMass
                    Fill 2D-histogram with updated time and E (Enu/Evis);
                2. NuE channel: similar operation with IBD besides slice histograms saved;                 
            - Garching Model snFluenceDetAtTime -> 时间能量微分谱，假设所有类型的中微子拥有相同的质量
                - No mass original model in SNGarchingIntegFcn(time, E, type)
                - Mass -> snFluenceDetAtTime, different MH has different osc param values. Oscillation considered here.
                    - NO:
                    Fnue = sin2theta13 * Fnue + (1-sin2theta13) * Fnux
                    Fnuebar = cos2theta12*cos2theta13*Fnuebar + (1-cos2theta12*cos2theta13)Fnuxbar
                    Fnux= 0.5*(1-sin2theta13)*Fnue + 0.5*(1+sin2theta13)*Fnux
                    Fnuxbar = 0.5*(1-cos2theta12*cos2theta13)*Fnuebar + 0.5*(1+cos2theta12*cos2theta13)*Fnuxbar
                    
                    - IO:
                    Fnue = sin2theta12*cos2theta13 * Fnue + (1-sin2theta12*cos2theta13) * Fnux
                    Fnuebar = sin2theta13*Fnuebar + (1-sin2theta13)Fnuxbar
                    Fnux= 0.5*(1-sin2theta12*cos2theta13)*Fnue + 0.5*(1+sin2theta12*cos2theta13)*Fnux
                    Fnuxbar = 0.5*(1-sin2theta13)*Fnuebar + 0.5*(1+sin2theta13)*Fnuxbar
                    
                    deltaT = 5.14e-3 * (nuMass*nuMass) * (100.0/E/E) * (dist/10)
                    time += time + DeltaT
                    
        - Detector response: 探测器响应的描述，以NuE弹性散射道为例

        - SNDetect.cc" consider cross section, Ev->Evis...  (quenching, resolution)




    - generateDataSets.py: 产生toyMC datasets 
        - Read PDF files from etSpec/fineSpec
        - Get total entries in current PDF: -> Poisson sampling, 500 times subsampling;
        - Load 2D and 1D projection (on time) histograms -> RooDataHist, then generate histPdf -> RootHistPdf
        - randomly sampling :
            - fixed stat. -> fixed nu event number
            - rndm  stat. -> poisson random sampling event number



    *问题：目前的拟合似乎是基于沉积能量Edep所填直方图进行的，而并没有加上quenching和resolution等一系列探测器响应。*


    - nllFit.py: configurations in cmd.py -> python cmd.py to see details; 
        - A likelihood fit;
        - Structure:
            1. load inputHist from etSpec/fineSpec rootFile; inputHist1D: as inputHist projectionX with Ethr and Evismax limits e.g. [0.1MeV, 10.0MeV]   (用于构建pdf)
            2. convert inputHist -> RooDataHist -> RooHistPdf
               inputHist1D -> RooDataHist -> RooHistPdf
               ExtendedPdf  ->  fitPdf1D
            3. load MC datasets from dataset/fineSpec root file -> only get nuTime1D and evt1D here
            4. loop in datasets:
                a. RooDataSet (nuTime1D with current evtID)
                b. fitdata construction: add(timeTmp, weight, weightError) -> 1D fitting
                c. fitTo dataset -> 5 times trys -> sort the minimum minNll as the final results
                d. fitting parameters: nsig, offset ( 2 physical parameters )

            - global time offset requires fine tuning ...

        - Outputs: 
            1. A summary results for "iSub, offset, nll, nsig, nEvtTrue, status"
            2. pdf for fitting plots


    - nllProfile.py

    - nllProfileMH.py



    - datasets: 
        - tFixedStat: nuTime1D->toy dataset of nu time profiles 
    - etSpec:
        - fitting Pdf: hET_xxx, hETBKG_xxx...
