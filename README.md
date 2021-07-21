# supernova
supernova study in JUNO

## packages from wen

    超新星模拟中的重要参数：模型名称，探测中微子类型，探测反应道，中微子绝对质量，超新星距离。

    -simulation: packages developed by Huiling Li;
        - data: all numerical models collected from different groups: SN, preSN, long cooling...
        - src: source codes for SN detection simulation;   
            - SNnumGarchingSrc(imodel)
            - SNGarchingIntegFcn(imodel)
            - SNsource: distance
            - SNeffectLS: energy_threshold, 

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
                    



    - nllFit.py: configurations in cmd.py -> python cmd.py to see details; 
        - A likelihood fit;
        - Read fitting pdf from etSpec root histogram (Th2D->ProjectX);
        - Read toy data from datasets;
        - loop 100 events
    - datasets: 
        - tFixedStat: nuTime1D->toy dataset of nu time profiles 
    - etSpec:
        - fitting Pdf: hET_xxx, hETBKG_xxx...
