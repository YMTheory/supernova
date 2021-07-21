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
            - EnergyBinning: NuP: Ethr - 5MeV; Nue: Ethre - 25; IBD: 1.022 - 25+1.022
        - Modify PDFs with different neutrino masses:
            - Garching Model snFluenceDetAtTime
                - No mass original model in SNGarchingIntegFcn(time, E, type)
                - Mass -> snFluenceDetAtTime, different MH has different osc param values.
                    deltaT = 5.14e-3 * (nuMass*nuMass) * (100.0/E/E) * (dist/10)
                    time += time + DeltaT
                    
            - IBD channel: 




    - nllFit.py: configurations in cmd.py -> python cmd.py to see details; 
        - A likelihood fit;
        - Read fitting pdf from etSpec root histogram (Th2D->ProjectX);
        - Read toy data from datasets;
        - loop 100 events
    - datasets: 
        - tFixedStat: nuTime1D->toy dataset of nu time profiles 
    - etSpec:
        - fitting Pdf: hET_xxx, hETBKG_xxx...
