import ROOT
import uproot as up
import numpy as np
import matplotlib.pyplot as plt

import sys

if __name__ == "__main__" :

    mod = int(sys.argv[1])
    mod_arr = []
    mod_arr.append(mod)

    # Load P.D.F:

    pdf_NO = []
    pdf_IO = []

    #mod_arr = [81120,81121,81122,81123,82500,82501,82502,82503,82700,82701,82702,82703,84000,84001,84002,84003,91120,91121,91122,91123,92500,92501,92502, 92503,92700,92701,92702,92703,94000,94001,94002,94003]
    Ethr, fitEmax = 0.1, 60
    binEthr, binEmax = 1, 300

    imod = 81120
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile0   = ROOT.TFile(filename, 'read')
    inputHist0 = pdfFile0.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist01D = inputHist0.ProjectionX('_pdfsig0px', binEthr, binEmax, '')
    pdf_NO.append(inputHist01D)
    inputHist90 = pdfFile0.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist901D = inputHist90.ProjectionX('_pdfsig90px', binEthr, binEmax, '')
    pdf_IO.append(inputHist901D)
    imod = 81121
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile1   = ROOT.TFile(filename, 'read')
    inputHist1 = pdfFile1.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist11D = inputHist1.ProjectionX('_pdfsig1px', binEthr, binEmax, '')
    pdf_NO.append(inputHist11D)
    inputHist91 = pdfFile1.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist911D = inputHist91.ProjectionX('_pdfsig91px', binEthr, binEmax, '')
    pdf_IO.append(inputHist911D)
    imod = 81122
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile2   = ROOT.TFile(filename, 'read')
    inputHist2 = pdfFile2.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist21D = inputHist2.ProjectionX('_pdfsig2px', binEthr, binEmax, '')
    pdf_NO.append(inputHist21D)
    inputHist92 = pdfFile2.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist921D = inputHist92.ProjectionX('_pdfsig92px', binEthr, binEmax, '')
    pdf_IO.append(inputHist921D)
    imod = 81123
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile3   = ROOT.TFile(filename, 'read')
    inputHist3 = pdfFile3.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist31D = inputHist3.ProjectionX('_pdfsig3px', binEthr, binEmax, '')
    pdf_NO.append(inputHist31D)
    inputHist93 = pdfFile3.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist931D = inputHist93.ProjectionX('_pdfsig93px', binEthr, binEmax, '')
    pdf_IO.append(inputHist931D)
    imod = 82500
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile4   = ROOT.TFile(filename, 'read')
    inputHist4 = pdfFile4.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist41D = inputHist4.ProjectionX('_pdfsig4px', binEthr, binEmax, '')
    pdf_NO.append(inputHist41D)
    inputHist94 = pdfFile4.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist941D = inputHist94.ProjectionX('_pdfsig94px', binEthr, binEmax, '')
    pdf_IO.append(inputHist941D)
    imod = 82501
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile5   = ROOT.TFile(filename, 'read')
    inputHist5 = pdfFile5.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist51D = inputHist5.ProjectionX('_pdfsig5px', binEthr, binEmax, '')
    pdf_NO.append(inputHist51D)
    inputHist95= pdfFile5.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist951D = inputHist95.ProjectionX('_pdfsig95px', binEthr, binEmax, '')
    pdf_IO.append(inputHist951D)
    imod = 82502
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile6   = ROOT.TFile(filename, 'read')
    inputHist6 = pdfFile6.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist61D = inputHist6.ProjectionX('_pdfsig6px', binEthr, binEmax, '')
    pdf_NO.append(inputHist61D)
    inputHist96 = pdfFile6.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist961D = inputHist96.ProjectionX('_pdfsig96px', binEthr, binEmax, '')
    pdf_IO.append(inputHist961D)
    imod = 82503
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile7   = ROOT.TFile(filename, 'read')
    inputHist7 = pdfFile7.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist71D = inputHist7.ProjectionX('_pdfsig7px', binEthr, binEmax, '')
    pdf_NO.append(inputHist71D)
    inputHist97 = pdfFile7.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist971D = inputHist97.ProjectionX('_pdfsig97px', binEthr, binEmax, '')
    pdf_IO.append(inputHist971D)
    imod = 82700
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile8   = ROOT.TFile(filename, 'read')
    inputHist8 = pdfFile8.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist81D = inputHist8.ProjectionX('_pdfsig8px', binEthr, binEmax, '')
    pdf_NO.append(inputHist81D)
    inputHist98 = pdfFile8.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist981D = inputHist98.ProjectionX('_pdfsig98px', binEthr, binEmax, '')
    pdf_IO.append(inputHist981D)
    imod = 82701
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile9   = ROOT.TFile(filename, 'read')
    inputHist9 = pdfFile9.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist91D = inputHist9.ProjectionX('_pdfsig9px', binEthr, binEmax, '')
    pdf_NO.append(inputHist91D)
    inputHist99 = pdfFile9.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist991D = inputHist99.ProjectionX('_pdfsig99px', binEthr, binEmax, '')
    pdf_IO.append(inputHist991D)
    imod = 82702
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile10   = ROOT.TFile(filename, 'read')
    inputHist10 = pdfFile10.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist101D = inputHist10.ProjectionX('_pdfsig10px', binEthr, binEmax, '')
    pdf_NO.append(inputHist101D)
    inputHist910 = pdfFile10.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9101D = inputHist910.ProjectionX('_pdfsig910px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9101D)
    imod = 82703
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile11   = ROOT.TFile(filename, 'read')
    inputHist11 = pdfFile11.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist111D = inputHist11.ProjectionX('_pdfsig11px', binEthr, binEmax, '')
    pdf_NO.append(inputHist111D)
    inputHist911 = pdfFile11.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9111D = inputHist911.ProjectionX('_pdfsig911px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9111D)
    imod = 84000
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile12   = ROOT.TFile(filename, 'read')
    inputHist12 = pdfFile12.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist121D = inputHist12.ProjectionX('_pdfsig12px', binEthr, binEmax, '')
    pdf_NO.append(inputHist121D)
    inputHist912 = pdfFile12.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9121D = inputHist912.ProjectionX('_pdfsig912px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9121D)
    imod = 84001
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile13   = ROOT.TFile(filename, 'read')
    inputHist13 = pdfFile13.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist131D = inputHist13.ProjectionX('_pdfsig13px', binEthr, binEmax, '')
    pdf_NO.append(inputHist131D)
    inputHist913 = pdfFile13.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9131D = inputHist913.ProjectionX('_pdfsig913px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9131D)
    imod = 84002
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile14   = ROOT.TFile(filename, 'read')
    inputHist14 = pdfFile14.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist141D = inputHist14.ProjectionX('_pdfsig14px', binEthr, binEmax, '')
    pdf_NO.append(inputHist141D)
    inputHist914 = pdfFile14.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9141D = inputHist914.ProjectionX('_pdfsig914px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9141D)
    imod = 84003
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile15   = ROOT.TFile(filename, 'read')
    inputHist15 = pdfFile15.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist151D = inputHist15.ProjectionX('_pdfsig15px', binEthr, binEmax, '')
    pdf_NO.append(inputHist151D)
    inputHist915 = pdfFile15.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9151D = inputHist915.ProjectionX('_pdfsig915px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9151D)
    imod = 91120
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile16   = ROOT.TFile(filename, 'read')
    inputHist16 = pdfFile16.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist161D = inputHist16.ProjectionX('_pdfsig16px', binEthr, binEmax, '')
    pdf_NO.append(inputHist161D)
    inputHist916 = pdfFile16.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9161D = inputHist916.ProjectionX('_pdfsig916px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9161D)
    imod = 91121
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile17   = ROOT.TFile(filename, 'read')
    inputHist17 = pdfFile17.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist171D = inputHist17.ProjectionX('_pdfsig17px', binEthr, binEmax, '')
    pdf_NO.append(inputHist171D)
    inputHist917 = pdfFile17.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9171D = inputHist917.ProjectionX('_pdfsig917px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9171D)
    imod = 91122
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile18   = ROOT.TFile(filename, 'read')
    inputHist18 = pdfFile18.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist181D = inputHist18.ProjectionX('_pdfsig18px', binEthr, binEmax, '')
    pdf_NO.append(inputHist181D)
    inputHist918 = pdfFile18.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9181D = inputHist918.ProjectionX('_pdfsig918px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9181D)
    imod = 91123
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile19   = ROOT.TFile(filename, 'read')
    inputHist19 = pdfFile19.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist191D = inputHist19.ProjectionX('_pdfsig19px', binEthr, binEmax, '')
    pdf_NO.append(inputHist191D)
    inputHist919 = pdfFile19.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9191D = inputHist919.ProjectionX('_pdfsig919px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9191D)
    imod = 92500
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile20   = ROOT.TFile(filename, 'read')
    inputHist20 = pdfFile20.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist201D = inputHist20.ProjectionX('_pdfsig20px', binEthr, binEmax, '')
    pdf_NO.append(inputHist201D)
    inputHist920 = pdfFile20.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9201D = inputHist920.ProjectionX('_pdfsig920px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9201D)
    imod = 92501
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile21   = ROOT.TFile(filename, 'read')
    inputHist21 = pdfFile21.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist211D = inputHist21.ProjectionX('_pdfsig21px', binEthr, binEmax, '')
    pdf_NO.append(inputHist211D)
    inputHist921 = pdfFile21.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9211D = inputHist921.ProjectionX('_pdfsig921px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9211D)
    imod = 92502
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile22   = ROOT.TFile(filename, 'read')
    inputHist22 = pdfFile22.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist221D = inputHist22.ProjectionX('_pdfsig22px', binEthr, binEmax, '')
    pdf_NO.append(inputHist221D)
    inputHist922 = pdfFile22.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9221D = inputHist922.ProjectionX('_pdfsig922px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9221D)
    imod = 92503
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile23   = ROOT.TFile(filename, 'read')
    inputHist23 = pdfFile23.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist231D = inputHist23.ProjectionX('_pdfsig23px', binEthr, binEmax, '')
    pdf_NO.append(inputHist231D)
    inputHist923 = pdfFile23.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9231D = inputHist923.ProjectionX('_pdfsig923px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9231D)
    imod = 92700
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile24   = ROOT.TFile(filename, 'read')
    inputHist24 = pdfFile24.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist241D = inputHist24.ProjectionX('_pdfsig24px', binEthr, binEmax, '')
    pdf_NO.append(inputHist241D)
    inputHist924 = pdfFile24.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9241D = inputHist924.ProjectionX('_pdfsig924px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9241D)
    imod = 92701
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile25   = ROOT.TFile(filename, 'read')
    inputHist25 = pdfFile25.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist251D = inputHist25.ProjectionX('_pdfsig25px', binEthr, binEmax, '')
    pdf_NO.append(inputHist251D)
    inputHist925 = pdfFile25.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9251D = inputHist925.ProjectionX('_pdfsig925px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9251D)
    imod = 92702
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile26   = ROOT.TFile(filename, 'read')
    inputHist26 = pdfFile26.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist261D = inputHist26.ProjectionX('_pdfsig26px', binEthr, binEmax, '')
    pdf_NO.append(inputHist261D)
    inputHist926 = pdfFile26.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9261D = inputHist926.ProjectionX('_pdfsig926px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9261D)
    imod = 92703
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile27   = ROOT.TFile(filename, 'read')
    inputHist27 = pdfFile27.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist271D = inputHist27.ProjectionX('_pdfsig27px', binEthr, binEmax, '')
    pdf_NO.append(inputHist271D)
    inputHist927 = pdfFile27.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9271D = inputHist927.ProjectionX('_pdfsig927px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9271D)
    imod = 94000
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile28   = ROOT.TFile(filename, 'read')
    inputHist28 = pdfFile28.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist281D = inputHist28.ProjectionX('_pdfsig28px', binEthr, binEmax, '')
    pdf_NO.append(inputHist281D)
    inputHist928 = pdfFile28.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9281D = inputHist928.ProjectionX('_pdfsig928px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9281D)
    imod = 94001
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile29   = ROOT.TFile(filename, 'read')
    inputHist29 = pdfFile29.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist291D = inputHist29.ProjectionX('_pdfsig29px', binEthr, binEmax, '')
    pdf_NO.append(inputHist291D)
    inputHist929 = pdfFile29.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9291D = inputHist929.ProjectionX('_pdfsig929px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9291D)
    imod = 94002
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile30   = ROOT.TFile(filename, 'read')
    inputHist30 = pdfFile30.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist301D = inputHist30.ProjectionX('_pdfsig30px', binEthr, binEmax, '')
    pdf_NO.append(inputHist301D)
    inputHist930 = pdfFile30.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9301D = inputHist930.ProjectionX('_pdfsig930px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9301D)
    imod = 94003
    filename = "/junofs/users/miaoyu/supernova/production/PDFs/2kpc/evistSpec_mod%d_cha1.root"%(imod)
    print(filename)
    pdfFile31   = ROOT.TFile(filename, 'read')
    inputHist31 = pdfFile31.Get("hET_mod%d_cha1_mh%d"%(imod, 1))
    inputHist311D = inputHist31.ProjectionX('_pdfsig31px', binEthr, binEmax, '')
    pdf_NO.append(inputHist311D)
    inputHist931 = pdfFile31.Get("hET_mod%d_cha1_mh%d"%(imod, 2))
    inputHist9311D = inputHist931.ProjectionX('_pdfsig931px', binEthr, binEmax, '')
    pdf_IO.append(inputHist9311D)

    
    lowBin, highBin = inputHist01D.GetBinLowEdge(1), inputHist01D.GetBinLowEdge(120)+inputHist01D.GetBinWidth(120)
    print(lowBin, highBin)


    # Load Data
    #fig, ax = plt.subplots()

    nEvent = 500

    testNO, testIO = [], []
    for imod in mod_arr:
        filename = "/junofs/users/miaoyu/supernova/production/Data/10kpc/data_mod%d_cha1_mh1_val.root"%(imod)
        print(filename)
        ff = up.open(filename)
        evtID = ff["evt"]["evtID"].array()
        time = ff['evt']['time'].array()
    
        for n in range(nEvent):
            h_data = ROOT.TH1D("h_data", "", 120, -0.02, 0.1)
            tmp_NO, tmp_IO = [], []
            for i, j in zip(evtID, time):
                if i == n:
                    h_data.Fill(j)
                if i > n:
                    break


            # KS TEST ...
            for ii in range(32):
                h_pdf = pdf_NO[ii]
                print(h_pdf)
                tmp_NO.append(h_data.KolmogorovTest(h_pdf))
                h_pdf = pdf_IO[ii]
                tmp_IO.append(h_data.KolmogorovTest(h_pdf))

            
            tmp_NO = np.array(tmp_NO)
            tmp_IO = np.array(tmp_IO)

            aveNO = np.mean(tmp_NO)
            aveIO = np.mean(tmp_IO)

            testNO.append(aveNO - aveIO)


    #ax.hist(testNO, bins=100, color="blue", edgecolor="black", label="NO")


    for imod in mod_arr:
        filename = "/junofs/users/miaoyu/supernova/production/Data/10kpc/data_mod%d_cha1_mh2_val.root"%(imod)
        print(filename)
        ff = up.open(filename)
        evtID = ff["evt"]["evtID"].array()
        time = ff['evt']['time'].array()
    
        for n in range(nEvent):
            h_data = ROOT.TH1D("h_data", "", 120, -0.02, 0.1)
            tmp_NO, tmp_IO = [], []
            for i, j in zip(evtID, time):
                if i == n:
                    h_data.Fill(j)
                if i > n:
                    break


            # KS TEST ...
            for ii in range(32):
                h_pdf = pdf_NO[ii]
                print(h_pdf)
                tmp_NO.append(h_data.KolmogorovTest(h_pdf))
                h_pdf = pdf_IO[ii]
                tmp_IO.append(h_data.KolmogorovTest(h_pdf))

            
            tmp_NO = np.array(tmp_NO)
            tmp_IO = np.array(tmp_IO)

            aveNO = np.mean(tmp_NO)
            aveIO = np.mean(tmp_IO)

            testIO.append(aveNO - aveIO)


    #ax.hist(testIO, bins=100, color="crimson", edgecolor="black", label="IO")

    #ax.set_xlabel("Test statistics", fontsize=13)

    #ax.legend(prop={"size":13})
    #
    #plt.tight_layout()
    #plt.savefig("KStest_2kpc.pdf")
    #plt.show()

    import uproot3
    mod = mod_arr[0]
    with uproot3.recreate("KStest_mod%d_10kpc.root"%mod) as f:
        f["KS"] = uproot3.newtree({"NO":"float64", "IO":"float64"})
        f["KS"].extend({"NO":testNO, "IO":testIO})








