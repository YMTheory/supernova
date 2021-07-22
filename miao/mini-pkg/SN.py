from SNnumGarchingSrc import SNnumGarchingSrc

if __name__ == "__main__":
    GarSrc = SNnumGarchingSrc(82500, 10)
    print( GarSrc.totalSNFluenceDetTimeShift(1, 1.0, 10, 1) )
