import sys
import argparse
from calculate_tau import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='tau calculation')
    parser.add_argument('-f', type=str, required=True, default=False)
    parser.add_argument('-o', type=str, required=True, default=False)
    args = parser.parse_args()

print(args.f)

print("Open tpm file")
tpm=openExpData(args.f)
print("Done")
print("Filter zero and low expressed genes")
spExpFil=removeZeroExpGenes(tpm)
print(spExpFil.shape)
print("Done")
spExpFil2=removeLowExpressedGenes(spExpFil)
print(spExpFil2.shape)
print("Done")
tauDf=calculate_tau(spExpFil2)
print("Done")
print("Writing output file")
tauDf.to_csv(args.f+"tau.csv")
