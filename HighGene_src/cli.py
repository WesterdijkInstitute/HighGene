import argparse
from .HighGene import *
from argparse import RawTextHelpFormatter

def main():
    argParser = argparse.ArgumentParser(description="HighGene - Highlighting gene content variations, \
For more information see: https://github.com/WesterdijkInstitute/HighGene", formatter_class=RawTextHelpFormatter)
    argParser.add_argument("-i", "--input", type=str, help="Input folder", required=True)
    argParser.add_argument("-o", "--output", type=str, help="Output folder + File name", required=True)
    argParser.add_argument("-m", "--mode", type=int, choices=range(0, 6), help="Analysis modes:\n\
    0: Custom File\n\
    1: Transcription Factors\n\
    2: Peptidase Clans\n\
    3: CAZymes\n\
    4: KEGG\n\
    5: KOG\n\n", required=True)

    args = argParser.parse_args()
    if args.mode == 0:
        settings = input("Allow multiple factors to occur on the same gene? (y/n): ")
        HeatGenerator(CustomData(args.input, settings).custom_reader(), args.output).generate()
    if args.mode == 1:
        HeatGenerator(TfData(args.input).tf_reader(), args.output).generate()
    if args.mode == 2:
        settings = input("Visualise Clans (0) or Families (1)? (0/1) :")
        HeatGenerator(PepData(args.input, settings).pep_reader(), args.output).generate()
    if args.mode == 3:
        HeatGenerator(CazyData(args.input).cazy_reader(), args.output).generate()
    if args.mode == 4:
        HeatGenerator(KeggData(args.input).kegg_reader(), args.output).generate()
    if args.mode == 5:
        HeatGenerator(KogData(args.input).kog_reader(), args.output).generate()

    print(f"Figure saved as {args.output}.html. Thank you for using GeCCo!")

if __name__ == "__main__":
    main()
