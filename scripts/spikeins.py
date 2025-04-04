import argparse as arg
import pandas as pd

def percent(values):
    """
    get a list of values return a list of percentages
    :param values: list of values whatever it may be
    :return: list of percentages
    """
    total=sum(values)
    percentages=[]
    for value in values:
        if total==0:
            percentages.append(0)
        else:
            pct=value/total
            percentages.append(pct)
    return percentages


if __name__ == "__main__":
    parser = arg.ArgumentParser(description='combine detected ERCC/SIRV along with actual values for ERCC just a '
                                            'comparison of actual vs detected concentraion but in different units'
                                            'using TPM for all these measurements')
    parser.add_argument('--gene', type=str, help='gene expression values', action="store")
    parser.add_argument('--iso', type=str, help='isoform expression values', action="store")
    parser.add_argument('--actual', type=str, help='actual concentrations', action="store")
    parser.add_argument('--samplename', type=str, help='samplename', action="store")
    parser.add_argument('--analysis_path', type=str, help='analysis path', action="store")

    args = parser.parse_args()


    gene=pd.read_csv(args.gene, header=0, sep="\t")
    gene=gene[gene["gene_id"].str.contains("ERCC")]

    iso=pd.read_csv(args.iso, header=0, sep="\t")
    iso=iso[iso["gene_id"].str.contains("SIRV")]


    actual = pd.read_csv(args.actual, header=0, sep=",")

    gene_merged = pd.merge(gene[["gene_id", "TPM"]],
                           actual, how='inner', on='gene_id')

    gene_name=args.analysis_path+"/"+args.samplename+".ercc.txt"
    gene_merged.to_csv(gene_name, sep="\t", index=False, header=True)

    iso_name = args.analysis_path+"/"+args.samplename + ".sirv.txt"
    iso_merged = pd.merge(iso[["transcript_id", "TPM", "IsoPct"]],
                           actual, how='inner', on='transcript_id')

    iso_merged.to_csv(iso_name, sep="\t", index=False, header=True)
