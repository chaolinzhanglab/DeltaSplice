from argparse import ArgumentParser
from sys import stdin as std_in
from sys import stdout as std_out
import logging
import pysam
from deltasplice.pred_utils import get_delta_prediction, Annotator

def get_options():
    parser=ArgumentParser(description="Version: 1.0.0")
    parser.add_argument('-I', metavar='input', nargs='?', default=std_in,
                        help='path to the input VCF file, defaults to standard in')
    parser.add_argument('-O', metavar='output', nargs='?', default=std_out,
                        help='path to the output VCF file, defaults to standard out')
    parser.add_argument('-R', metavar='reference', required=True,
                        help='path to the reference genome fasta file')
    parser.add_argument('-A', metavar='annotation', required=True,
                        help='"grch37" (GENCODE V24lift37 canonical annotation file in '
                             'package), "grch38" (GENCODE V24 canonical annotation file in '
                             'package), or path to a similar custom gene annotation file')
    parser.add_argument('-D', metavar='distance', nargs='?', default=50,
                        type=int, choices=range(0, 5000),
                        help='maximum distance between the variant and gained/lost splice '
                             'site, defaults to 50')
    parser.add_argument('-U', metavar='use_reference_info', nargs='?', default=1,
                        type=int, choices=[0, 1],
                        help='whether to use reference information, default to 1')
    parser.add_argument('-M', metavar='mask', nargs='?', default=0,
                        type=int, choices=[0, 1],
                        help='mask scores representing annotated acceptor/donor gain and '
                             'unannotated acceptor/donor loss, defaults to 0')
    args = parser.parse_args()

    return args

def main():
    args=get_options()
    if None in [args.I, args.O, args.D]:
        logging.error('Usage: spliceai [-h] [-I [input]] [-O [output]] -R reference -A annotation'
                      '[-D [distance]][-U [use_reference_info]]') 
        exit()
    
    vcf=pysam.VariantFile(args.I)
    header = vcf.header
    header.add_line('##INFO=<ID=DeltaSplice,Number=.,Type=String,Description="DeltaSplicev1.0.0 variant '
                    'annotation. These include delta scores (DS) and delta positions (DP) for '
                    'acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). '
                    'Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">')
    output = pysam.VariantFile(args.O, mode='w', header=header)
    ann = Annotator(args.R, args.A)

    for record in vcf:
        scores = get_delta_prediction(record, args.D, args.R, args.U, ann, args.M)
        if len(scores) > 0:
            record.info['DeltaSplice'] = scores
        output.write(record)

    vcf.close()
    output.close()

if __name__=="__main__":
    main()