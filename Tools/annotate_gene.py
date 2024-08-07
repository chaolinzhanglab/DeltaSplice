from loguru import logger
import os
import numpy as np
from deltasplice.utils import parse_bed_line
import json


def annotate(Lines, WDict):
    for line in Lines:
        Fa, Chr = line[:2]
        if Fa not in WDict:
            WDict[Fa] = {}
        if Chr not in WDict[Fa]:
            WDict[Fa][Chr] = {}
            WDict[Fa][Chr]["+"] = {}
            WDict[Fa][Chr]["-"] = {}

        fa, chrom, _, strand, _, _, ss3s, ss5s, ss3vs, ss5vs, _ = line
        assert len(ss3s) == len(ss3vs)
        assert len(ss5s) == len(ss5vs)
        assert chrom == Chr
        if strand == "+":
            for idx, v in zip(ss3s, ss3vs):
                assert np.isnan(v) or (v > 0 and v <= 1.0)
                assert idx not in WDict[fa][chrom]["+"] or WDict[fa][chrom]["+"][idx][2] == 0 or WDict[fa][chrom]["+"][idx][2] == v or np.isnan(v) or np.isnan(WDict[fa][chrom]["+"][idx][2])
                if idx not in WDict[fa][chrom]["+"] or WDict[fa][chrom]["+"][idx][2] == 0 or np.isnan(WDict[fa][chrom]["+"][idx][2]):
                    if idx not in WDict[fa][chrom]["+"]:
                        WDict[fa][chrom]["+"][idx] = [0, 0, 0]
                    if v==0:
                        v=1e-3
                    WDict[fa][chrom]["+"][idx][2] = v
            for idx, v in zip(ss5s, ss5vs):
                assert np.isnan(v) or (v > 0 and v <= 1.0)
                assert idx not in WDict[fa][chrom]["+"] or WDict[fa][chrom]["+"][idx][1] == 0 or WDict[fa][chrom]["+"][idx][1] == v or np.isnan(v) or np.isnan(WDict[fa][chrom]["+"][idx][1])
                if idx not in WDict[fa][chrom]["+"] or WDict[fa][chrom]["+"][idx][1] == 0 or np.isnan(WDict[fa][chrom]["+"][idx][1]):
                    if idx not in WDict[fa][chrom]["+"]:
                        WDict[fa][chrom]["+"][idx] = [0, 0, 0]
                    if v==0:
                        v=1e-3
                    WDict[fa][chrom]["+"][idx][1] = v
        else:
            assert strand == "-"
            for idx, v in zip(ss3s, ss3vs):
                assert np.isnan(v) or (v > 0 and v <= 1.0)
                assert idx not in WDict[fa][chrom]["-"] or WDict[fa][chrom]["-"][idx][1] == 0 or WDict[fa][chrom]["-"][idx][1] == v or np.isnan(v) or np.isnan(WDict[fa][chrom]["-"][idx][1])
                if idx not in WDict[fa][chrom]["-"] or WDict[fa][chrom]["-"][idx][1] == 0 or np.isnan(WDict[fa][chrom]["-"][idx][1]):
                    if idx not in WDict[fa][chrom]["-"]:
                        WDict[fa][chrom]["-"][idx] = [0, 0, 0]
                    if v==0:
                        v=1e-3
                    WDict[fa][chrom]["-"][idx][1] = v
            for idx, v in zip(ss5s, ss5vs):
                assert np.isnan(v) or (v > 0 and v <= 1.0)
                assert idx not in WDict[fa][chrom]["-"] or WDict[fa][chrom]["-"][idx][2] == 0 or WDict[fa][chrom]["-"][idx][2] == v or np.isnan(v) or np.isnan(WDict[fa][chrom]["-"][idx][2])
                if idx not in WDict[fa][chrom]["-"] or WDict[fa][chrom]["-"][idx][2] == 0 or np.isnan(WDict[fa][chrom]["-"][idx][2]):
                    if idx not in WDict[fa][chrom]["-"]:
                        WDict[fa][chrom]["-"][idx] = [0, 0, 0]
                    if v==0:
                        v=1e-3
                    WDict[fa][chrom]["-"][idx][2] = v

    logger.info("Finish writeing {}/{}".format(Fa, Chr))


def main():
    with open(TableFile, "r") as f:
        Table = f.readlines()
    Table = [parse_bed_line(_) for _ in Table]
    WDict = {}
    annotate(Table, WDict)
    with open(os.path.join(AnnoPath, "data.json"), "w") as f:
        f.write(json.dumps(WDict))


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser=ArgumentParser()
    parser.add_argument("--save_path", required=True, help="path to save output files", type=str)
    parser.add_argument("--input_file", required=True, help="path to input file", type=str)
    args=parser.parse_args()
    AnnoPath=args.save_path
    TableFile=args.input_file
    if not os.path.exists(AnnoPath):
        os.mkdir(AnnoPath)
    main()
