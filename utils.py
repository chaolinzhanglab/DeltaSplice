import numpy as np
from scipy.stats import spearmanr, pearsonr
import os
from loguru import logger
from sklearn.metrics import average_precision_score
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.interpolate import interpn
from torch.utils.data import Dataset
import torch
import json
from config import Fapath, IN_MAP, repdict, SeqTable
from pyfasta import Fasta


def parse_bed_line(line):
    eps = 1e-3
    line = line.replace("\n", "").split("\t")
    fa, hgchrom, name, _, chrom, strand, start, end, ss3s, ss5s, ss3vs, ss5vs = line
    ss3s = [int(x)-1 for x in ss3s.split(",") if len(x) > 0]  # with 0 base
    ss5s = [int(x)-1 for x in ss5s.split(",") if len(x) > 0]
    ss3vs = [max(float(x), eps) for x in ss3vs.split(",") if len(x) > 0]
    ss5vs = [max(float(x), eps) for x in ss5vs.split(",") if len(x) > 0]
    return fa, chrom, name, strand, int(start)-1, int(end)-1, ss3s, ss5s, ss3vs, ss5vs, hgchrom


def get_anno_list():
    with open("AnnoList.txt", "r") as f:
        content = f.readlines()
    content = [x.replace("\n", "") for x in content]
    return content


def get_species_list():
    with open("SpeciesList.txt", "r") as f:
        content = f.readlines()
    content = [x.replace("\n", "") for x in content]
    content = [x[0].upper()+x[1:] for x in content]
    return content


def density_scatter(x, y, x_name, y_name, savepath, bins=50, sort=True):
    plt.figure()

    fig, ax = plt.subplots()
    data, x_e, y_e = np.histogram2d(x, y, bins=bins, density=True)
    z = interpn((0.5*(x_e[1:] + x_e[:-1]), 0.5*(y_e[1:]+y_e[:-1])), data, np.vstack([x, y]).T, method="splinef2d", bounds_error=False)
    z[np.where(np.isnan(z))] = 0.
    if sort:
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    ax.scatter(x, y, c=z)
    norm = Normalize(vmin=np.min(z), vmax=np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm), ax=ax)
    cbar.ax.set_ylabel('Density')

    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.savefig(savepath)


class MutGenerator(Dataset):
    def __init__(self,  EL, jsonfile):
        with open(jsonfile, "r") as f:
            self.data = json.load(f)
        self.fafiles = {}
        self.EL = EL

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        d = self.data[index]
        if "mutseq" in d:
            seq, mutseq = d["seq"], d["mutseq"]
            strand = d["strand"]
            start, end = int(d["start"]), int(d["end"])
        else:
            species, chrom, start, end, strand, mutpos, oriv, mutv = d["species"], d["chrom"], int(d["start"]), int(d["end"]), d["strand"], int(d["mutpos"]), d["oriv"].upper(), d["mutv"].upper()
            if species not in self.fafiles:
                self.fafiles[species] = Fasta(os.path.join(Fapath, species+".fa"))
            seq = self.fafiles[species][chrom][start:end].upper()
            assert seq[mutpos] == oriv
            mutseq = seq[:mutpos]+mutv+seq[mutpos+1:]
        assert len(seq)>=self.EL+5000, "the length of sequence {} should be large than {}".format(len(seq), self.EL+5000)
        assert len(mutseq)==len(seq)
        label = np.zeros([end-start-self.EL, 3])
        mutlabel = np.zeros([end-start-self.EL, 3])

        for v in d["label"]:
            idx, value = v
            idx = int(idx)
            if idx-self.EL//2 >= 0 and idx< end-start-self.EL//2:
                value = np.array([float(_) for _ in value])
                label[idx-self.EL//2] = value

        for v in d["mutlabel"]:
            idx, value = v
            idx = int(idx)
            assert idx >= self.EL//2 and idx < end-start-self.EL//2
            value = np.array([float(_) for _ in value])
            mutlabel[idx-self.EL//2] = value

        if strand == "-":
            seq = [repdict[_] for _ in seq][::-1]
            mutseq = [repdict[_] for _ in mutseq][::-1]
            label = np.copy(label[::-1])
            mutlabel = np.copy(mutlabel[::-1])

        seq = IN_MAP[[SeqTable[_] for _ in seq]][:, :4]
        mutseq = IN_MAP[[SeqTable[_] for _ in mutseq]][:, :4]

        label[:, 0] = 1.-label[:, 1:].sum(-1)
        return {"X": seq, "mutX": mutseq, "single_pred_psi": label, "mutY": mutlabel, "chrom": d["chrom"], "strand": d["strand"], "species": "hg19", "txstart": d["start"], "txend": d["end"], "name": "mut"}


class DataGenerator(Dataset):
    def __init__(self, EL, jsonfile=None):
        with open(jsonfile, "r") as f:
            self.data = json.load(f)
        self.fafiles = {}
        self.EL = EL

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        d = self.data[index]
        species, chrom, start, end, strand, name = d["species"], d["chrom"], int(d["start"]), int(d["end"]), d["strand"], d["name"]
        if species not in self.fafiles:
            self.fafiles[species] = Fasta(os.path.join(Fapath, species+".fa"))

        seq = self.fafiles[species][chrom][start:end].upper()
        label = np.zeros([end-start-self.EL, 3])
        for v in d["label"]:
            idx, value = v
            idx = int(idx)
            if idx-self.EL//2 >= 0 and idx < end-start-self.EL//2:
                value = np.array([float(_) for _ in value])
                label[idx-self.EL//2] = value
        if strand == "-":
            seq = [repdict[_] for _ in seq][::-1]
            label = np.copy(label[::-1])
        seq = IN_MAP[[SeqTable[_] for _ in seq]][:, :4]
        label[:, 0] = 1.-label[:, 1:].sum(-1)

        return {"X": seq, "single_pred_psi": label, "species": species, "chrom": chrom, "name": name, "strand": strand, "txstart": start, "txend": end}


class GetSummaryStatisticsCallback():
    def __init__(
        self,
        model,
        train_data,
        validation_data,
        test_data=None,
        mut_data=None,
        summary_fout_path=None,
        model_save_path=None,
    ):
        super().__init__()
        if isinstance(test_data, list):
            if len(test_data) == 0:
                test_data = None
        self.model = model
        self.train_data = train_data
        self.validation_data = validation_data
        self.test_data = test_data
        self.mut_data = mut_data
        self.model_save_path = model_save_path

        # create files to save eval results in the training process
        if not os.path.exists(summary_fout_path):
            os.mkdir(summary_fout_path)
        self.train_fout = open(
            os.path.join(summary_fout_path, "train_data.eval.log"), "w"
        )
        self.validation_fout = open(
            os.path.join(summary_fout_path, "validation_data.eval.log"), "w"
        )
        if not self.test_data is None:
            self.test_fout = open(
                os.path.join(summary_fout_path, "test_data.eval.log"), "w"
            )
        else:
            self.test_fout = None

        # def on_train_begin(self, epoch, logs=None):
        # write headers to files
        print_summary_file_header(self.train_fout)
        print_summary_file_header(self.validation_fout)
        if not self.test_fout is None:
            print_summary_file_header(self.test_fout)
        if self.model_save_path is not None and not os.path.exists(self.model_save_path):
            os.mkdir(self.model_save_path)

    def train_one_epoch(self):
        tloss = 0

        for i, d in enumerate(self.train_data):
            loss = self.model.train_onestep(d)
            logout = "|".join([k+":"+str('%.3g' % loss[k]) for k in loss])
            logger.info(logout)

            if i % 500 == 0 and self.model_save_path is not None:
                torch.save(self.model.state_dict(), os.path.join(self.model_save_path, "model.ckpt-temp"))
                if self.mut_data is not None:
                    logger.info("Mut data")
                    Pred, Mutv = [], []

                    for index, d in enumerate(self.mut_data):
                        if index % 100 == 99:
                            print(index)
                            print(spearmanr(Pred, Mutv), pearsonr(Pred, Mutv), min(Pred), min(Mutv))
                        Y = d["mutY"]
                        TPRED = self.model.predict(d)
                        gpred = TPRED["single_pred_psi"]
                        pred = TPRED["mutY"]-gpred
                        if len(Y.shape) == 4:
                            Y = Y[:, :, 0]
                        assert len(pred.shape) == 3
                        assert len(Y.shape) == 3
                        mutv = Y[:, :, 1:].reshape(-1)
                        gpred = gpred[:, :, 1:].reshape(-1)
                        pred = pred[:, :, 1:].reshape(-1)
                        position = np.nonzero(mutv)
                        if len(position[0]) == 0:
                            print("skip {}".format(index))
                            continue
                        mutv = mutv[position].mean()
                        pred = pred[position].mean()
                        gpred = gpred[position].mean()

                        if abs(mutv) > 1e-10:
                            Pred.append(pred)
                            Mutv.append(mutv)
                    logger.info("SpearManR is {} PearsonR is {}".format(spearmanr(Pred, Mutv), pearsonr(Pred, Mutv)))
        if self.model.scheduler is not None:
            for v in self.model.scheduler:
                v.step()
                logger.info("scheduler learning rate to {}".format(self.model.optimizer[0].param_groups[0]['lr']))

        return {"loss": tloss, "val_loss": float("inf")}

    def fit(self, num_epoch):
        for epoch in range(num_epoch):
            logs = self.train_one_epoch()
            self.on_epoch_end(epoch, logs)

    def _write_line(self, data, fout, epoch, logs):
        line = "\t".join(
            map(
                str,
                [
                    epoch,
                    self.model.optimizer[0].param_groups[0]['lr'],
                    logs["loss"],
                    logs["val_loss"],
                ],
            )
        )
        Y_true_1, Y_true_2, Y_pred_1, Y_pred_2 = collect_predictions(
            self.model, data
        )

        topk1_acc_a, auprc_a, threshold_a, num_idx_true_a = get_top1_statistics(
            np.asarray(Y_true_1), np.asarray(Y_pred_1)
        )
        logger.info(
            "Acceptor: topk1_acc {} auprc {} threshold {} num_idx_true {}".format(
                topk1_acc_a, auprc_a, threshold_a, num_idx_true_a
            )
        )
        topk1_acc_d, auprc_d, threshold_d, num_idx_true_d = get_top1_statistics(
            np.asarray(Y_true_2), np.asarray(Y_pred_2)
        )
        logger.info(
            "Donor: topk1_acc {} auprc {} threshold {} num_idx_true {}".format(
                topk1_acc_d, auprc_d, threshold_d, num_idx_true_d
            )
        )
        line += "\t" + "\t".join(
            map(
                str,
                [
                    topk1_acc_a,
                    auprc_a,
                    threshold_a,
                    num_idx_true_a,
                    topk1_acc_d,
                    auprc_d,
                    threshold_d,
                    num_idx_true_d,
                ],
            )
        )

        rho_a, spa, num_idx_true_a2 = get_correlation(
            np.asarray(Y_true_1), np.asarray(Y_pred_1)
        )
        logger.info(
            "Acceptor correlation: rho {} pearson {} num_idx_true {}".format(
                rho_a, spa, num_idx_true_a2
            )
        )
        rho_d, spd, num_idx_true_d2 = get_correlation(
            np.asarray(Y_true_2), np.asarray(Y_pred_2)
        )
        logger.info(
            "Donor correlation: rho {} pearson {} num_idx_true {}".format(
                rho_d, spd, num_idx_true_d2
            )
        )
        line += "\t" + "\t".join(map(str, rho_a + rho_d))

        fout.write(line + os.linesep)
        fout.flush()

    def on_epoch_end(self, epoch, logs={}):
        # save model

        if self.model_save_path is not None:
            torch.save(self.model.state_dict(), os.path.join(self.model_save_path, "model.ckpt-{}".format(epoch)))
        # write logs
        logger.info("Epoch: {}".format(epoch))
        logger.info("Validation data:")
        self._write_line(self.validation_data, self.validation_fout, epoch, logs)

        if not self.test_data is None:
            logger.info("Test data")
            self._write_line(self.test_data, self.test_fout, epoch, logs)


def collect_predictions(model, data):
    Y_true_1 = []
    Y_true_2 = []
    Y_pred_1 = []
    Y_pred_2 = []

    for d in data:
        X, Y = d["X"], d["single_pred_psi"].numpy()
        if len(Y.shape) == 4:
            Y = Y[:, :, 0]
        pred = model.predict(d)
        pred = pred["single_pred_psi"]

        assert len(pred.shape) == len(Y.shape)
        assert len(pred.shape) == 3

        Y_true_1.extend(Y[:, :, 1].flatten())
        Y_true_2.extend(Y[:, :, 2].flatten())
        Y_pred_1.extend(pred[:, :, 1].flatten())
        Y_pred_2.extend(pred[:, :, 2].flatten())
    return Y_true_1, Y_true_2, Y_pred_1, Y_pred_2


def get_top1_statistics(y_true, y_pred, eps=1e-10):
    y_true = np.copy(y_true)
    y_true[np.isnan(y_true)] = 1.0

    idx_true = np.nonzero(y_true > eps)[0]
    argsorted_y_pred = np.argsort(y_pred)
    sorted_y_pred = y_pred[argsorted_y_pred]
    idx_pred = argsorted_y_pred[-len(idx_true):]
    topk1_acc = (
        np.size(np.intersect1d(idx_true, idx_pred))
        * 1.0
        / (min(len(idx_pred), len(idx_true)) + eps)
    )
    threshold = sorted_y_pred[-(len(idx_true)-1)]
    auprc = average_precision_score((y_true > eps).astype(int), y_pred)
    return topk1_acc, auprc, threshold, len(idx_true)


def print_top1_statistics(topk1_accuracy, auprc, threshold, num_idx_true):
    print(
        ("\033[91m%.4f\t\033[0m\033[94m%.4f\t\033[0m%.4f\t%d")
        % (topk1_accuracy, auprc, threshold, num_idx_true)
    )


def get_correlation(y_true, y_pred, eps=1e-10):
    y_true = np.copy(y_true)
    y_true[np.isnan(y_true)] = -1

    rho = []
    pearson_r = []
    num_idx_true = []

    for psi_t in [0, 0.1, 0.2, 0.3]:
        idx_true = np.nonzero(
            np.logical_and(y_true >= psi_t + eps, y_true < 1.0 - psi_t)
        )[0]
        rho1, pval1 = spearmanr(y_true[idx_true], y_pred[idx_true])
        rho.append(rho1)
        pearson_r.append(pearsonr(y_true[idx_true], y_pred[idx_true])[0])
        num_idx_true.append(np.size(idx_true))
    return rho, pearson_r, num_idx_true


def print_correlation(rho, num_idx_true):
    print(
        ("\033[91m%.4f(%d)\033[0m\t%.4f(%d)\t%.4f(%d)\t%.4f(%d)")
        % (
            rho[0],
            num_idx_true[0],
            rho[1],
            num_idx_true[1],
            rho[2],
            num_idx_true[2],
            rho[3],
            num_idx_true[3],
        )
    )


def print_summary_file_header(fout):
    line = "\t".join(["epoch", "learning_rate", "train_loss", "val_loss"])

    topk_cols = [
        "acceptor_topk_accuracy",
        "acceptor_auprc",
        "acceptor_threshold",
        "acceptor_site_num",
        "donor_topk_accuracy",
        "donor_auprc",
        "donor_threshold",
        "donor_site_num",
    ]
    cor_cols = [
        "acceptor_cor_all",
        "validation_acceptor_cor_0.1",
        "validation_acceptor_cor_0.2",
        "acceptor_cor_0.3",
        "donor_cor_all",
        "donor_cor_0.1",
        "validation_donor_cor_0.2",
        "donor_cor_0.3",
    ]

    cols = topk_cols + cor_cols

    line += "\t" + "\t".join(cols)

    fout.write(line + os.linesep)
    fout.flush()


def write_splice_site_file_header(fout):
    yt = ["Yt"]
    yp = ["Yp"]
    header = "\t".join(
        ["#species", "chrom", "chromStart", "chromEnd", "name", "score", "strand", "type"]
        + yt
        + yp
    )
    fout.write(header + os.linesep)
    fout.flush()


def write_splice_sites(fouts, CHROM, NAME, STRAND, TX_START, TX_END, SPECIES, Yp, Yt):
    output_class_labels = ["Null", "acceptor", "donor"]
    # The three neurons per output correspond to no splicing, splice acceptor (AG)
    # and splice donor (GT) respectively.

    num_row = 0

    if num_row == 0:
        num_row = Yt.shape[0]
    else:
        assert num_row == Yt.shape[0]
        assert num_row == Yp.shape[0]

    for i in range(num_row):
        chrom = CHROM[
            i
        ]  # .decode()  #TODO: there is compatibility issue with more recent version of h5py
        name = NAME[i]  # .decode()
        strand = STRAND[i]  # .decode()
        tx_start = TX_START[i]
        tx_end = TX_END[i]
        species = SPECIES[i]

        for c in [1, 2]:
            fo = fouts[c - 1]

            y_t = np.copy(Yt[i, :, c])
            y_t[np.isnan(y_t)] = 1
            idx = np.nonzero(y_t > 1e-10)[0]  # positions of true splice sites

            for j in idx:
                pos = (tx_start + tx_end)//2

                yt = [str(Yt[i, j, c])]
                yp = [str(Yp[i, j, c])]

                line = "\t".join(
                    map(
                        str,
                        [
                            species,
                            chrom,
                            pos,
                            pos+1,
                            name,
                            0,
                            strand,
                            output_class_labels[c],
                        ]
                        + yt
                        + yp,
                    )
                )
                fo.write(line + os.linesep)
                fo.flush()
