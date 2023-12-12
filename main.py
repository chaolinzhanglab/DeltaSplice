from loguru import logger
import os
import argparse
import importlib
import random
import numpy as np
from utils import MutGenerator, GetSummaryStatisticsCallback, write_splice_site_file_header, write_splice_sites, get_top1_statistics, get_correlation, collect_predictions, density_scatter
import torch
from torch.utils.data import DataLoader
import copy


def main():
    # load config file
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", help="the path to the config file")
    args = parser.parse_args()
    config_path = args.config.replace("/", ".")
    config = importlib.import_module(config_path)
    config = config.config
    logger.add(config.log_path)

    # set seeds
    os.environ["PYTHONHASHSEED"] = str(1+config.seed)
    os.environ["TF_DETERMINISTIC_OPS"] = "1"
    random.seed(1+config.seed)
    np.random.seed(1+config.seed)
    torch.random.manual_seed(1+config.seed)
    print("finish set env")
    from utils import DataGenerator

    if config.is_train:
        # prepare model
        train_data = DataLoader(DataGenerator(EL=config.EL, jsonfile=config.trainjsonfile),
                                shuffle=True, batch_size=config.batch_size, drop_last=True)
        validation_data = DataLoader(DataGenerator(EL=config.EL, jsonfile=config.validjsonfile, ), shuffle=False,
                                     batch_size=config.batch_size, drop_last=True)  # only test on human data
        if hasattr(config, "mut_data") and config.mut_data is not None:
            mut_data = DataLoader(MutGenerator(
                jsonfile=config.mut_data, EL=config.EL,
            ), batch_size=1, shuffle=False)
        else:
            mut_data = None

        test_data = []  # in train mode, do not consider test data
        summary = GetSummaryStatisticsCallback(
            config.model,
            train_data, validation_data, test_data=test_data, mut_data=mut_data,
            summary_fout_path=os.path.join(config.save_path, "summary_log"),
            model_save_path=os.path.join(config.save_path, "models")
        )

        logger.info("Finish loading data ...")
        logger.info(
            "Train data size {} validation data size {} test data size {}".format(len(train_data), len(validation_data),
                                                                                  len(test_data)))

        summary.fit(config.epoch_num)
        logger.info("Finish training")
        return

    else:
        # load model
        assert config.model_path is not None
        assert config.testjsonfile is not None or config.mut_data is not None
        if config.testjsonfile is None and config.mut_data is not None:
            config.testjsonfile = config.mut_data
            DataGenerator = MutGenerator
        else:
            from utils import DataGenerator
        if not isinstance(config.testjsonfile, list):
            config.testjsonfile = [config.testjsonfile]
        if not isinstance(config.model_path, list):
            config.model_path = [config.model_path]
        save_dir = os.path.join(config.save_path, "test_results/")
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)

        Models = [copy.deepcopy(config.model) for _ in config.model_path]
        [m.load_state_dict(torch.load(b)) for m, b in zip(Models, config.model_path)]

        for path in config.testjsonfile:
            save_path_ = "{}+{}".format(path.replace("/", "_"), "eval"+config.model_path[0].replace("/", "_"))
            save_path_ = os.path.join(save_dir, save_path_)
            save_path = []
            test_data = DataLoader(DataGenerator(EL=config.EL, jsonfile=path), batch_size=config.batch_size, shuffle=False, num_workers=5)
            if config.save_files:
                Y_true_1 = []
                Y_true_2 = []
                Y_pred_1 = []
                Y_pred_2 = []
                save_path.append(open(save_path_+"_acceptor", "w"))
                save_path.append(open(save_path_+"_donor", "w"))
                [write_splice_site_file_header(_) for _ in save_path]
                for index, d in enumerate(test_data):
                    if index % 100 == 0:
                        print(index)
                    X, Y = d["X"], d["single_pred_psi"].numpy()
                    if "mutX" in d:
                        d.pop("mutX")

                    pred = sum([m.predict(d)["single_pred_psi"] for m in Models])/len(Models)
                    if len(Y.shape) == 4:
                        Y = Y[:, :, 0]
                    assert len(Y.shape) == len(pred.shape)
                    assert len(pred.shape) == 3

                    # note that tx_start/tx_end refer to the start/end of input sequence instead of gene start/end
                    CHROM = d["chrom"]
                    NAME = d["name"]
                    STRAND = d["strand"]
                    TX_START = d["txstart"].numpy()
                    TX_END = d["txend"].numpy()
                    SPECIES = d["species"]
                    write_splice_sites(save_path, CHROM, NAME, STRAND,
                                       TX_START, TX_END, SPECIES, pred, Y)
                    pred_sites = np.nonzero(Y[:, :, 1:].sum(-1))

                    Y_true_1.extend(Y[pred_sites][:, 1].flatten())
                    Y_true_2.extend(Y[pred_sites][:, 2].flatten())
                    Y_pred_1.extend(pred[pred_sites][:, 1].flatten())
                    Y_pred_2.extend(pred[pred_sites][:, 2].flatten())
                [_.close() for _ in save_path]

                topk1_acc_a, auprc_a, threshold_a, num_idx_true_a = get_top1_statistics(
                    np.asarray(Y_true_1), np.asarray(Y_pred_1)
                )
                logger.info("File {} Model {} Acceptor: topk1_acc {} auprc {} threshold {} num_idx_true {}".format(
                    path, config.model_path[0], topk1_acc_a, auprc_a, threshold_a, num_idx_true_a
                ))
                rho_a, spa,  num_idx_true_a2 = get_correlation(
                    np.asarray(Y_true_1), np.asarray(Y_pred_1)
                )

                logger.info("Acceptor correlation: rho {} pearson {} num_idx_true {}".format(
                    rho_a, spa,  num_idx_true_a2
                ))
                topk1_acc_d, auprc_d, threshold_d, num_idx_true_d = get_top1_statistics(
                    np.asarray(Y_true_2), np.asarray(Y_pred_2)
                )
                logger.info("File {} Model{} Donor: topk1_acc {} auprc {} threshold {} num_idx_true {}".format(
                    path, config.model_path[0], topk1_acc_d, auprc_d, threshold_d, num_idx_true_d
                ))
                rho_d, spd,  num_idx_true_d2 = get_correlation(
                    np.asarray(Y_true_2), np.asarray(Y_pred_2)
                )

                logger.info("Donor correlation: rho {} pearson {} num_idx_true {}".format(
                    rho_d, spd, num_idx_true_d2
                ))
                y_true1, y_true2 = np.copy(Y_true_1), np.copy(Y_true_2)
                y_true1[np.isnan(y_true1)] = -1
                y_true2[np.isnan(y_true2)] = -1

                idx_true1, idx_true2 = np.nonzero(y_true1 > 0+1e-10), np.nonzero(y_true2 > 0+1e-10)
                print(len(idx_true1[0]), len(idx_true2[0]))

                density_scatter(np.asarray(Y_true_1)[idx_true1], np.asarray(Y_pred_1)[idx_true1], "GT", "Pred", os.path.join(
                    save_dir, "Acceptor_{}.png".format("{}+{}".format(path.replace("/", "_"), config.model_path[0].replace("/", "_")))))

                density_scatter(np.asarray(Y_true_2)[idx_true2], np.asarray(Y_pred_2)[idx_true2], "GT", "Pred", os.path.join(
                    save_dir, "Donor_{}.png".format("{}+{}".format(path.replace("/", "_"), config.model_path[0].replace("/", "_")))))
                logger.info("finish {}".format(path))

            elif config.print_summary:
                Pred = [collect_predictions(m, test_data) for m in Models]
                Y_true_1, Y_true_2, Y_pred_1, Y_pred_2 = sum([p[0] for p in Pred])/len(Pred), sum([p[1] for p in Pred])/len(Pred), sum([p[2]
                                                                                                                                        for p in Pred])/len(Pred), sum([p[3] for p in Pred])/len(Pred)

                topk1_acc_a, auprc_a, threshold_a, num_idx_true_a = get_top1_statistics(
                    np.asarray(Y_true_1), np.asarray(Y_pred_1)
                )
                logger.info("File {} Model {} Acceptor: topk1_acc {} auprc {} threshold {} num_idx_true {}".format(
                    path, config.model_path[0], topk1_acc_a, auprc_a, threshold_a, num_idx_true_a
                ))
                rho_a, spa,  num_idx_true_a2 = get_correlation(
                    np.asarray(Y_true_1), np.asarray(Y_pred_1)
                )
                logger.info("Acceptor correlation: rho {} pearson {} num_idx_true {}".format(
                    rho_a, spa,  num_idx_true_a2
                ))
                topk1_acc_d, auprc_d, threshold_d, num_idx_true_d = get_top1_statistics(
                    np.asarray(Y_true_2), np.asarray(Y_pred_2)
                )
                logger.info("File {} Model{} Donor: topk1_acc {} auprc {} threshold {} num_idx_true {}".format(
                    path, config.model_path[0], topk1_acc_d, auprc_d, threshold_d, num_idx_true_d
                ))
                rho_d, spd,  num_idx_true_d2 = get_correlation(
                    np.asarray(Y_true_2), np.asarray(Y_pred_2)
                )
                logger.info("Donor correlation: rho {} pearson {} num_idx_true {}".format(
                    rho_d, spd, num_idx_true_d2
                ))


if __name__ == "__main__":
    main()
