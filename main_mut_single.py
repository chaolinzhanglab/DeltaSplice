from loguru import logger
import os
import argparse
import importlib
import random
import numpy as np
import torch
from utils import (
    MutGenerator,
    density_scatter
)
from scipy.stats import spearmanr, pearsonr
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

    # set seeds
    os.environ["PYTHONHASHSEED"] = str(1 + config.seed)
    os.environ["TF_DETERMINISTIC_OPS"] = "1"
    random.seed(1 + config.seed)
    np.random.seed(1 + config.seed)
    torch.random.manual_seed(1+config.seed)

    print("finish set env")
    # prepare model
    # train model

    assert config.model_path is not None, "model_path (the path to the trained models) must be set in config"
    assert config.mut_data is not None, "mut_data must be set in config"
    if not isinstance(config.mut_data, list):
        config.mut_data = [config.mut_data]
    if not isinstance(config.model_path, list):
        config.model_path = [config.model_path]

    save_dir = os.path.join(config.save_path, "test_results/")
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    Models = []
    for model_path in config.model_path:
        Models.append(copy.deepcopy(config.model))
        Models[-1].load_state_dict(torch.load(model_path))

    for path in config.mut_data:
        logger.info("Eval {} ".format(path))
        # load data
        test_data = DataLoader(MutGenerator(
            EL=config.EL, jsonfile=path
        ), batch_size=1, shuffle=False, num_workers=5)
        # evaluation mutations
        Pred_ref = []
        Pred_delta = []
        Gt_delta = []
        for index, d in enumerate(test_data):
            if index % 100 == 99:
                print("finish {}/{}".format(index, len(test_data)))
            d.pop("single_pred_psi")
            gt_delta = d["mutY"]  # GT of delta usage
            pred = [m.predict(d, use_ref=False) for m in Models]
            pred_ref = sum([v["single_pred_psi"] for v in pred])/len(pred)
            pred_delta = sum([v["mutY"] for v in pred])/len(pred)-pred_ref

            if len(gt_delta.shape) == 4:
                gt_delta = gt_delta[:, :, 0]
            assert len(pred_delta.shape) == 3
            assert len(gt_delta.shape) == 3

            gt_delta = gt_delta[:, :, 1:].reshape(-1)
            pred_ref = pred_ref[:, :, 1:].reshape(-1)
            pred_delta = pred_delta[:, :, 1:].reshape(-1)
            position = np.nonzero(gt_delta)
            if len(position) == 0:
                logger.warning(
                    "Encounting mutations without non-zero labels, will be skipped")
                continue
            gt_delta = gt_delta[position].mean()
            pred_ref = pred_ref[position].mean()
            pred_delta = pred_delta[position].mean()

            Pred_delta.append(pred_delta)
            Gt_delta.append(gt_delta)
            Pred_ref.append(pred_ref)
        """if not np.isnan(sum(Gt_delta)):
            spearmanR=spearmanr(Pred_delta, Gt_delta)
            pearsonR=pearsonr(Pred_delta, Gt_delta)
            logger.info("SpearmanR & PearsonR between predictions and GT are {}&{}".format(spearmanR, pearsonR))
            density_scatter(np.array(Gt_delta), np.array(Pred_delta), "GT", "Pred", os.path.join(
                save_dir, "Mutation_{}.png".format("{}+{}".format(path.replace("/", "_"), config.model_path[0].replace("/", "_")))))
        else:
            logger.warning("There are Nan in labels")"""
        with open(os.path.join(save_dir, "Single_Mutation_{}.txt".format("{}+{}".format(path.replace("/", "_"), config.model_path[0].replace("/", "_")))), "w") as f:
            for a, b, c in zip(Pred_delta, Gt_delta, Pred_ref):
                f.writelines("{}\t{}\t{}\n".format(c, a, b))


if __name__ == "__main__":
    main()
