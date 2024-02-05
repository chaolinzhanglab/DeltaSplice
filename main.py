from loguru import logger
import os
import argparse
import random
import numpy as np
from deltasplice.utils import MutGenerator, GetSummaryStatisticsCallback
from deltasplice.models.model_utils import (
    get_available_gpus,
)
import torch
from torch.utils.data import DataLoader
import copy
import deltasplice.constant as d_constant
from deltasplice.pred_utils import eval_test_data, eval_mut_data
def main():
    # load config file
    parser = argparse.ArgumentParser()
    parser.add_argument("--save_path", required=True, help="path to save files")
    parser.add_argument("--is_train", default=False, help="whether to train", type=bool)
    parser.add_argument("--use_reference", default=False, help="whether to use reference in the prediction process", type=bool)
    parser.add_argument("--train_data_path", default=None, help="the path to train data", type=str)
    parser.add_argument("--test_data_path", default=None, help="the path to test data", type=str, nargs="*")
    parser.add_argument("--mut_data_path", default=None, help="the path to mutation data", type=str, nargs="*")
    parser.add_argument("--valid_data_path", default=None, help="the path to validation data", type=str)
    parser.add_argument("--train_num_epochs", default=5, help="the path to validation data", type=int)
    parser.add_argument("--num_workers", default=5, help="the number of workers for data loading")
    parser.add_argument("--batch_size_per_gpu", default=8, help="batch size for each gpu")
    parser.add_argument("--seed", default=321, help="random seed", type=int)
    parser.add_argument("--load_model_path", default=None, help="the path to load data in the eval mode", type=str, nargs="*")
    args = parser.parse_args()
    
    logger.add(os.path.join(args.save_path, "log"))

    # set seeds
    seed=args.seed+1
    os.environ["PYTHONHASHSEED"] = str(seed)
    os.environ["TF_DETERMINISTIC_OPS"] = "1"
    random.seed(seed)
    np.random.seed(seed)
    torch.random.manual_seed(seed)
    batch_size=args.batch_size_per_gpu*get_available_gpus()
    print("finish set env")
    

    if args.is_train:
        from deltasplice.utils import DataGenerator
        assert args.train_data_path is not None and args.valid_data_path is not None, "when using training mode, train_data_path and valid_data_path must be set"
        assert args.test_data_path is None and args.mut_data_path is None and args.use_reference is False, "in the training mode, test_data_path, mut_data_path and use_reference should not be set"
        
        # prepare model
        train_data = DataLoader(DataGenerator(EL=d_constant.EL, jsonfile=args.train_data_path),
                                shuffle=True, batch_size=batch_size, drop_last=True, num_workers=args.num_workers)
        validation_data = DataLoader(DataGenerator(EL=d_constant.EL, jsonfile=args.valid_data_path), shuffle=False,
                                     batch_size=batch_size, drop_last=True, num_workers=args.num_workers)  

        summary = GetSummaryStatisticsCallback(
            d_constant.model,
            train_data, validation_data, test_data=[], mut_data=None,
            summary_fout_path=os.path.join(args.save_path, "summary_log"),
            model_save_path=os.path.join(args.save_path, "models")
        )

        logger.info("Finish loading data ...")
        logger.info(
            "Train data size {} validation data size {}".format(len(train_data), len(validation_data),
                                                                                  ))
        summary.fit(args.train_num_epochs)
        logger.info("Finish training")
        return

    else:
        # load model
        assert args.load_model_path is not None, "the paths to load data should be set in the evaluation mode"
        assert args.test_data_path is not None or args.mut_data_path is not None, "one of test_data_path and mut_data_path should be set"
        assert not (args.test_data_path is not None and args.mut_data_path is not None), "only one of test_data_path and mut_data_path should be set"
        
        if args.mut_data_path is not None: 
            DataGenerator = MutGenerator
            data_path=args.mut_data_path
            func=eval_mut_data
        else:
            from deltasplice.utils import DataGenerator
            data_path=args.test_data_path
            func=eval_test_data
       
        save_dir = os.path.join(args.save_path, "test_results/")
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)

        Models = [copy.deepcopy(d_constant.model) for _ in args.load_model_path]
        [m.load_state_dict(torch.load(b))
         for m, b in zip(Models, args.load_model_path)]

        for path in data_path:
            save_path= "{}+{}".format(path.replace("/", "_"),
                                        args.load_model_path[0].replace("/", "_"))
            save_path = os.path.join(save_dir, save_path)
            
            test_data = DataLoader(DataGenerator(
                EL=d_constant.EL, jsonfile=path), batch_size=batch_size, shuffle=False, num_workers=args.num_workers)
            func(test_data, Models, save_path, args.use_reference)
            
          

if __name__ == "__main__":
    main()
