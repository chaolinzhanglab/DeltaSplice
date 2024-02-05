from deltasplice.utils import  write_splice_sites, write_splice_site_file_header
import numpy as np
from loguru import logger
def get_prediction(inp, models, use_ref, ys=None):
    models[0].predict(inp)
    pred=sum([m.predict(inp, use_ref=use_ref)["single_pred_psi"]
                            for m in models])/len(models)
    if ys is None:
        return pred
    if len(ys.shape)==4:
        ys=ys[:, : , 0]
    assert len(ys.shape)==len(pred.shape)
    ly=ys.shape[1]
    bias=pred.shape[1]-ly
    assert bias%2==0
    pred=pred[:, bias//2:bias//2+ly]
    assert ys.shape==pred.shape
    
    return pred

def eval_test_data(data, models, save_path, use_ref):
    save_paths=(open(save_path+"_acceptor", "w"),open(save_path+"_donor", "w"))
    [write_splice_site_file_header(_) for _ in save_paths]
    Y_true_1 = []
    Y_true_2 = []
    Y_pred_1 = []
    Y_pred_2 = []
    for i, d in enumerate(data):
        if i%100==0:
            print("finish eval ", i)
        X, Y = d["X"], d["single_pred_psi"].numpy()
        pred=get_prediction(d, models, False, Y)
        CHROM = d["chrom"]
        NAME = d["name"]
        STRAND = d["strand"]
        TX_START = d["txstart"].numpy()
        TX_END = d["txend"].numpy()
        SPECIES = d["species"]
        write_splice_sites(save_paths, CHROM, NAME, STRAND,
                                    TX_START, TX_END, SPECIES, pred, Y)
        
        pred_sites = np.nonzero(Y[:, :, 1:].sum(-1))

        Y_true_1.extend(Y[pred_sites][:, 1].flatten())
        Y_true_2.extend(Y[pred_sites][:, 2].flatten())
        Y_pred_1.extend(pred[pred_sites][:, 1].flatten())
        Y_pred_2.extend(pred[pred_sites][:, 2].flatten())
        
    [_.close() for _ in save_paths]
    
def eval_mut_data(data, models, save_path, use_ref):
    Pred_ref = []
    Pred_delta = []
    Gt_delta = []
    for i, d in enumerate(data):
        if i%100==0:
            print("finish eval ", i)
        gt_delta = d["mutY"]  # GT of delta usage
        pred = [m.predict(d, use_ref=use_ref) for m in models]
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
   
    with open(save_path, "w") as f:
        f.writelines("Pred_delta\tGt_delta\tPred_ref\n")
        for a, b, c in zip(Pred_delta, Gt_delta, Pred_ref):
            f.writelines("{}\t{}\t{}\n".format(c, a, b))