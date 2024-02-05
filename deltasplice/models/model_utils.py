import random
import numpy as np
import torch
import sys


def neg_loss(eps=1e-10):
    def loss_func(y_true: torch.tensor, y_pred: torch.tensor):
        assert len(y_pred.shape) == len(y_true.shape)
        #y_true = y_true*tf.math.sign(y_pred)
        loss = (y_true*y_pred).mean()

        return -loss
    return loss_func


def l2_norm_dist():
    def loss_func(y_true: torch.tensor, y_pred: torch.tensor):
        #print(y_true.shape, y_pred.shape)
        assert len(y_pred.shape) == len(y_true.shape)
        y_true[torch.isnan(y_true)] = 0
        mask = (y_true.abs() > 0).to(float)

        loss = (mask*(y_true - y_pred*mask) ** 2).mean()
        return loss
    return loss_func


def get_available_gpus():
    return torch.cuda.device_count()


def cross_entropy(eps=1e-10):
    def loss_func(y_true: torch.tensor, y_pred: torch.tensor):
        assert len(y_pred.shape) == len(y_true.shape)
        if not y_true.requires_grad:
            y_true[torch.isnan(y_true)] = 0
            y_true = y_true.clamp(0.)
        loss = (-y_true*torch.log(y_pred.clamp(eps, 1.))).sum(-1).mean()
        return loss
    return loss_func


def cross_entropy_with_clamp(eps=1e-10, lower_bound=-1, upperbound=float("inf")):
    def loss_func(y_true: torch.tensor, y_pred: torch.tensor):
        assert len(y_pred.shape) == len(y_true.shape)
        y_true[torch.isnan(y_true)] = 0

        loss = (-y_true*torch.log(y_pred+eps)).sum(-1).clamp(lower_bound, upperbound).mean()
        return loss
    return loss_func
