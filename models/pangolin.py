import numpy as np
import torch
import torch.nn.functional as F
import torch.nn as nn
from .model_utils import cross_entropy
from loguru import logger
# code from https://github.com/tkzeng/Pangolin
L = 32
# convolution window size in residual units
W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11,
                21, 21, 21, 21, 41, 41, 41, 41])
# atrous rate in residual units
AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4,
                 10, 10, 10, 10, 25, 25, 25, 25])


class ResBlock(nn.Module):
    def __init__(self, L, W, AR, pad=True):
        super(ResBlock, self).__init__()
        self.bn1 = nn.BatchNorm1d(L)
        s = 1
        # padding calculation: https://discuss.pytorch.org/t/how-to-keep-the-shape-of-input-and-output-same-when-dilation-conv/14338/2
        if pad:
            padding = int(1 / 2 * (1 - L + AR * (W - 1) - s + L * s))
        else:
            padding = 0
        self.conv1 = nn.Conv1d(L, L, W, dilation=AR, padding=padding)
        self.bn2 = nn.BatchNorm1d(L)
        self.conv2 = nn.Conv1d(L, L, W, dilation=AR, padding=padding)

    def forward(self, x):
        out = self.bn1(x)
        out = torch.relu(out)
        out = self.conv1(out)
        out = self.bn2(out)
        out = torch.relu(out)
        out = self.conv2(out)
        out = out + x
        return out


class Pangolin(nn.Module):
    def __init__(self, L, W, AR):
        super(Pangolin, self).__init__()
        self.n_chans = L
        self.conv1 = nn.Conv1d(4, L, 1)
        self.skip = nn.Conv1d(L, L, 1)
        self.resblocks, self.convs = nn.ModuleList(), nn.ModuleList()
        for i in range(len(W)):
            self.resblocks.append(ResBlock(L, W[i], AR[i]))
            if (((i + 1) % 4 == 0) or ((i + 1) == len(W))):
                self.convs.append(nn.Conv1d(L, L, 1))
        self.conv_last1 = nn.Conv1d(L, 3, 1)
        self.conv_last2 = nn.Conv1d(L, 3, 1)

    def forward(self, x):
        x = x.permute(0, 2, 1)
        conv = self.conv1(x)
        skip = self.skip(conv)
        j = 0
        for i in range(len(W)):
            conv = self.resblocks[i](conv)
            if (((i + 1) % 4 == 0) or ((i + 1) == len(W))):
                dense = self.convs[j](conv)
                j += 1
                skip = skip + dense
        CL = 2 * np.sum(AR * (W - 1))
        skip = F.pad(skip, (-CL // 2, -CL // 2))
        sites = F.softmax(self.conv_last1(skip), dim=1).permute(0, 2, 1)
        psi = F.softmax(self.conv_last2(skip), dim=1).permute(0, 2, 1)

        return sites, psi


class MainModel(nn.Module):
    def __init__(self, L, W, AR, dropout=None, EL=10000, optimizer=None, scheduler=None):
        super().__init__()
        assert len(W) == len(AR)
        self.encode = nn.parallel.DataParallel(Pangolin(L, W, AR).cuda())

        assert optimizer
        self.runoptimizer = optimizer(self.parameters())
        self.optimizer = [self.runoptimizer]

        if scheduler is not None:
            self.runscheduler = scheduler(self.runoptimizer)
            self.scheduler = [self.runscheduler]
        else:
            self.scheduler = None

        self.exp_loss = cross_entropy()

        logger.info("the number of parameters is {}".format(sum(p.numel() for p in self.encode.parameters())))

    def train_exp_step(self, d):
        seq, true_expression, sites = d["X"].float(), d["single_pred_psi"].float(), d["single_pred_psi"].float()

        self.train()
        self.runoptimizer.zero_grad()

        true_expression[torch.isnan(true_expression)] = 0
        sites[torch.isnan(sites)] = 1.
        sites = (sites > 0).float()
        sites[:, :, 0] = 1-sites[:, :, 1:].sum(-1)

        psite, psi = self.encode(seq.cuda())

        sloss = self.exp_loss(sites.cuda(), psite)
        uloss = self.exp_loss(true_expression.cuda(), psi)

        loss = uloss+sloss
        loss.backward()
        self.runoptimizer.step()

        return uloss, sloss

    def train_onestep(self, d):
        uloss,  sloss = self.train_exp_step(d)
        return {"uloss": uloss.item(), "sloss": sloss.item()}

    def get_eval_res(self, x1):
        ret = self.encode(x1)[1]
        ret = ret.detach().cpu().numpy()
        return ret

    def predict(self, d):
        self.eval()
        with torch.no_grad():
            X = d["X"].cuda().float()
            Ret = self.get_eval_res(X)
            ret = {"single_pred_psi": Ret}
            return ret
