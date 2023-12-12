"""
This file has the functions necessary to create the SpliceAI and related resnet models.
"""
import torch
import torch.nn as nn
from functools import partial
from .model_utils import cross_entropy
from loguru import logger


def residualunit(l, w, ar, dropout=None, norm="BN"):
    w = int(w)
    ar = int(ar)
    if norm.lower() == "bn":
        Norm = partial(nn.BatchNorm1d, l)
    else:
        print("GN is used")
        Norm = partial(nn.GroupNorm, num_groups=1, num_channels=l)
    if dropout:
        return nn.Sequential(
            *[
                Norm(),
                nn.ReLU(),
                nn.Conv1d(in_channels=l, out_channels=l,
                          kernel_size=w, dilation=ar),
                Norm(),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Conv1d(in_channels=l, out_channels=l,
                          kernel_size=w, dilation=ar),
            ]
        )
    else:
        return nn.Sequential(
            *[

                Norm(),
                nn.ReLU(),
                nn.Conv1d(in_channels=l, out_channels=l,
                          kernel_size=w, dilation=ar),
                Norm(),
                nn.ReLU(),
                nn.Conv1d(in_channels=l, out_channels=l,
                          kernel_size=w, dilation=ar),
            ]
        )


class ResidualUnit(nn.Module):
    def __init__(self, base_unit):
        super().__init__()
        self.net = base_unit()

    def forward(self, x):
        px = self.net(x)
        pl = px.shape[-1]
        xl = x.shape[-1]
        assert (xl-pl) % 2 == 0
        x = x[:, :, (xl-pl)//2:pl+(xl-pl)//2]+px
        return x


class Encode(nn.Module):
    def __init__(self, EL, Ls, Ws, ARs, dropout):
        super().__init__()
        net = []
        self.exlength = 0
        for i in range(len(Ls)):
            base_unit = partial(
                residualunit, l=Ls[i], w=Ws[i], ar=ARs[i], dropout=dropout, norm="BN")
            net.append(ResidualUnit(base_unit))
            self.exlength += Ws[i]//2*ARs[i]*2

        self.encodenet = nn.Sequential(*net)
        self.EL = EL
        self.targetconv1d = nn.Conv1d(4, Ls[0], 1)
        self.bn = nn.BatchNorm1d(Ls[-1], affine=False)
        self.usagelinear = nn.Sequential(
            nn.Linear(2+Ls[-1]*2, 256), nn.ReLU(), nn.Linear(256, Ls[-1]))
        self.out_usage_net = nn.Sequential(nn.Linear(
            Ls[-1], 256), nn.ReLU(), nn.Linear(256, 256), nn.ReLU(), nn.Linear(256, 3))
        self.out_delta_net = nn.Sequential(nn.Linear(
            Ls[-1], 256), nn.ReLU(), nn.Linear(256, 256), nn.ReLU(), nn.Linear(256, 3))
        self.out_site_net = nn.Sequential(nn.Linear(
            Ls[-1], 256), nn.ReLU(), nn.Linear(256, 256), nn.ReLU(), nn.Linear(256, 3))
        self.soft = nn.Softmax(-1)

    def pred_delta(self, mut, ref, refusage):
        reffeature = self.usagelinear(
            torch.cat([torch.log(refusage[:, :, 1:].clamp(1e-10)), mut, ref], -1))
        ret = self.soft(self.out_delta_net(mut-ref+reffeature))
        return ret

    def forward(self, targetseq, seq=None, exp=None,  targetusage=None, rusage=None):
        length = targetseq.shape[1]
        bs = len(targetseq)
        outlength = length-self.EL

        if targetusage is None:  # eval mode
            targetfeature = self.targetconv1d(targetseq.permute(0, 2, 1))
            if seq is not None:
                seqfeature = self.targetconv1d(seq.permute(0, 2, 1))
                targetfeature = torch.cat([targetfeature, seqfeature], 0)

            targetfeature = self.bn(self.encodenet(
                targetfeature)).permute(0, 2, 1)
            padsize = targetfeature.shape[1]-outlength
            assert padsize % 2 == 0
            padsize = padsize//2
            targetfeature = targetfeature[:, padsize:padsize+outlength]

            # ref usage values is None
            if exp is None or torch.isnan(exp.sum()):
                if seq is None:
                    return self.soft(self.out_usage_net(targetfeature[:bs]))
                else:
                    exp = self.soft(self.out_usage_net(targetfeature[bs:]))
                    return self.pred_delta(targetfeature[:bs], targetfeature[bs:], exp)

            else:
                return self.pred_delta(targetfeature[:bs], targetfeature[bs:], exp)
        else:  # train mode, where exp is the junction sites
            length = targetseq.shape[1]
            targetfeature = self.targetconv1d(targetseq.permute(0, 2, 1))
            targetfeature = self.bn(self.encodenet(
                targetfeature)).permute(0, 2, 1)
            padsize = targetfeature.shape[1]-outlength
            assert padsize % 2 == 0
            padsize = padsize//2
            targetfeature = targetfeature[:, padsize:padsize+outlength]
            # pred psi
            psite = self.soft(self.out_site_net(targetfeature))
            pusage = self.soft(self.out_usage_net(targetfeature))

            # prepare data
            if rusage is None:
                rusage = torch.rand_like(targetusage[:, :, :1])*2-1.
                racc = rusage.clamp(0,)
                rdon = (-rusage).clamp(0,)
                rusage = torch.round(
                    torch.cat([1-racc-rdon, racc, rdon], -1), decimals=3)

            rlabel = targetusage
            rpusage = self.pred_delta(targetfeature, torch.flip(
                targetfeature, (0,)), torch.flip(targetusage.detach(), (0,)))  # random pair

            r1pusage = self.pred_delta(targetfeature, targetfeature, rusage)
            r1label = rusage
            return psite, exp, rpusage, rlabel, r1pusage, r1label, pusage, targetusage


class MainModel(nn.Module):
    def __init__(self, L, W, AR, dropout=None, EL=10000, optimizer=None, scheduler=None):
        super().__init__()
        assert len(W) == len(AR)
        self.encode = nn.parallel.DataParallel(
            Encode(EL, [L]*(len(W)), W, AR, dropout).cuda())

        assert optimizer

        self.runoptimizer = optimizer(self.parameters())
        self.optimizer = [self.runoptimizer]
        self.ref_seq, self.ref_usage = None, None

        if scheduler is not None:
            self.runscheduler = scheduler(self.runoptimizer)
            self.scheduler = [self.runscheduler]
        else:
            self.scheduler = None

        self.exp_loss = cross_entropy()
        logger.info("the number of parameters is {}".format(
            sum(p.numel() for p in self.encode.parameters())))

    def train_exp_step(self, d):
        seq, true_expression, sites = d["X"].float(
        ), d["single_pred_psi"].float(), d["single_pred_psi"].float()

        self.train()
        self.runoptimizer.zero_grad()

        true_expression[torch.isnan(true_expression)] = 0
        sites[torch.isnan(sites)] = 1.
        sites = (sites > 0).float()
        sites[:, :, 0] = 1-sites[:, :, 1:].sum(-1)
        sitepred, sitelabel, rpred, rlabel, r1pred, r1label, brainpred, brainlabel = self.encode(
            seq.cuda(), targetusage=true_expression.cuda(), exp=sites.cuda())

        sloss = self.exp_loss(sitelabel, sitepred)
        rloss = self.exp_loss(rlabel, rpred) - self.exp_loss(rlabel, rlabel)
        r1loss = self.exp_loss(r1label, r1pred)-self.exp_loss(r1label, r1label)
        brain_loss = self.exp_loss(brainlabel, brainpred)

        loss = sloss+rloss*1e-5+r1loss*1e-5+brain_loss
        loss.backward()
        self.runoptimizer.step()

        return sloss, rloss, r1loss, brain_loss

    def train_onestep(self, d):
        sloss, rloss, r1loss, brain_loss = self.train_exp_step(d)
        return {"sloss": sloss.item(),  "rloss": rloss.item(),  "r1loss": r1loss.item(), "brain_loss": brain_loss.item()}

    def get_eval_res(self, x1, x2=None, exp=None):
        ret = self.encode(x1, seq=x2, exp=exp)
        ret = ret.detach().cpu().numpy()
        return ret

    def predict(self, d):
        self.eval()
        with torch.no_grad():
            X = d["X"].cuda().float()
            if "mutX" in d:
                mutX = d["mutX"].cuda().float()
                if "single_pred_psi" in d:
                    exp = d["single_pred_psi"].cuda().float()
                    exp = torch.cat([exp, exp], 0)
                else:
                    exp = None
                Pred = self.get_eval_res(
                    torch.cat([X, mutX], 0), torch.cat([X, X], 0), exp)
                Pred, MutPred = Pred[:1], Pred[1:]

                ret = {"single_pred_psi": Pred, "mutY": MutPred}
            else:
                Ret = self.get_eval_res(X)
                ret = {"single_pred_psi": Ret}
            return ret
