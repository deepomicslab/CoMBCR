import warnings
warnings.filterwarnings('ignore')
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim.lr_scheduler import ReduceLROnPlateau, MultiStepLR
from torch import optim
from torch.utils.data import Dataset, DataLoader
import os
import numpy as np
import re
import pickle
import pandas as pd
import random
import csv
from sklearn.preprocessing import LabelEncoder
from transformers import ( 
        RobertaTokenizer, 
        RoFormerModel,
        pipeline
    )

def seed_torch(seed=0):
	random.seed(seed)
	os.environ['PYTHONHASHSEED'] = str(seed) # 为了禁止hash随机化，使得实验可复现
	np.random.seed(seed)
	torch.manual_seed(seed)
	torch.cuda.manual_seed(seed)
	torch.cuda.manual_seed_all(seed) # if you are using multi-GPU.
	torch.backends.cudnn.benchmark = False
	torch.backends.cudnn.deterministic = True
	os.environ['CUBLAS_WORKSPACE_CONFIG']=':16:8'
	torch.use_deterministic_algorithms(True, warn_only=True)


class AverageMeter(object):
    """Computes and stores the average and current value"""

    def __init__(self):
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count

class Logger(object):

    def __init__(self, path, header):
        self.log_file = open(path, 'w')
        self.logger = csv.writer(self.log_file, delimiter='\t')

        self.logger.writerow(header)
        self.header = header

    def __del(self):
        self.log_file.close()

    def log(self, values):
        write_values = []
        for col in self.header:
            assert col in values
            write_values.append(values[col])

        self.logger.writerow(write_values)
        self.log_file.flush()
        
class EarlyStopping:
    """Early stops the training if validation loss doesn't improve after a given patience."""
    def __init__(self, save_path, checkpoint, patience=7, verbose=False, delta=0):
        """
        
            save_path : ???????
            patience (int): How long to wait after last time validation loss improved.
                            Default: 7
            verbose (bool): If True, prints a message for each validation loss improvement. 
                            Default: False
            delta (float): Minimum change in the monitored quantity to qualify as an improvement.
                            Default: 0
        """
        self.save_path = save_path
        self.checkpoint = checkpoint
        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.val_loss_min = np.Inf
        self.delta = delta

    def __call__(self, val_loss, model):

        score = -val_loss

        if self.best_score is None:
            self.best_score = score
            self.save_checkpoint(val_loss, model)
        elif score < self.best_score * (1-self.delta):
            self.counter += 1
            print(f'EarlyStopping counter: {self.counter} out of {self.patience}')
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            self.save_checkpoint(val_loss, model)
            self.counter = 0

    def save_checkpoint(self, val_loss, model):
        '''Saves model when validation loss decrease.'''
        if self.verbose:
            print(f'Validation loss decreased ({self.val_loss_min:.6f} --> {val_loss:.6f}).  Saving model ...')
        path = os.path.join(self.save_path, self.checkpoint)
        torch.save(model.state_dict(), path)	# ??????????????
        self.val_loss_min = val_loss
        
class ProjectionHead(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(ProjectionHead, self).__init__()
        self.layers = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, output_dim)
        )

    def forward(self, x):
        return self.layers(x)

class BCRencoder(nn.Module):
    def __init__(self):
        super(BCRencoder, self).__init__()
        # Find the directory of the current file (combcr.py)
        current_dir = os.path.dirname(os.path.abspath(__file__))
        self.berta = RoFormerModel.from_pretrained(os.path.join(current_dir, "BCRencoder"))
        #for p in self.parameters():
        #    p.requires_grad = False
        
    def forward(self, x):
        bertaoutput = self.berta(**x).last_hidden_state
        embs = []
        for i, tokens in enumerate(x['input_ids']):
            embs.append(bertaoutput[i, 1:len(tokens)-1].mean(0).reshape(1,-1))
        x = torch.concat(embs)
        #x = bertaoutput[:,0,:].reshape(-1, bertaoutput.shape[-1])
        return x
        

class Encoder_profile(nn.Module):
    """
    this is the profile encoder class, which contains the 3 linear transformation
    layers and 1 layer normalization
    
    Parameters:
        param in_dim: the input profile dimension
        param hid_dim: the hidden dimension of first transformation layer
        param hid_dim2: the hidden dimension of second transformation layer
        param out_dim: the output dimension
    
    Returns:
        new profile compressed embedding(bs, out_dim)
    """
    
    def __init__(self, in_dim = 5, hid_dim = 5, hid_dim2 = 5, out_dim = 5):
        super(Encoder_profile, self).__init__()
        self.layer1 = nn.Linear(in_dim, hid_dim)
        self.layer2 = nn.Linear(hid_dim, hid_dim2)
        self.layer3 = nn.Linear(hid_dim2, out_dim)
        
        # define the tanh as activation function
        self.activation = nn.ReLU()
        
        # define the layer normalization
        self.layer_norm = nn.LayerNorm(out_dim)
        
    def forward(self, profile):
        compressed_profile = self.activation(self.layer1(profile))
        compressed_profile = self.activation(self.layer2(compressed_profile))
        compressed_profile = self.activation(self.layer3(compressed_profile))
        compressed_profile = self.layer_norm(compressed_profile)
        return compressed_profile

class CustomDataset(Dataset):
    def __init__(self, bcrfile_path, rnafile_path, separatebatch=False):
        super(CustomDataset, self).__init__()
        rnafile = pd.read_csv(rnafile_path, sep=",", index_col="barcode", low_memory=False)
        self.rnamatrix = torch.from_numpy(rnafile.values).float()
        
        bcrfile = pd.read_csv(bcrfile_path, sep=",", index_col="barcode", low_memory=False)
        le = LabelEncoder()
        bcrfile = bcrfile[(bcrfile["chain"] == "IGH")]
        bcrfile["whole_seq"] = bcrfile["fwr1"].str.cat([bcrfile["cdr1"], bcrfile["fwr2"], bcrfile["cdr2"], bcrfile["fwr3"], bcrfile["cdr3"], bcrfile["fwr4"]])
        bcrfile['sample'] = bcrfile['sample'].astype(str)
        if separatebatch:
            print("keep invariant of a batch")
            bcrfile["new_whole_seq"] = bcrfile["sample"].str.cat(bcrfile["whole_seq"], sep=":")
            self.bcridentity = le.fit_transform(bcrfile.new_whole_seq.tolist()) #map cdr3 into integer labels
        else:
            self.bcridentity = le.fit_transform(bcrfile.whole_seq.tolist())
        self.bcridentity = torch.from_numpy(self.bcridentity)
        self.bcr = bcrfile.whole_seq.tolist()
        
        assert(len(rnafile)==len(bcrfile))
    
    def __len__(self):
        return self.rnamatrix.shape[0]
    
    def __getitem__(self, idx):
        identity = self.bcridentity[idx]
        gex = self.rnamatrix[idx]
        bcr = self.bcr[idx]
        
        return bcr, gex, identity, idx

class CoMBCR_model(nn.Module):
    def __init__(self, encoderBCR_out_dim,
                 encoderprofile_in_dim, encoderprofile_hid_dim, encoderprofile_hid_dim2, encoderprofle_out_dim):
        super(CoMBCR_model, self).__init__()
        
        # define the BCRencoder
        
        self.encoder_BCR = BCRencoder()
        self.BCR_project = ProjectionHead(1024, encoderBCR_out_dim, encoderBCR_out_dim)
        # define the temperature (default set as 1)
        self.encoder_profile = Encoder_profile(encoderprofile_in_dim, encoderprofile_hid_dim, 
                                              encoderprofile_hid_dim2, encoderprofle_out_dim)
        # define the project
        self.profile_proj = ProjectionHead(encoderprofle_out_dim, encoderprofle_out_dim, encoderprofle_out_dim)
        
        
        
    def forward(self, bcr_token, exp):
             
        # extract the features of BCR
        encoderBCR_feature = self.encoder_BCR(bcr_token)
        encoderBCR_embedding = self.BCR_project(encoderBCR_feature)
        # normalization of BCR sequence embedding
        encoderBCR_embedding = F.normalize(encoderBCR_embedding, dim = -1)
        
        # extract the profile features based on the encoder_profile
        encoderprofile_feature = self.encoder_profile(exp)
        encoderprofile_embedding = self.profile_proj(encoderprofile_feature)
        # normalization of expression embedding
        encoderprofile_embedding = F.normalize(encoderprofile_embedding, dim = -1)
        
        return encoderBCR_embedding, encoderprofile_embedding
    
    
class mydefine_loss(nn.Module):
    def __init__(self, temperature=0.2, lam=0.1):
        super().__init__()
        self.temperature = temperature
        self.lam = lam
        
    def forward(self, encoderBCR_embedding, encoderprofile_embedding, batch_BCRsimilar, batch_expdist, idx_BCR):
        '''
        output: graphs.num x emb
        batch_BCRdist/batch_expdist: batch_size x batch_size
        '''
        
        # pMHC encoder and TCR encoder similarity score calculation
        sim_p2b = encoderprofile_embedding @ encoderBCR_embedding.T / self.temperature
        sim_b2p = encoderBCR_embedding @ encoderprofile_embedding.T / self.temperature
        
        # same TCR with mulitiple differently profile
        with torch.no_grad():
            
            idx_BCR = idx_BCR.view(-1, 1)
            pos_idx_BCR = torch.eq(idx_BCR, idx_BCR.T)
            
            sim_targets = (pos_idx_BCR).float().to(encoderBCR_embedding.device)
            sim_targets = sim_targets / sim_targets.sum(1, keepdim = True)
        #p2b loss
        loss_p2t = -torch.sum(F.log_softmax(sim_p2b, dim = 1) * sim_targets, dim = 1).mean()
        loss_b2p = -torch.sum(F.log_softmax(sim_b2p, dim = 1) * sim_targets, dim = 1).mean()
        loss_cmc = (loss_p2t + loss_b2p) / 2
        
        #b2b loss
        sim_b2b = encoderBCR_embedding @ encoderBCR_embedding.T
        loss_b2b = ((sim_b2b - batch_BCRsimilar)**2).mean()
        #loss_b2b = -torch.sum(F.log_softmax(1 - sim_b2b, dim = 1) * torch.softmax(batch_BCRsimilar, dim = 1), dim = 1).mean()
        
        # profile similarity score calculation
        sim_p2p = encoderprofile_embedding @ encoderprofile_embedding.T
        loss_p2p_inner = -torch.sum(F.log_softmax(sim_p2p, dim = 1) * sim_targets, dim = 1).mean()
        #p2p loss
        loss_p2p = -torch.sum(F.log_softmax(1 - sim_p2p, dim = 1) * torch.softmax(batch_expdist, dim = 1), dim = 1).mean() + self.lam*loss_p2p_inner
        
        return loss_cmc, loss_b2b, loss_p2p 


def CoMBCR_main(bcrpath, rnapath, bcroriginal, outdir, checkpoint="best_network.pth", 
                lr=1e-6, lam=1e-1, batch_size=256, epochs=200, patience=15, lr_step=[40,100], encoderprofile_in_dim=5000, separatebatch = False):
    
    
    # Find the directory of the current file (combcr.py)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    if not os.path.exists(os.path.join(outdir)):
        os.mkdir(os.path.join(outdir))
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    print('learning rate is ', lr)
    bcrfile_path = bcrpath
    rnafile_path = rnapath
    dataset = CustomDataset(bcrfile_path = bcrfile_path, rnafile_path=rnafile_path, separatebatch=separatebatch)
    
    # prepare intra-distance
    bcrsimilar = pd.read_csv(bcroriginal, index_col="barcode", low_memory=False).values
    assert(bcrsimilar.shape[0] == len(dataset))
    bcrsimilar = torch.from_numpy(bcrsimilar)
    bcrsimilar = F.normalize(bcrsimilar, dim = -1)
    bcrsimilar = bcrsimilar @ bcrsimilar.T
    rnavalues = pd.read_csv(rnafile_path, sep=",", index_col="barcode").values
    rnavalues = torch.from_numpy(rnavalues)
    rnadist = torch.cdist(rnavalues, rnavalues)
    
    seed_torch()
    tokenizer = RobertaTokenizer.from_pretrained(os.path.join(current_dir, "tokenizer"), max_len=150)
    early_stopping = EarlyStopping(outdir, patience=patience, checkpoint=checkpoint)
    loader = DataLoader(dataset, batch_size=batch_size,shuffle=True)
    runmodel = CoMBCR_model(encoderBCR_out_dim=256,
                     encoderprofile_in_dim=encoderprofile_in_dim, encoderprofile_hid_dim=1024, 
                     encoderprofile_hid_dim2=512, encoderprofle_out_dim=256)
    loss_cmc_epoch = AverageMeter()
    loss_b2b_epoch = AverageMeter()
    loss_p2p_epoch = AverageMeter()
    loss_all = AverageMeter()
    logger = Logger(os.path.join(outdir, "{}.log".format(checkpoint)), ['Epoch', 'loss', 'loss_cmc', 'loss_p2p', 'loss_b2b'])
    calculate_loss = mydefine_loss(temperature=0.2, lam=lam)
    calculate_loss.to(device)
    runmodel.to(device)
    #optimizer=optim.Adam(filter(lambda p: p.requires_grad, runmodel.parameters()), lr=lr)
    optimizer=optim.Adam(runmodel.parameters(), lr=lr)
    #scheduler = ReduceLROnPlateau(optimizer, patience=10, verbose=True)
    scheduler = MultiStepLR(optimizer, milestones=lr_step, verbose=True)
    epochs = epochs
    seed_torch()
    
    runmodel.train()
    for epoch in np.arange(epochs):
        loss_cmc_epoch.reset()
        loss_b2b_epoch.reset()
        loss_p2p_epoch.reset()
        loss_all.reset()
        for batch, (bcr, gex, identity, idx) in enumerate(loader):
            tokenized_input = tokenizer(list(bcr), return_tensors='pt', padding=True)
            tokenized_input = tokenized_input.to(device)
            gex = gex.to(device)
            identity_BCR = identity.to(device)
            idxs = idx
            optimizer.zero_grad()
           
            encoderBCR_embedding, encoderprofile_embedding = runmodel(tokenized_input, gex)
            batch_bcr_simiar = bcrsimilar[:, idxs][idxs,:].to(device)
            batch_gex_dist = rnadist[:, idxs][idxs,:].to(device)
            loss_cmc, loss_b2b, loss_p2p = calculate_loss(encoderBCR_embedding, encoderprofile_embedding, batch_bcr_simiar, batch_gex_dist, identity_BCR)
            loss = loss_cmc + loss_p2p/(loss_p2p.detach()/loss_cmc.detach()) + loss_b2b/(loss_b2b.detach()/loss_cmc.detach())
            #loss = loss_cmc/(loss_cmc.detach()) + loss_p2p/loss_p2p.detach() + loss_b2b/loss_b2b.detach()
            loss.backward()
            optimizer.step()
            
            loss_cmc_epoch.update(loss_cmc.item(), gex.shape[0]*gex.shape[0]) #embs.shape[0]*(embs.shape[0]-1)/2
            loss_b2b_epoch.update(loss_b2b.item(), gex.shape[0]*gex.shape[0])
            loss_p2p_epoch.update(loss_p2p.item(), gex.shape[0]*gex.shape[0])
            loss_all.update(loss_cmc.item()+loss_b2b.item()+loss_p2p.item(), 1)
            #loss_all.update(loss.item(), 1)

            # if batch %10 == 0:
            #     print('steps:{}\tloss:{:.5f}\tloss_cmc:{:.5f}\tloss_p2p:{:.5f}\tloss_b2b:{:.5f}'.format(batch, loss, loss_cmc, loss_p2p, loss_b2b))
            torch.cuda.empty_cache()
        early_stopping(loss_all.avg, runmodel)
        #scheduler.step(loss_all.avg)
        scheduler.step()
        if early_stopping.early_stop:
            print("Early stopping")
            break
        logger.log({'Epoch':epoch, 'loss': round(loss_all.avg, 6), 'loss_cmc':round(loss_cmc_epoch.avg, 6), 'loss_p2p':round(loss_p2p_epoch.avg, 6), 'loss_b2b': round(loss_b2b_epoch.avg, 6)})
        if epoch %1 ==0:
            #print(embs)
            print('Epoch:[{}/{}]\tloss:{:.5f}\tloss_cmc:{:.6f}\tloss_p2p:{:.6f}\tloss_b2b:{:.6f}'.format(epoch, epochs, loss_all.avg, loss_cmc_epoch.avg, loss_p2p_epoch.avg, loss_b2b_epoch.avg))

    runmodel = CoMBCR_model(encoderBCR_out_dim=256,
                 encoderprofile_in_dim=encoderprofile_in_dim, encoderprofile_hid_dim=1024, 
                 encoderprofile_hid_dim2=512, encoderprofle_out_dim=256)
    loader = DataLoader(dataset, batch_size=256,shuffle=False)
    runmodel.load_state_dict(torch.load(os.path.join(outdir, checkpoint), weights_only=True))
    runmodel.eval()
    runmodel.to(device)
    BCRallembs = []
    gexallembs = []
    for epoch in np.arange(1):
        with torch.no_grad():
            for batch, (bcr, gex, identity, idx) in enumerate(loader):
                tokenized_input = tokenizer(list(bcr), return_tensors='pt', padding=True)
                tokenized_input = tokenized_input.to(device)
                encoderBCR_embedding = runmodel.encoder_BCR.eval()(tokenized_input)
                encoderBCR_embedding = runmodel.BCR_project.eval()(encoderBCR_embedding)
                # normalization of BCR sequence embedding
                encoderBCR_embedding = F.normalize(encoderBCR_embedding, dim = -1)
                BCRallembs.append(encoderBCR_embedding.cpu().numpy())
                
                gex = gex.to(device)
                encoderprofile_embedding = runmodel.encoder_profile.eval()(gex)
                encoderprofile_embedding = runmodel.profile_proj.eval()(encoderprofile_embedding)
                # normalization of BCR sequence embedding
                encoderprofile_embedding = F.normalize(encoderprofile_embedding, dim = -1)
                gexallembs.append(encoderprofile_embedding.cpu().numpy())
    BCRembeddings = np.concatenate(BCRallembs, axis=0)
    gexembeddings = np.concatenate(gexallembs, axis=0)
    
    if not os.path.exists(os.path.join(outdir, "Embeddings")):
        os.mkdir(os.path.join(outdir, "Embeddings"))
    pd.DataFrame(BCRembeddings).to_csv(os.path.join(outdir, "Embeddings", "bcrembeddings.csv"), index=False) 
    pd.DataFrame(gexembeddings).to_csv(os.path.join(outdir, "Embeddings", "gexembeddings.csv"), index=False) 
    
    return BCRembeddings, gexembeddings

 

