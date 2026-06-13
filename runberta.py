import warnings
warnings.filterwarnings('ignore')
from transformers import ( 
        RobertaTokenizer, 
        RoFormerModel,
    )
import pandas as pd
import argparse
import torch
import os
import numpy as np
from torch.utils.data import Dataset, DataLoader
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import optim
import argparse
import random
import gc
import CoMBCR

def seed_torch(seed=0):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

def get_package_abspath():
    # Get the directory path of the package
    package_path = os.path.dirname(CoMBCR.__file__)
    # Get the absolute path
    abs_path = os.path.abspath(package_path)
    return abs_path

# Print the absolute path of the CoMBCR package
print(get_package_abspath())
    
def batch_process_sequences(sequences, tokenizer, model, device, batch_size=32):
    berta_embs = []
    
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i + batch_size]
        tokenized_input = tokenizer(batch, return_tensors='pt', padding=True, truncation=True)
        
        # Move tensors to GPU if available
        tokenized_input = {key: val.to(device) for key, val in tokenized_input.items()}
        
        with torch.no_grad():
            emb = model(**tokenized_input).last_hidden_state
            attention_mask = tokenized_input['attention_mask']
            
            # Calculate the actual lengths of sequences
            sequence_lengths = attention_mask.sum(dim=1).cpu().numpy() - 2  # Exclude [CLS] and [SEP]
            
            for j, seq_len in enumerate(sequence_lengths):
                # Collect embeddings excluding [CLS] and [SEP] tokens for each sequence
                berta_embs.append(emb[j, 1:seq_len+1, :].mean(0).cpu().numpy().reshape(1,-1))
    
    # Concatenate results from all batches
    berta_emb = np.concatenate(berta_embs, axis=0)
    
    return berta_emb

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--datapath', type=str)
    parser.add_argument('--outdir', type=str, default="./")
    parser.add_argument('--outfilename', type=str, default="antiberta_embedding.csv")
    args = parser.parse_args()
    
    datapath = args.datapath
    outdir = args.outdir
    outfilename = args.outfilename
    
    seed_torch()
    # Find the directory of the current file (combcr.py)
    current_dir = get_package_abspath()
     # Check for GPU availability
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    tokenizer = RobertaTokenizer.from_pretrained(os.path.join(current_dir, "tokenizer"), max_len=150)
    model = RoFormerModel.from_pretrained(os.path.join(current_dir, "BCRencoder")).to(device)
    
    bcr_file = pd.read_csv(datapath, index_col="barcode")
    bcr_file["whole_seq"] = bcr_file["fwr1"].str.cat([bcr_file["cdr1"], bcr_file["fwr2"], bcr_file["cdr2"], bcr_file["fwr3"], bcr_file["cdr3"], bcr_file["fwr4"]])
    
    model.to(device)
    # Process sequences in batches
    berta_emb  = batch_process_sequences(bcr_file["whole_seq"].tolist(), tokenizer, model, device=device)
    del model
    gc.collect()
    torch.cuda.empty_cache()

    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    assert(berta_emb.shape[0] == bcr_file.shape[0])
    berta_emb_df = pd.DataFrame(berta_emb)
    berta_emb_df['barcode'] = bcr_file.index.tolist()
    berta_emb_df.to_csv(os.path.join(outdir, outfilename), index=False)
    
    
   