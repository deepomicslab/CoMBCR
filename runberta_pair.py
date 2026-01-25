import warnings
warnings.filterwarnings('ignore')
from transformers import ( 
        RobertaTokenizer,
        RoFormerTokenizer, 
        RoFormerModel,
        pipeline
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

def format_bcr_for_tokenizer(df, igh_col='IGH', igl_col='IGL'):
    """
    Convert B cell sequences to tokenizer format.
    For each B cell, combine IGH and IGL sequences with [SEP] token.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing B cell data with IGH and IGL columns
    igh_col : str
        Column name for heavy chain sequences (default: 'IGH')
    igl_col : str
        Column name for light chain sequences (default: 'IGL')
    
    Returns:
    --------
    pandas.DataFrame
        DataFrame with additional 'tokenized_seq' column
    """
    
    def format_sequence(row):
        """Format a single B cell sequence"""
        igh_seq = row[igh_col]
        igl_seq = row[igl_col]
        
        # Check if sequences are valid
        if pd.isna(igh_seq) or pd.isna(igl_seq):
            return None
        
        # Convert sequences to space-separated amino acids
        igh_formatted = ' '.join(list(igh_seq))
        igl_formatted = ' '.join(list(igl_seq))
        
        # Combine with [SEP] token
        tokenized_seq = f"{igh_formatted} [SEP] {igl_formatted}"
        
        return tokenized_seq
    
    tokenized_seqs = df.apply(format_sequence, axis=1)
    
    return list(tokenized_seqs.values)
    
def batch_process_sequences(sequences, tokenizer, model, device, batch_size=32):
    berta_embs = []
    
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i + batch_size]
        tokenized_input = tokenizer(batch, return_tensors='pt', padding='max_length', truncation=True, max_length=255)
 
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
    
    tokenizer = RoFormerTokenizer.from_pretrained(os.path.join(current_dir, "BCRencoder"))
    model = RoFormerModel.from_pretrained(os.path.join(current_dir, "BCRencoder")).to(device)
    
    bcr_file = pd.read_csv(datapath, index_col="barcode")
    # Process sequences in batches
    tokenized_sequences = format_bcr_for_tokenizer(bcr_file)
    print(f"Found {len(tokenized_sequences)} valid sequences")
    # Process sequences in batches
    model.to(device)
    berta_emb  = batch_process_sequences(tokenized_sequences, tokenizer, model, device=device)
    del model
    gc.collect()
    torch.cuda.empty_cache()

    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    assert(berta_emb.shape[0] == bcr_file.shape[0])
    berta_emb_df = pd.DataFrame(berta_emb)
    berta_emb_df['barcode'] = bcr_file.index.tolist()
    berta_emb_df.to_csv(os.path.join(outdir, outfilename), index=False)
    
    
   