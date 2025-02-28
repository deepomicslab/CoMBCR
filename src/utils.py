import torch
from transformers import RoFormerModel
import os
from huggingface_hub import snapshot_download

def download_BCRencoder():
    # Find the directory of the current file (combcr.py)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    if not os.path.exists(os.path.join(current_dir, 'BCRencoder')):
        os.mkdir(os.path.join(current_dir, 'BCRencoder'))
    snapshot_download(repo_id="alchemab/antiberta2-cssp", local_dir=os.path.join(current_dir, 'BCRencoder'))
    #model = RoFormerModel.from_pretrained("alchemab/antiberta2-cssp")
    print("Download Finished. Path {}".format(os.path.join(current_dir, 'BCRencoder')))