#!/usr/bin/env python3
import torch
import torch.nn as nn
import torch.optim as optim
from transformers import Adafactor
import numpy as np
from torch.autograd import grad
from torch.utils.data import DataLoader, TensorDataset, DistributedSampler
import os
from tqdm import tqdm
from torch.optim.lr_scheduler import StepLR
from torch.cuda.amp import autocast, GradScaler
import torch.distributed as dist
import torch.multiprocessing as mp
import argparse
from torch.nn.parallel import DistributedDataParallel as DDP
import logging

def ddp_setup():
    dist.init_process_group(backend="nccl")
    torch.cuda.set_device(int(os.environ["LOCAL_RANK"]))

def cleanup():
    dist.destroy_process_group()
    
# weight_history = []

class TransformerWithPhysics(nn.Module):
    def __init__(self, input_dim, output_dim, hidden_dim, num_layers, num_heads, dropout=0.1):
        super(TransformerWithPhysics, self).__init__()
        self.embedding = nn.Linear(input_dim, hidden_dim)
        encoder_layer = nn.TransformerEncoderLayer(d_model=hidden_dim, nhead=num_heads,\
                                                   dropout=dropout)
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        self.fc_out = nn.Linear(hidden_dim, output_dim)
        
#     def get_weights(self):
#         weights = {}

#         # 获取嵌入层的权重
#         weights['embedding_weight'] = self.embedding.weight.clone().detach().cpu().numpy()

#         # 获取 Transformer Encoder 层的权重
#         transformer_weights = []
#         for i, layer in enumerate(self.transformer.layers):
#             transformer_weights.append({
#                 'self_attn_in_proj_weight': layer.self_attn.in_proj_weight.clone().detach().cpu().numpy(),
#                 'self_attn_out_proj_weight': layer.self_attn.out_proj.weight.clone().detach().cpu().numpy(),
#                 'linear1_weight': layer.linear1.weight.clone().detach().cpu().numpy(),
#                 'linear2_weight': layer.linear2.weight.clone().detach().cpu().numpy()
#             })
#         weights['transformer_weights'] = transformer_weights

#         # 获取输出层的权重
#         weights['fc_out_weight'] = self.fc_out.weight.clone().detach().cpu().numpy()

#         return weights


    def forward(self, x):
        x = self.embedding(x)
        x = x.permute(1, 0, 2)  # Transformer expects (seq_len, batch_size, hidden_dim)
        x = self.transformer(x)
        x = x.permute(1, 0, 2)  # Convert back to (batch_size, seq_len, hidden_dim)
        x = self.fc_out(x)
        return x

    def physics_loss(self, x, t, Y):
        input_seq = torch.cat((t, x), dim=2)
        
        Y_pred = self.forward(input_seq)
        
        if Y.dim() == 2:
            Y = Y.unsqueeze(1)
        
        dY_t_pred = []
        dY_x_pred = []

        for m in range(Y_pred.shape[2]):
            y_pred = Y_pred[:, :, m].reshape(-1)
            dy_t_pred = grad(
                y_pred, t, 
                grad_outputs=torch.ones_like(y_pred),
                create_graph=True,
                retain_graph=True
            )[0].reshape(Y_pred.shape[0], Y_pred.shape[1], -1)
            dy_x_pred = grad(
                y_pred, x, 
                grad_outputs=torch.ones_like(y_pred),
                create_graph=True,
                retain_graph=True
            )[0].reshape(Y_pred.shape[0], Y_pred.shape[1], -1)
            dY_t_pred.append(dy_t_pred)
            dY_x_pred.append(dy_x_pred)
        
        loss_eq = 0
        for batch in range(Y_pred.shape[0]):
            for seq in range(Y_pred.shape[1]):
                Y_pred_seq = Y_pred[batch, seq, :]
                dY_t_pred_seq = torch.cat([dY_t_pred[m][batch, seq, :].unsqueeze(0) \
                                           for m in range(Y_pred.shape[2])],\
                                          dim=0)
                dY_x_pred_seq = torch.cat([dY_x_pred[m][batch, seq, :].unsqueeze(0) \
                                           for m in range(Y_pred.shape[2])],\
                                          dim=0)
                
                n_pred, u_pred, p_pred, q_pred, E_pred = Y_pred_seq[0], Y_pred_seq[1],\
                Y_pred_seq[2], Y_pred_seq[3], Y_pred_seq[4]
                dn_t_pred, du_t_pred, dp_t_pred, dq_t_pred, dE_t_pred = dY_t_pred_seq[0],\
                dY_t_pred_seq[1], dY_t_pred_seq[2], dY_t_pred_seq[3], dY_t_pred_seq[4]
                dn_x_pred, du_x_pred, dp_x_pred, dq_x_pred, dE_x_pred = dY_x_pred_seq[0],\
                dY_x_pred_seq[1], dY_x_pred_seq[2], dY_x_pred_seq[3], dY_x_pred_seq[4]
                
                e1 = dn_t_pred + u_pred * dn_x_pred + n_pred * du_x_pred
                e2 = du_t_pred + u_pred * du_x_pred + (1 / n_pred) * dp_x_pred + E_pred
                e3 = dp_t_pred + u_pred * dp_x_pred + 3 * p_pred * du_x_pred + dq_x_pred
                e4 = dE_x_pred - (1.0 - n_pred)

                loss_eq += nn.MSELoss()(e1, torch.zeros_like(e1)) + \
                nn.MSELoss()(e2, torch.zeros_like(e2)) + \
                nn.MSELoss()(e3, torch.zeros_like(e3)) + \
                nn.MSELoss()(e4, torch.zeros_like(e4))
        
        loss_eq /= (Y_pred.shape[0] * Y_pred.shape[1] * Y_pred.shape[2])
        
        loss_mse = nn.MSELoss()(Y_pred[:, :, 3], Y[:, :, 3]) + nn.MSELoss()(Y_pred[:, :, 2], Y[:, :, 2])
        loss = loss_eq + loss_mse
        
        return loss

    
logging.basicConfig(filename='/lustre/home/sztu_pg_chenyue/jupyter/TwoStream/logfile_pq.log', level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
scaler = GradScaler()
class distribued_trainer:
    def __init__(
            self,
            model: TransformerWithPhysics,
            train_data: DataLoader,
            optimizer: optim.Optimizer,
            save_every: int,
            snapshot_path: str,
    ) -> None:
        self.gpu_id = int(os.environ["LOCAL_RANK"])
        input_dim = 2
        hidden_dim = 128
        output_dim = 5
        num_layers = 4
        num_heads = 8
        self.model = TransformerWithPhysics(input_dim, output_dim, hidden_dim, \
                                            num_layers, num_heads) 
        self.model = model.to(self.gpu_id)
        self.train_data = train_data
        self.optimizer = optimizer
        self.save_every = save_every
        self.epochs_run = 0
        self.snapshot_path = snapshot_path
        self.loss_values = [] 
        self.epoch_col = []
        if os.path.exists(snapshot_path):
            print("Loading snapshot")
            self._load_snapshot(snapshot_path)

        self.model = DDP(self.model, device_ids=[self.gpu_id])

    def _load_snapshot(self, snapshot_path):
        snapshot = torch.load(snapshot_path)
        model_state_dict = snapshot["model_state_dict"]
    
        # Remove 'module.' prefix if it exists
        if 'module.' in next(iter(model_state_dict.keys())):
            model_state_dict = {k.replace('module.', ''): v for k, v in model_state_dict.items()}
    
        self.model.load_state_dict(model_state_dict)
        self.optimizer.load_state_dict(snapshot["optimizer_state_dict"])
        self.epochs_run = snapshot["epoch"]
        print(f"Resuming training from snapshot at Epoch {self.epochs_run}")
        return model_state_dict
    
   

    def _run_batch(self, x, t, Y):
        self.optimizer.zero_grad()
        with autocast():
            loss = self.model.module.physics_loss(x, t, Y)
        scaler.scale(loss).backward()
        scaler.step(self.optimizer)
        scaler.update()
        
        return loss

    def _run_epoch(self, epoch):
        b_sz = len(next(iter(self.train_data))[0])
        print(f"[GPU{self.gpu_id}] Epoch {epoch} | Batchsize: {b_sz} | Steps: {len(self.train_data)}")
        self.train_data.sampler.set_epoch(epoch)
        epoch_loss = 0.0
        for batch in self.train_data:
            t_batch, x_batch, y_batch = [data.to(self.gpu_id) for data in batch]
            t_batch.requires_grad = True
            x_batch.requires_grad = True
            batch_loss = self._run_batch(x_batch, t_batch, y_batch)
#             print(f"Batch loss: {batch_loss}")  # 日志记录
        
#             if batch_loss is None:
#                 raise ValueError("Batch loss is None")
        
#             epoch_loss += batch_loss
#         return epoch_loss
            epoch_loss += batch_loss
        avg_epoch_loss = epoch_loss / len(self.train_data)
        
        return avg_epoch_loss

    def _save_snapshot(self, epoch):
        os.makedirs(os.path.dirname(self.snapshot_path), exist_ok=True)
        snapshot = {
            'epoch': epoch,
            'model_state_dict': self.model.module.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict() 
        }
#         with torch.no_grad():
#             weight_history.append(self.model.module.get_weights())
        torch.save(snapshot, self.snapshot_path)
        print(f"Epoch {epoch} | Training snapshot saved at {self.snapshot_path}")
        
#         if epoch % 10 == 0:
#             weight_file_path = os.path.join(os.path.dirname(self.snapshot_path), f'weights_epoch_{epoch}.npz')
#             np.savez(weight_file_path, **self.model.module.get_weights())
#             print(f"Epoch {epoch} | Model weights saved at {weight_file_path}")

    def train(self, max_epochs: int, loss_threshold: float = 1e-7):
        for epoch in range(self.epochs_run, max_epochs):
            avg_epoch_loss = self._run_epoch(epoch)
            logging.info(f'Epoch: {epoch}, Loss: {avg_epoch_loss}')
            
            if avg_epoch_loss < loss_threshold:
                print(f"Loss has reached the threshold of {loss_threshold}. Stopping training.")
                break
            
            if self.gpu_id == 2 and epoch % self.save_every == 0:
                self._save_snapshot(epoch)
            if epoch % 10 == 0:
                self.loss_values.append(avg_epoch_loss)
                self.epoch_col.append(epoch)
        np.savez('/lustre/home/sztu_pg_chenyue/jupyter/TwoStream/models/ld_tf_ddp_3/pq/loss_values.npz', loss=self.loss_values, epoch=self.epoch_col)
        print(f'All loss values saved to loss_values.npz')
            

def prepare_dataloader(dataset: TensorDataset, batch_size: int):
    return DataLoader(
        dataset,
        batch_size=batch_size,
        pin_memory=True,
        shuffle=False,
        sampler=DistributedSampler(dataset)
    )

def main(save_every: int, total_epochs: int, batch_size: int, snapshot_path: str = "/lustre/home/sztu_pg_chenyue/jupyter/TwoStream/models/ld_tf_ddp_3/q/snapshot.pt"):
    ddp_setup()
    
    # 加载数据
    data = np.load('/lustre/home/sztu_pg_chenyue/jupyter/TwoStream/simulation/data_LD_6.npz')
    X, T = data["X"].flatten(), data["T"].flatten()
    Y_full = np.column_stack([(data[key]).flatten() for key in "nupqE"])

    ind = np.random.randint(0,X.shape,size=4000)
    X_train = X[ind]
    T_train = T[ind]
    Y_train = Y_full[ind]

    X_train = torch.tensor(X_train[:, None], dtype=torch.float32)
    T_train = torch.tensor(T_train[:, None], dtype=torch.float32)
#     p_train = torch.tensor(Y_train[:, 2], dtype=torch.float32)
    Y_train = torch.tensor(Y_train, dtype=torch.float32)
    
    train_set1 = TensorDataset(T_train.unsqueeze(-1), X_train.unsqueeze(-1), Y_train.unsqueeze(1))  # load your dataset
    input_dim = 2
    hidden_dim = 128
    output_dim = 5
    num_layers = 4
    num_heads = 8
    model = TransformerWithPhysics(input_dim, output_dim, hidden_dim, num_layers, num_heads)  # load your model
    lr = 0.001
#     optimizer = Adafactor(
#         model.parameters(),
#         lr=1e-3,
#         scale_parameter=True,
#         relative_step=False,  # 手动控制学习率时设置为 False
#         warmup_init=False
#      )
#     scheduler = AdafactorSchedule(optimizer)
    optimizer = optim.AdamW(model.parameters(), lr=lr) # 使用Transformers提供的AdaW优化器，实现了梯度偏差校正和权重衰减
    
   
    batch_size = 64
    train_data = prepare_dataloader(train_set1, batch_size)
    trainer = distribued_trainer(model, train_data, optimizer, save_every, snapshot_path)
    trainer.train(total_epochs)
    cleanup()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='simple distributed training job')
    parser.add_argument('total_epochs', type=int, help='Total epochs to train the model')
    parser.add_argument('save_every', type=int, help='How often to save a snapshot')
    parser.add_argument('--batch_size', default=64, type=int, help='Input batch size on each device (default: 32)')
    args = parser.parse_args()
    
    main(args.save_every, args.total_epochs, args.batch_size)
