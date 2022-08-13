import pdb
import pickle
import sqlite3
import numpy as np
import torch
import pytorch_lightning as pl

from contextlib import closing
from torch.utils.data import Dataset, DataLoader

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# TODO: Change this to load data from CSV or other simpler format
class EmbeddingDataset(Dataset):

    def __init__(self, datasource):
        self._datasource = datasource
        self._conn = sqlite3.connect(datasource,
                                    uri=True,
                                    check_same_thread=False)

    def __len__(self):
        with closing(self._conn.cursor()) as cur:
            cur.execute('SELECT count(*) FROM smiles')
            rec = cur.fetchone()
        return rec[0]

    def __getitem__(self, idx):
        with closing(self._conn.cursor()) as cur:
            cur.execute(
                '''
                SELECT smiles, embedding, embedding_dim,
                   logp, wt, hdonors, hacceptors, rbonds, qed
                FROM smiles
                WHERE id = ?
                ''', [idx])
            rec = cur.fetchone()
        emb = np.reshape(pickle.loads(rec[1]), pickle.loads(rec[2])).squeeze()
        emb = torch.from_numpy(emb).to(device)
        emb = torch.flatten(emb).float()
        sample = (emb, torch.tensor(rec[8]).to(device))
        # sample = (emb,
        #           torch.tensor([rec[3], rec[4], rec[5], rec[6], rec[7], rec[8]]))
        return sample


class Generator(torch.nn.Module):

    def __init__(self,
                 emb_dim=512*512,
                 latent_dim=512):
        super().__init__()
        self.layer1 = torch.nn.Sequential(torch.nn.Linear(in_features=latent_dim + 1, out_features=256),
                                          torch.nn.LeakyReLU())
        self.layer2 = torch.nn.Sequential(torch.nn.Linear(in_features=256, out_features=512),
                                          torch.nn.LeakyReLU())
        self.output = torch.nn.Sequential(torch.nn.Linear(in_features=512, out_features=emb_dim),
                                          torch.nn.Tanh())

    def forward(self, z, y):
        x = torch.cat([z, torch.reshape(y, (y.shape[0], 1))], dim=1)

        x = self.layer1(x)
        x = self.layer2(x)
        x = self.output(x)
        return x


class Discriminator(torch.nn.Module):

    def __init__(self, emb_dim=512*512):
        super().__init__()
        self.layer1 = torch.nn.Sequential(torch.nn.Linear(in_features=emb_dim + 1, out_features=512),
                                          torch.nn.LeakyReLU())
        self.layer2 = torch.nn.Sequential(torch.nn.Linear(in_features=512, out_features=256),
                                          torch.nn.LeakyReLU())
        self.output = torch.nn.Sequential(torch.nn.Linear(in_features=256, out_features=1),
                                          torch.nn.Sigmoid())

    def forward(self, x, y):
        x = torch.cat([x, torch.reshape(y, (y.shape[0], 1))], dim=1)
        x = self.layer1(x)
        x = self.layer2(x)
        x = self.output(x)
        return x


class CGAN(pl.LightningModule):

    def __init__(self, latent_dim=512):
        super().__init__()
        self.latent_dim = latent_dim
        self.generator = Generator(latent_dim=latent_dim)
        self.discriminator = Discriminator()

    def forward(self, z, y):
        return self.generator(z, y)

    def generator_step(self, x):
        z = torch.randn(x.shape[0], self.latent_dim, device=device)
        y = torch.randn((x.shape[0]), device=device)
        generated_embs = self(z, y)

        d_output = torch.squeeze(self.discriminator(generated_embs, y))
        g_loss = torch.nn.BCELoss()(d_output,
                                    torch.ones(x.shape[0], device=device))

        return g_loss

    def discriminator_step(self, x, y):
        # Real emb
        # x = torch.flatten(x)
        d_output = torch.squeeze(self.discriminator(x, y))
        loss_real = torch.nn.BCELoss()(d_output,
                                       torch.ones(x.shape[0], device=device))

        # Fake emb
        z = torch.randn(x.shape[0], self.latent_dim, device=device)
        # pdb.set_trace()
        y = torch.randn((x.shape[0]), device=device)

        generated_embs = self(z, y)
        d_output = torch.squeeze(self.discriminator(generated_embs, y))
        loss_fake = torch.nn.BCELoss()(d_output,
                                       torch.zeros(x.shape[0], device=device))

        return loss_real + loss_fake

    def training_step(self, batch, batch_idx, optimizer_idx):
        X, y = batch

        # train generator
        if optimizer_idx == 0:
            loss = self.generator_step(X)

        # train discriminator
        if optimizer_idx == 1:
            loss = self.discriminator_step(X, y)

        return loss

    def configure_optimizers(self):
        g_optimizer = torch.optim.Adam(self.generator.parameters(), lr=0.0002)
        d_optimizer = torch.optim.Adam(self.discriminator.parameters(), lr=0.0002)
        return [g_optimizer, d_optimizer], []


if __name__ == "__main__":
    data = EmbeddingDataset('/clara/virtual_screening/dataset.db')
    dataloader = DataLoader(data, batch_size=32, shuffle=True, num_workers=0)

    model = CGAN(latent_dim=255)

    trainer = pl.Trainer(max_epochs=50, gpus=1 \
        if torch.cuda.is_available() else 0, progress_bar_refresh_rate=50)
    trainer.fit(model, dataloader)
