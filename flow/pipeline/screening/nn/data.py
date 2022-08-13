import sys
import logging
import argparse
import sqlite3
import pickle
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Lipinski

from flow.utils.megamolbart import smiles_to_embedding


logging.basicConfig(format='%(asctime)s %(name)s [%(levelname)s]: %(message)s',
                    level=logging.DEBUG)
log = logging.getLogger(__name__)


def argparser(argv):
    parser = argparse.ArgumentParser(description='Create training data')

    parser.add_argument('--dataset_file',
                        dest='dataset_file',
                        type=str,
                        default='/content/dataset.h5',
                        required=False,
                        help='Database connection string.')

    args = parser.parse_args(argv)
    return args

def generate_dataset(dataset_file,
                     service_port='localhost:50051',
                     ds_size=11000):

    conn_chembl = sqlite3.connect('/data/chembl.db',
                                  uri=True,
                                  check_same_thread=False)
    rec = 1
    df = pd.DataFrame(columns=['smi', 'emb', 'dim',
                               'logp', 'wt', 'hdonors',
                               'hacceptors', 'rbonds', 'qed'])
    while rec < ds_size:
        cur_chembl = conn_chembl.execute('''
                                    SELECT canonical_smiles
                                    FROM compound_structures
                                    ORDER BY RANDOM()
                                    LIMIT 1
                                    ''')
        smi = cur_chembl.fetchone()[0]

        try:
            emb_result = smiles_to_embedding(smi,
                                             padding_size=512,
                                             service_port=service_port)
        except Exception as e:
            pass

        emb = emb_result.embedding
        dim = emb_result.dim

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        emb = pickle.dumps(np.reshape(np.array(list(emb)), list(dim)))
        dim = pickle.dumps(list(dim))

        logp = Descriptors.MolLogP(mol)
        wt = Descriptors.MolWt(mol)
        hdonors = Lipinski.NumHDonors(mol)
        hacceptors = Lipinski.NumHAcceptors(mol)
        rbonds = Lipinski.NumRotatableBonds(mol)
        qed = QED.qed(mol)

        try:
            df = df.append({'smi': smi, 'emb': emb, 'dim': dim,
                       'logp': logp, 'wt': wt, 'hdonors': hdonors,
                       'hacceptors': hacceptors, 'rbonds': rbonds, 'qed': qed}, ignore_index=True)
            rec += 1
        except Exception as e:
            print(e)

        if rec % 100 == 0:
            print(f'{rec} records processed')
            df.reset_index(drop=True, inplace=True)
            df.to_hdf(dataset_file, key='data')

    df.to_hdf(dataset_file, key='data')


if __name__ == '__main__':
    args = argparser(sys.argv[1:])
    dataset_file = args.dataset_file

    generate_dataset(dataset_file, service_port='localhost:50052')
