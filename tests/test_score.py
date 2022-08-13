import pytest
import logging

from flow.pipeline.screening.task import GenerateMolecule

log = logging.getLogger(__name__)


def test_create_task():
    generator = GenerateMolecule()

    mol_str = ['O=C(N1[C@@H](CCC1)C(=O)[O-])[C@@H](N[C@H](C(=O)[O-])CCc1ccccc1)CCCC[NH3+]',
               'CCCCC[C@H](N[C@@H](CCc1ccccc1)C(=O)[N+](=O)[O-])C(=O)N1CCCC1',
               'CCCCC[C@H](N[C@@H](CCc1ccccc1)[N+](=O)[O-])C(=O)N1CCC[C@H]1C(N)=O',
               'CCCCC[C@H](N[C@@H](CCc1ccccc1)C(=O)N1CCCC1)C(=O)N1CCC[C@@H]1C(=O)[N+](=O)[O-]',
               'CCCC[C@H](N[C@@H](CCCC)C(=O)N1CCN(C)C[C@H]1C[N+][N+](=O)[O-])C(=O)N1CCN(C)C[C@H]1N',
               'CCCC[C@H](N[C@@H](CCc1ccccc1)C(=O)N1C[C@@H]2CCC[C@H]21)C(=O)NCCc1ccccc1',
               'CCCC[C@H](N[C@@H](CCc1ccccc1)C(=O)OC)C(=O)N1C[C@@H]2CC[C@H]1[C@@H]([N+](=O)[O-])C2',
               'CCCCC1CC1',
               'CCCC[C@H](N[C@H]1CCc2ccccc2C1)C(=O)N1CCC[C@@H]1C',
               'CCCC[C@H](N[C@@H](CCc1ccccc1)N=[N+]([O-])[O-])C(=O)N1CCC[C@H]1C[N+](=O)[O-]',
               'CCCCC[C@H](N[C@@H](CCCc1ccccc1)C(=O)N1CCCC[C@]1(C)Cc1ccccc1)C(=O)O',
               'CCCCC[C@H](N[C@@H](CCc1ccccc1)C(=O)[N+](=O)[O-])C(=O)N1CCC[C@@H]1C(=O)[N+](=O)[O-]',
               '[CH]CCCC[C@H](N[C@@H](CCc1ccccc1)C(N)=O)C(=O)N1[C@@H]2CCC[C@H]2C[C@H]1[N+](=O)[O-]',
               'CCCC[C@H]1CCCN1C(=O)[C@H](CCc1ccccc1)N[C@@H](CCc1ccccc1)C(=O)O',
               'CCCC[C@H](N=C(O)[C@@H](N)[N+](=O)[O-])C(=O)N1CCC[C@H]1C[N+]=O',
               'CCCCC[C@H](N[C@@H](C/C=C/[N+](=O)C=CO)C(=O)OCC)C(=O)N1CCC[C@@H](CO)[C@@H]1C[N+](=O)[O-]',
               'O=C(NCCN1C[C@@H]2CCC[C@H]21)[C@H](CCc1ccccc1)N1C[C@@H]2CCC[C@H]21',
               'CCCC[C@H](N[C@@H](CCCC)C(=O)N1C[C@@H]2CCC[C@@H]2[C@H]1CCc1ccccc1)C(=O)NCCc1ccccc1',
               'CCCCNC(=O)[C@H](CCCC)N[C@@H](CCc1ccccc1)C(=O)N1CCc2ccccc2C1',
               'CCCCNC1CN(C(=O)[C@H](CCCc2ccccc2)N[C@@H](CCC)C(=O)NCCc2ccccc2)CCN[C@@H]1C',
               'CCCC[C@H](N[C@@H](CCc1ccccc1)C(=O)NC(=O)CC)C(=O)NCCc1ccccc1',
               'CC[C@H](N[C@@H](CN1CC[C@@H]2OCC[C@H]21)C(=O)NC)C(=O)N[C@@H](CCc1ccccc1)C(=O)NCCc1ccccc1',
               'CCCC[C@H](N[C@@H](CCC)C(=O)NCCc1ccccc1)C(=O)NCCc1ccccc1',
               'CCCC[C@H](N[C@H]1CN[C@@H](CCc2ccccc2)C1)C(=O)NCCc1ccccc1',
                ]

    min_score = 2
    filtered_smiles = generator._score_smiles(mol_str, min_score)
    scores = dict(filter(lambda a: a[1] < min_score, filtered_smiles.items()))

    assert len(scores) == 0, f'Score function should return smiles with score greater than {min_score}'
