import sys
import logging
import itertools
import argparse
import pdb

from cuml import Ridge
import cuml
import cupy as cp

from rdkit import Chem
from rdkit.Chem import QED

from generativesampler_pb2 import EmbeddingList

from flow.pipeline.screening.data import ScreeningData
from flow.pipeline.screening.pose_generate import score_molecule, generate_conformers
from flow.utils.megamolbart import interpolate, smiles_to_embedding, embedding_to_smiles, sample as mm_sample

logging.basicConfig(format='%(asctime)s %(name)s [%(levelname)s]: %(message)s',
                    level=logging.DEBUG)
log = logging.getLogger(__name__)


def argparser(argv):
    parser = argparse.ArgumentParser(description='Extraplote embedding')

    parser.add_argument('-g', '--generation',
                        dest='generation',
                        type=int,
                        required=True,
                        help='Generation of the ligand library')

    parser.add_argument('-n', '--num_mol',
                        dest='num_mol',
                        type=int,
                        default=4,
                        required=False,
                        help='Number of molecules to be generated.')

    parser.add_argument('--db',
                        dest='db_url',
                        type=str,
                        default=None,
                        required=False,
                        help='Database connection string.')

    args = parser.parse_args(argv)
    return args


def add_jitter(embedding, radius, cnt):
    distorteds = []
    for _ in range(cnt):
        noise = cp.random.normal(0, radius, embedding.shape)
        distorted = noise + embedding
        distorteds.append(cp.asnumpy(distorted))

    return distorteds


class Regression():

    input_spec = ['generation']
    output_spec = ['curr_generation']

    def optimize(self, **kwargs):

        generation = kwargs['generation']
        num_mol = kwargs['num_mol']
        db_url = kwargs['db_url']

        data = ScreeningData(datasource=db_url)
        # all_smiles = data.fetch_gsmiles()
        # all_smiles = data.fetch_smiles_by_generation_range(generation - 3,
        #                                                    generation + 1)
        all_smiles = data.fetch_gsmiles_by_generation(generation)

        assert len(all_smiles) > 0, 'No molecules to optimize'
        log.debug(f'Input data containes {len(all_smiles)} geneated molecules')

        orig_embs = []
        embs = []
        dims = []
        scores = []
        for details in all_smiles.values():
            scores.append(details['score'])
            emb = smiles_to_embedding(details['smiles'])
            orig_embs.append(emb)
            embs.append(cp.reshape(cp.array(list(emb.embedding)), list(emb.dim)))
            dims.append(list(emb.dim))

        scores = cp.array(scores)
        embs = cp.array(embs)

        clf = Ridge(alpha=1.0)
        clf.fit(embs, scores);
        # Predict scores
        pred_scores = clf.predict(embs)
        r2 = cuml.metrics.regression.r2_score(scores, pred_scores)
        err = cuml.metrics.regression.mean_squared_error(scores, pred_scores)
        data.add_stat('r2', generation, str(r2))
        data.add_stat('mse', generation, str(err))

        # error = scores - pred_scores
        min_idx = int(cp.argmin(pred_scores))
        pred_min_emb = cp.array(embs[min_idx])

        #TODO: Verify coef_ can be used like this
        endpoint = cp.random.normal(clf.coef_, 1, pred_min_emb.shape)
        samples = cp.linspace(cp.zeros(pred_min_emb.shape),
                              endpoint,
                              endpoint=False,
                              num=(num_mol + 1))

        orig_emb = orig_embs[min_idx]
        dim = orig_emb.dim
        mask = orig_emb.pad_mask

        mols = []
        gsmiles = []
        for jittered_emb in samples[1:]:
            jittered_emb = jittered_emb.flatten().tolist()
            m_gsmiles = embedding_to_smiles(jittered_emb, dim, mask)
            try:
                mol = Chem.MolFromSmiles(m_gsmiles.generatedSmiles[0])
                if mol is not None:
                    mols.append(mol)
                    gsmiles.append(m_gsmiles.generatedSmiles[0])
            except Exception as ex:
                pass

        gsmiles = set(gsmiles)
        # Add the predicted smiles to next generation
        data.add_smiles(gsmiles, generation + 1)
        log.info(f'Generated {gsmiles} molecules')
        return gsmiles


class SimpleQEDScreen():

    input_spec = ['generation']
    output_spec = ['curr_generation']

    def optimize(self, **kwargs):
        generation = kwargs['generation']
        db_url = kwargs['db_url']
        num_mol = kwargs['num_mol']

        data = ScreeningData(datasource=db_url)
        all_smiles = data.fetch_gsmiles_by_generation(generation)

        assert len(all_smiles) > 0, 'No molecules to optimize'
        log.debug(f'Input data containes {len(all_smiles)} geneated molecules')

        scores = {}
        for details in all_smiles.values():
            smi = details['smiles']
            docking_score = details['score']
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    continue

                qed = QED.qed(mol)
                scores[smi] = (1 - qed) * docking_score
            except Exception as ex:
                pass

        # Sort by score
        scores = {k: v for k, v in sorted(scores.items(), key=lambda item: item[1])}
        selected_smiles = list(scores.keys())[: num_mol]
        for smi in selected_smiles:
            log.info(f'Selected molecules {smi}, score: {scores[smi]}')
            data.add_smiles(smi, generation + 1)

        return selected_smiles


class RadialBasisScreen():

    def multiquadric_rbf(self, radius, epsilon):
            return cp.sqrt(1 + (epsilon * radius)**2)

    def optimize(self, **kwargs):
        generation = kwargs['generation']
        db_url = kwargs['db_url']
        num_mol = kwargs['num_mol']
        service_port = kwargs.get('service_port', 'localhost:50051')
        receptor_file = kwargs.get('receptor_file', '/pipeline/inputs/rec.pdbqt')
        score_config = kwargs.get('score_config', '/pipeline/inputs/config')
        epsilon = kwargs.get('epsilon', 0.01)

        data = ScreeningData(datasource=db_url)
        all_smiles = data.fetch_gsmiles_by_generation(generation, top_k=100)

        assert len(all_smiles) > 0, 'No molecules to optimize'
        log.debug(f'Input data containes {len(all_smiles)} geneated molecules')

        x0_smis = []
        x0_dims = []
        x0_embs = []
        y0_scrs = []

        x1_smis = []
        x1_dims = []
        x1_embs = []
        y1_scrs = []

        seed_size = int(len(all_smiles) * 0.8)
        print(f'Seed size: {seed_size}')
        for _, smi_details in list(all_smiles.items())[:seed_size]:
            dim = smi_details['emb_dim']
            emb = smi_details['emb']
            emb = cp.reshape(cp.asarray(emb), dim).squeeze().flatten()

            x0_dims.append(dim)
            x0_embs.append(emb)
            x0_smis.append(smi_details['smiles'])
            y0_scrs.append(smi_details['score'])

        for _, smi_details in list(all_smiles.items())[seed_size:]:
            dim = smi_details['emb_dim']
            emb = smi_details['emb']
            emb = cp.reshape(cp.asarray(emb), dim).squeeze().flatten()

            x1_dims.append(dim)
            x1_embs.append(emb)
            x1_smis.append(smi_details['smiles'])
            y1_scrs.append(smi_details['score'])

        r = cuml.metrics.pairwise_distances(
            cp.asarray(x0_embs),
            cp.asarray(x1_embs),
            metric='euclidean').T

        # TODO: Revisit for bigger dataset
        r1 = []
        for emb1 in x1_embs:
            row = []
            for emb0 in x0_embs:
                row.append(emb0.squeeze().flatten() - emb1.squeeze().flatten())
            r1.append(cp.asarray(row))
        r1 = cp.asarray(r1)

        # Multiquadric
        A = self.multiquadric_rbf(r, epsilon)
        w = cp.matmul(cp.linalg.pinv(A).T, cp.asarray(y0_scrs))

        tmp = ((epsilon**2) / cp.sqrt(1 + (epsilon**2 * r**2)))
        tmp = cp.multiply(tmp.reshape(tmp.shape + (1,)), r1)
        grad = cp.einsum("i,ijk->jk", w, tmp)

        for i in range(len(x0_embs)):
            try:
                emb = x0_embs[i]
                dim = x0_dims[i]
                smi = x0_smis[i]
                projected_emb  = emb - (3 * grad[i])
                projected_emb = projected_emb.tolist()
                mask = list(itertools.repeat(False, dim[0]))
                result = embedding_to_smiles(projected_emb,
                                             list(dim),
                                             mask,
                                             service_port=service_port)

                if result.generatedSmiles[0] == smi:
                    continue

                generated_smiles = interpolate([smi,
                                                result.generatedSmiles[0]],
                                                service_port=service_port)

                valid_smis = []
                for g_smiles in generated_smiles:
                    op_mol = Chem.MolFromSmiles(g_smiles)
                    if op_mol:
                        embedding = EmbeddingList(embedding=projected_emb,
                                                dim=dim,
                                                pad_mask=None)
                        valid_smis.append({'embedding': embedding,
                                        'smiles': g_smiles,
                                        'score': 0})

                data.add_ensamble(smi, valid_smis, generation + 1)
            except Exception as ex:
                pass


if __name__ == '__main__':
    args = argparser(sys.argv[1:])
    db_url = args.db_url
    generation = args.generation
    num_mol = args.num_mol

    opti = RadialBasisScreen()
    opti.optimize(generation=generation,
                  num_mol=num_mol,
                  db_url=db_url)
