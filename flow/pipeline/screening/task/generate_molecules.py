import logging
import concurrent.futures

import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

from flow.pipeline.screening.data import ScreeningData
from flow.utils.megamolbart import sample


log = logging.getLogger(__name__)

LIPINSKI_RULE = {'MolLogP': 5,
                 'MolWt': 500,
                 'NumHDonors': 10,
                 'NumHAcceptors': 5}

# Code to disable rdkit errors and warning
logger = rkl.logger()
logger.setLevel(rkl.ERROR)
rkrb.DisableLog('rdApp.error')


class GenerateMolecule():

    input_spec = ['smiles', 'generation']
    output_spec = ['curr_generation']

    def _score_smiles(self, gen_smiles, min_score):
        '''
        Score generated smiles and returns the ones with score above min_score.
        Returns:
            list: [{smiles: smiles,
                    embedding: embedding,
                    score: lipinski_score}]
        '''
        log.debug(
            f'Applying Lipinski score to filter {len(gen_smiles)} smiles')

        scored_list = []
        for smiles_details in gen_smiles:
            try:
                score = 0
                smi = smiles_details['smiles']
                m = Chem.MolFromSmiles(smi)

                score += 1 if Descriptors.MolLogP(m) >= LIPINSKI_RULE['MolLogP'] else 0
                score += 1 if Descriptors.MolWt(m) >= LIPINSKI_RULE['MolWt'] else 0
                score += 1 if Lipinski.NumHDonors(m) >= LIPINSKI_RULE['NumHDonors'] else 0
                score += 1 if Lipinski.NumHAcceptors(m) >= LIPINSKI_RULE['NumHAcceptors'] else 0
                log.debug(f'Lipinski score {smi}: {score}')
                smiles_details['score'] = score
                scored_list.append(smiles_details)
            except Exception as e:
                pass

        # Filter smiles which score for at least two of the Lipinski rules
        scored_list = list(
            filter(lambda a: a['score'] >= min_score, scored_list))
        log.debug(f'After filtering {len(scored_list)}, min_score: {min_score}')
        return scored_list

    def _generate_smiles(self,
                         smile_id,
                         smi,
                         num_sample,
                         service_port,
                         data,
                         min_lipinsky_score):
        generatedSmiles = []
        bal_sample = num_sample

        if type(smi) == str:
            smi = [smi]

        for smiles in smi:
            while True:
                batch_num_sample = min(bal_sample, 10)
                try:
                    generatedSmiles.extend(sample(smiles,
                                                  radius=1,
                                                  num_sample=batch_num_sample,
                                                  service_port=service_port))
                except Exception as e:
                    log.warning(f'{smiles} generated an exception: {e}')
                    pass
                bal_sample -= batch_num_sample
                if bal_sample <= 0:
                    break

        generatedSmiles = self._score_smiles(
            generatedSmiles, min_lipinsky_score)
        log.debug(f'Generated {len(generatedSmiles)}')

        # The service returns input smiles, therefore len should be at-least two
        assert len(generatedSmiles) >= 2, 'Generated smiles should be at-least 2'
        data.add_generated_smiles(smile_id, generatedSmiles)
        return smile_id, generatedSmiles

    def __call__(self, **kwargs):
        self.min_jitter_radius = 1
        generation = kwargs['generation']
        num_sample = kwargs['num_sample']
        service_port = kwargs['service_port']
        span = kwargs['span']
        min_lipinsky_score = kwargs.get('min_lipinsky_score', 0)

        log.info(f'Generating {num_sample} samples around gen {generation} molecules...')
        data = ScreeningData()
        all_smiles = data.fetch_smiles_by_generation(generation)

        with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
            for curr_span in range(1, span):
                log.info(f'Generating span {curr_span} {all_smiles}')
                futures = {executor.submit(
                    self._generate_smiles,
                    smile_id, smiles, num_sample,
                    service_port, data, min_lipinsky_score):
                    (smile_id, smiles) for (smile_id, smiles) in all_smiles.items()}

                all_smiles = {}
                for future in concurrent.futures.as_completed(futures):
                    smiles = futures[future]
                    try:
                        smi_id, generatedSmiles = future.result()

                        g_smis = [gs['smiles'] for gs in generatedSmiles]
                        all_smiles[smi_id] = g_smis

                    except Exception as exc:
                        log.warning(f'{smiles} generated an exception: {exc}')

                if len(all_smiles) == 0:
                    log.warning(f'Generation during span{curr_span} is empty')
                    break
        return generation
