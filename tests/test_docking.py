import io
import pytest
import logging
import random
import sqlite3
import asyncio
from contextlib import closing

from flow.dsl import task, pipeline, loop, condition
from flow.config import Configuration
from flow.utils.stringutils import substitute
from flow.pipeline.screening.task import (PrepareProtein,
                                         ScoreMolecule)


log = logging.getLogger(__name__)


class MoleculeQueue():

    def __init__(self, **kwargs):
        self.chembl_db = kwargs['chembl_db']
        self.conn = sqlite3.connect(self.chembl_db, uri=True)

    def __call__(self, **kwargs):
        generation = kwargs['_loop_count']
        row_id = random.randint(1, 1941411)
        log.debug(f'Fetching row {row_id} - {generation} from the database')
        self._create_table()
        cur = self.conn.execute('''
                                SELECT canonical_smiles
                                FROM compound_structures
                                WHERE rowid = ?
                                ''', [row_id])
        recs = cur.fetchone()
        smiles = recs[0]
        log.info(f'SMILES from ChEMBL {smiles}')
        with closing(self.wf_db.cursor()) as cur:
            cur.execute('''
                        INSERT INTO smiles(smiles, generation)
                        VALUES(?, ?)''', [smiles, generation])
            generation = cur.lastrowid
            log.info(f'ID {generation}')
            cur.execute('''
                        INSERT INTO generated_smiles(input_id, smiles)
                        VALUES(?, ?)''', [generation, smiles])
            self.wf_db.commit()
        return generation

    def _create_table(self):
        # TODO: Revisit. This is a hack to get the db_file to be substituted
        db_file = Configuration().get('platform_db')
        db_file = substitute(db_file, Configuration().config)
        log.info(f'Creating database tables requred for the task {db_file} {Configuration().config}')

        self.wf_db = sqlite3.connect(db_file, uri=True)
        ddl_stmt = f'''
                    CREATE TABLE IF NOT EXISTS smiles (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        smiles TEXT NULL,
                        generation INTEGER NOT NULL,
                        UNIQUE(smiles)
                    );

                    CREATE TABLE IF NOT EXISTS generated_smiles (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        input_id INTEGER NOT NULL,
                        smiles TEXT NULL,
                        score REAL,
                        score_model INTEGER,
                        cid TEXT NULL,
                        UNIQUE(smiles)
                    );

                    CREATE TABLE IF NOT EXISTS conformer_score (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        gs_id INTEGER NOT NULL,
                        cid TEXT NULL,
                        score REAL,
                        score_model INTEGER,
                        UNIQUE(gs_id, cid)
                    );
                    '''
        self.wf_db.executescript(ddl_stmt)


def next_generation(**kwargs):
    kwargs['generation'] = kwargs['generation'] + 1
    return kwargs


def terminal_condition(**kwargs):
    if kwargs['score'] is None:
        return True
    return kwargs['score'] > -13.8


@pytest.mark.skip(reason="This is used for verification and benchmarking only")
def test_vs_pipeline():
    '''
    THIS IS NOT A TEST. IT IS DEVELOPED FOR BENCHMARKING PURPOSES.
    '''

    config = '''
    platform:
        type: containers
        ws_root: /tmp/shared
        db: /tmp/shared/{{pipeline_id}}/common.sqlite3
        db_link: /pipeline/common.sqlite3
    '''
    data = io.StringIO()
    data.write(config)
    config = Configuration()
    config.load_config(data)
    log.info(f'Starting pipeline with config: {config.config}')

    virtual_screening =\
        pipeline(name='virtual_screening')(
            task(name='prepare_protein',
                 inputs=['protein_file', 'pp_repair_type', 'score_config'],
                 outputs=['pdbqt_file'])(PrepareProtein()),
            loop(name='Check Pose Score')(
                task(name='MoleculeQueue',
                     inputs=['_loop_count'],
                     outputs=['generation'])(
                         MoleculeQueue(chembe_db='/clara/testData/ChEMBL/chembl_27.db')),

                task(name='ScoreMolecule',
                     inputs=['generation', 'pdbqt_file'],
                     outputs=['score'])(ScoreMolecule()),
                condition(name='Check Pose Score',
                          inputs='score')(terminal_condition),
            )
        )

    pipeline_input = {'protein_file': 'file://./tests/data/6y2g/6y2g_clean.pdb',
                      'score_config': 'file://./tests/data/6y2g/config',
                      'pp_repair_type': 'checkhydrogens',
                      'generation': 0}
    asyncio.run(virtual_screening(**pipeline_input))
