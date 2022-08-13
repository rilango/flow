import logging
import sqlite3
import pickle

from contextlib import closing

from flow.config import Configuration
from flow.utils.singleton import Singleton
from flow.utils.stringutils import substitute

# Code to disable rdkit errors and warning
import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl

import warnings
warnings.filterwarnings('ignore')

log = rkl.logger()
log.setLevel(rkl.ERROR)
rkrb.DisableLog('rdApp.error')

log = logging.getLogger(__name__)


class ScreeningData(metaclass=Singleton):

    def __init__(self, datasource=None):
        self._conn = None
        self.datasource = datasource
        if datasource is not None:
            self.initialize_db()

    def initialize_db(self):

        # TODO: Revisit. This is a hack to get the db_file to be substituted
        if self.datasource is None:
            db_file = Configuration().get('platform_db')
            db_file = substitute(db_file, Configuration().config)
            self.datasource = db_file
        log.info(f'Creating database at {self.datasource}...')

        self._conn = sqlite3.connect(self.datasource,
                                     uri=True,
                                     check_same_thread=False)
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
                        embedding TEXT NOT NULL,
                        embedding_dim TEXT NOT NULL,
                        score REAL,
                        score_model INTEGER,
                        lipinski_score REAL,
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

                    CREATE TABLE IF NOT EXISTS stats (
                        id INTEGER PRIMARY KEY AUTOINCREMENT,
                        stat_type TEXT NOT NULL,
                        stat_key TEXT NOT NULL,
                        value TEXT
                    );
                    '''
        self._conn.executescript(ddl_stmt)

    def add_stat(self, stat_type, stat_key, value):
        '''
        Add a stat to the database.
        '''
        with closing(sqlite3.connect(self.datasource, uri=True)) as con, con, \
                closing(con.cursor()) as cur:
            cur.execute(
                'INSERT INTO stats(stat_type, stat_key, value) VALUES(?, ?, ?)',
                [stat_type, stat_key, value])

    def add_smiles(self, smiles, generation):
        '''
        Adds new smiles to a generation. If the smiles already exists, then
        an exception is raised.
        '''
        all_smiles = []
        if isinstance(smiles, str):
            all_smiles = [smiles]
        else:
            all_smiles = smiles

        added_smiles = {}
        with closing(sqlite3.connect(self.datasource, uri=True)) as con, con, \
                closing(con.cursor()) as cur:
            for new_smiles in all_smiles:
                try:
                    cur.execute(
                        'INSERT INTO smiles(smiles, generation) VALUES(?, ?)',
                        [new_smiles, generation])
                    added_smiles[cur.lastrowid] = new_smiles
                except sqlite3.IntegrityError:
                    log.warning(f'{new_smiles} already exists in the database')
            # self._conn.commit()

        return added_smiles

    def add_generated_smiles(self, smile_id, gsmiles_details):
        '''
        Add a list of generated smiles.
        '''
        with closing(sqlite3.connect(self.datasource, uri=True)) as con, con, \
                closing(con.cursor()) as cur:
            for gsmiles_details in gsmiles_details:
                generated_smile = gsmiles_details['smiles']
                score = gsmiles_details['score']
                embedding = gsmiles_details['embedding']
                embedding_dim = pickle.dumps(list(embedding.dim))
                embedding = pickle.dumps(list(embedding.embedding))

                try:
                    # The first smile is the original one. This must be
                    # accounted for while using this information.
                    cur.execute('''
                                INSERT INTO generated_smiles(
                                    input_id, smiles, lipinski_score, embedding, embedding_dim)
                                VALUES(?, ?, ?, ?, ?)
                                ''', [smile_id, generated_smile, score, embedding, embedding_dim])
                except sqlite3.IntegrityError:
                    log.warning(f'Generated smiles {generated_smile} already exists in the database')


    def add_ensamble(self, smi, gsmi_details, generation):
        '''
        Adds new smiles to a generation. If the smiles already exists, then
        an exception is raised.
        '''
        added_smiles = {}
        with closing(sqlite3.connect(self.datasource, uri=True)) as con, con, \
                closing(con.cursor()) as cur:
            try:
                cur.execute(
                    'INSERT INTO smiles(smiles, generation) VALUES(?, ?)',
                    [smi, generation])
                smiles_id = cur.lastrowid

                for gsmi in gsmi_details:
                    generated_smile = gsmi['smiles']
                    score = gsmi['score']
                    embedding = gsmi['embedding']
                    embedding_dim = pickle.dumps(list(embedding.dim))
                    embedding = pickle.dumps(list(embedding.embedding))
                    try:
                        cur.execute('''
                                    INSERT INTO generated_smiles(
                                        input_id, smiles, lipinski_score, embedding, embedding_dim)
                                    VALUES(?, ?, ?, ?, ?)
                                    ''', [smiles_id, generated_smile, score, embedding, embedding_dim])
                    except sqlite3.IntegrityError:
                        log.warning(f'Generated smiles {generated_smile} already exists in the database')

            except sqlite3.IntegrityError:
                log.warning(f'{smi} already exists in the database')


    def add_conformer_score(self, gsmiles_id, cid, min_score, score_model):
        '''
        Add a conformer score to the database.
        '''
        with closing(sqlite3.connect(self.datasource, uri=True)) as con, con, \
                closing(con.cursor()) as cur:
            log.debug(f'Inserting score for {gsmiles_id}_{cid} {min_score}')
            cur.execute('''
                INSERT INTO conformer_score(gs_id, cid, score, score_model)
                VALUES(?, ?, ?, ?)''',
                [gsmiles_id, cid, min_score, score_model])

            log.debug(f'Updating score for {gsmiles_id}_{cid} {min_score}')
            cur.execute('''
                UPDATE generated_smiles
                SET score = ?,
                    score_model = ?,
                    cid = ?
                WHERE id = ?
                and (score is NULL or score > ?)''',
                [min_score, score_model, cid, gsmiles_id, min_score])

    def fetch_smiles_by_generation(self, generation):
        '''
        Fetch smiles from the database by generation. These smiles are either
        the original input or the generated ones selected from the previous
        generation.
        '''
        assert self._conn is not None, 'Database connection is not initialized'

        with closing(self._conn.cursor()) as cur:
            cur.execute('''
                        SELECT s.id, s.smiles
                        FROM smiles s
                        WHERE s.generation = ?
                        ''', [generation])
            recs = cur.fetchall()

        all_smiles = {}
        for rec in recs:
            all_smiles[rec[0]] = rec[1]

        return all_smiles

    def fetch_gsmiles(self, include_unscored=False):
        '''
        Fetch all generated smiles from the database.
        '''
        assert self._conn is not None, 'Database connection is not initialized'

        score_condition = 'AND gs.score is not NULL'
        if include_unscored is False:
            score_condition = ''

        with closing(self._conn.cursor()) as cur:
            cur.execute(f'''
                        SELECT gs.id as gsmiles_id,
                               gs.smiles as smiles,
                               gs.score as score,
                               gs.embedding,
                               gs.embedding_dim as dim
                        FROM generated_smiles gs, smiles s
                        WHERE gs.input_id = s.id
                            AND gs.score is not NULL
                            {score_condition}
                        ''')
            recs = cur.fetchall()

        all_smiles = {}
        for rec in recs:
            all_smiles[rec[0]] = {
                'smiles': rec[1],
                'score': rec[2],
                'emb': pickle.loads(rec[3]),
                'emb_dim': pickle.loads(rec[4])
            }
        return all_smiles

    def fetch_gsmiles_by_generation(self,
                                    generation,
                                    include_unscored=False,
                                    top_k=None):
        '''
        Fetch generated smiles from the database by generation. These smiles are
        generated ones created using previous generation.
        '''
        assert self._conn is not None, 'Database connection is not initialized'

        score_condition = 'AND gs.score is not NULL'
        if include_unscored:
            score_condition = ''

        with closing(self._conn.cursor()) as cur:
            cur.execute(f'''
                        SELECT gs.id as gsmiles_id,
                               gs.smiles as smiles,
                               gs.score as score,
                               gs.embedding,
                               gs.embedding_dim as dim
                        FROM generated_smiles gs, smiles s
                        WHERE gs.input_id = s.id
                            AND s.generation = ?
                            {score_condition}
                        ORDER BY gs.score
                        ''', [generation])
            recs = cur.fetchall()
        if top_k is not None:
            recs = recs[:top_k]

        all_smiles = {}
        for rec in recs:
            all_smiles[rec[0]] = {
                'smiles': rec[1],
                'score': rec[2],
                'emb': pickle.loads(rec[3]),
                'emb_dim': pickle.loads(rec[4])
            }

        return all_smiles

    def fetch_smiles_by_generation_range(self,
                                         start_generation,
                                         end_generation,
                                         include_unscored=False):
        '''
        Fetch generated smiles from the database by generation. These smiles are
        generated ones created using previous generation.
        '''
        assert self._conn is not None, 'Database connection is not initialized'
        score_condition = 'AND gs.score is not NULL'
        if include_unscored is False:
            score_condition = ''

        with closing(self._conn.cursor()) as cur:
            cur.execute(f'''
                        SELECT gs.id as gsmiles_id,
                               gs.smiles as smiles,
                               gs.score as score,
                               gs.embedding,
                               gs.embedding_dim as dim
                        FROM generated_smiles gs, smiles s
                        WHERE gs.input_id = s.id
                            AND s.generation between ? and ?
                            {score_condition}
                        ''', [start_generation, end_generation])
            recs = cur.fetchall()

        all_smiles = {}
        for rec in recs:
            all_smiles[rec[0]] = {
                'smiles': rec[1],
                'score': rec[2],
                'emb': pickle.loads(rec[3]),
                'emb_dim': pickle.loads(rec[4])
            }

        return all_smiles

    def fetch_min_score(self, generation):
        '''
        Fetch the minimum score for a generation.
        '''
        with closing(self._conn.cursor()) as cur:
            cur.execute('''
                        SELECT min(gs.score)
                        FROM smiles s, generated_smiles gs
                        WHERE gs.input_id = s.id
                            AND s.generation = ?
                            AND gs.score is not null
                        ''', [generation])
            recs = cur.fetchone()

        return recs[0]

    def promote_smiles(self, generation, num_smiles):
        '''
        Promote smiles from the cuurent generation to next generation.
        '''
        with closing(sqlite3.connect(self.datasource, uri=True)) as con, con, \
                closing(con.cursor()) as cur:
            cur.execute('''
                        SELECT count(*)
                        FROM smiles s, generated_smiles gs
                        WHERE gs.input_id = s.id
                            AND s.generation = ?
                        ''', [generation])
            cnt_mol_scored = cur.fetchone()

            assert cnt_mol_scored[0] > 1, 'Not enough molecules scored'

            log.info(f'Screening {num_smiles} molecules for next generation')
            cur.execute('''
                        SELECT gs.smiles
                        FROM smiles s, generated_smiles gs
                        WHERE gs.input_id = s.id
                            AND s.generation = ?
                            AND gs.score is not null
                            AND s.smiles <> gs.smiles
                        ORDER by gs.score ASC
                        LIMIT ?
                        ''', [generation, num_smiles])
            selected_mols = cur.fetchall()

            for smiles in selected_mols:
                log.info(f'Adding {smiles} to next generation')
                try:
                    cur.execute('INSERT INTO smiles(smiles, generation) VALUES(?, ?)',
                                [smiles[0], generation + 1])
                except sqlite3.IntegrityError:
                    log.warning(f'{smiles} already in a previous generation')

        return generation