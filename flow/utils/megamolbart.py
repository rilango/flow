import grpc
import logging

from typing import List

from generativesampler_pb2_grpc import GenerativeSamplerStub
from generativesampler_pb2 import GenerativeSpec, GenerativeModel, EmbeddingList

log = logging.getLogger(__name__)


def sample(smi,
           radius=1,
           num_sample=10,
           padding_size=None,
           service_port='localhost:50051'):
    '''
    Sample a set of smiles from a given set of smiles
    '''
    log.debug(f'Generating samples around {smi}...')

    spec = GenerativeSpec(model=GenerativeModel.MegaMolBART,
                          smiles=smi,
                          radius=radius,
                          padding=padding_size,
                          numRequested=num_sample,
                          forceUnique=False,
                          sanitize=True)

    with grpc.insecure_channel(f'{service_port}') as channel:
        stub = GenerativeSamplerStub(channel)
        result = stub.FindSimilars(spec)

    log.debug(f'Generated SMILES {result.generatedSmiles}')

    # generatedSmiles = list(result.generatedSmiles)
    # embeddings = list(result.embeddings)
    all_gsmiles = []
    for i in range(len(result.generatedSmiles)):
        all_gsmiles.append({'embedding': result.embeddings[i],
                            'smiles': result.generatedSmiles[i]})
    return all_gsmiles


def interpolate(smis: List[str],
                radius=1,
                num_sample=10,
                service_port='localhost:50051'):
    '''
    Interpolate a set of smiles from a given set of smiles
    '''

    spec = GenerativeSpec(model=GenerativeModel.MegaMolBART,
                          smiles=smis,
                          radius=radius,
                          numRequested=num_sample,
                          forceUnique=False,
                          sanitize=True)

    with grpc.insecure_channel(f'{service_port}') as channel:
        stub = GenerativeSamplerStub(channel)
        result = stub.Interpolate(spec)

    return result.generatedSmiles


def smiles_to_embedding(smi,
                        padding_size=512,
                        service_port='localhost:50051'):
    '''
    Fetch embedding for a smiles string in MegaMolBART.
    '''
    spec = GenerativeSpec(model=GenerativeModel.MegaMolBART,
                          smiles=smi,
                          padding=padding_size)

    with grpc.insecure_channel(f'{service_port}') as channel:
        stub = GenerativeSamplerStub(channel)
        return stub.SmilesToEmbedding(spec)


def embedding_to_smiles(embedding,
                        dim,
                        mask,
                        service_port='localhost:50051'):
    '''
    Fetch embedding for a smiles string in MegaMolBART.
    '''
    spec = EmbeddingList(embedding=embedding,
                         dim=dim,
                         pad_mask=mask)

    with grpc.insecure_channel(f'{service_port}') as channel:
        stub = GenerativeSamplerStub(channel)
        return stub.EmbeddingToSmiles(spec)
