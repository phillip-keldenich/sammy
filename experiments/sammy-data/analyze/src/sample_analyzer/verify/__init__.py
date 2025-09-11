import itertools
import typing

from ._verify import (have_equal_coverage, count_interactions)


def _minimize_to_concrete_features(
        sample: typing.List[typing.Dict[str, bool]], concrete_features: typing.Set[str]
) -> typing.List[typing.Dict[str, bool]]:
    return [
        {k: v for k, v in conf.items() if k in concrete_features} for conf in sample
    ]


def _unique_interaction(f1, f1_val, f2, f2_val):
    if f1 < f2:
        return (f1, f1_val, f2, f2_val)
    else:
        return (f2, f2_val, f1, f1_val)

def retrieve_interactions(
        sample: typing.List[typing.Dict[str, bool]],
        concrete_features: list):
    if sample is None:
        return None

    concrete_features = set(concrete_features)
    sample = _minimize_to_concrete_features(sample, concrete_features)
    interactions = {
        _unique_interaction(f0[0], f0[1], f1[0], f1[1])
        for conf in sample
        for f0, f1 in itertools.combinations(conf.items(), 2)
    }
    return interactions

def count_interactions_py(
        sample: typing.List[typing.Dict[str, bool]],
        concrete_features: list,
) -> int:
    if sample is None:
        return None

    return len(retrieve_interactions(sample, concrete_features))