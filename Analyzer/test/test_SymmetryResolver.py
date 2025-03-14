import pytest
from SymmetryResolverClass import SymmetryResolver


@pytest.fixture
def symResolver():
    return SymmetryResolver(
        nNeighbors=3,
        nNextNeighbors=6,
        runsPath="",
        matchPattern="",
        sublattices=2,
        subbands=1,
    )


class TestSymmetryResolverSuite:

    def setup_method(self, symResolver):
        self.symResolver = symResolver
        print("Initialized sym resolver")

    def test_CalculateSymmetryGamma(self):
        pass
