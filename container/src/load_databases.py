"""
Download ensemble databases
"""

from pyensembl import EnsemblRelease

data = EnsemblRelease(109)

data.download()
data.index()