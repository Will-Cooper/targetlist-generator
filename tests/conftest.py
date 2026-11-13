import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def target_csv(tmp_path):
    """25 targets spread over the sky, mixed coordinate formats and numeric/string spectral types."""
    rng = np.random.default_rng(1)
    ra = rng.uniform(0, 360, 25)
    dec = np.degrees(np.arcsin(rng.uniform(-0.9, 0.9, 25)))
    df = pd.DataFrame({'Name': [f'J{i:04d} test' for i in range(25)], 'RA': ra, 'Dec': dec,
                       'SpT': rng.choice([65.0, 72.5, 'T6', 'L1'], 25), 'TMASSK': rng.uniform(10, 15, 25)})
    df[['RA', 'Dec']] = df[['RA', 'Dec']].astype(object)
    df.loc[0, ['RA', 'Dec']] = ['23:59:59.9999', '-00:30:00']  # wrap and sign edge cases
    path = tmp_path / 'targets.csv'
    df.to_csv(path, index=False)
    return path
