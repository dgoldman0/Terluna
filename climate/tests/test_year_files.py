import numpy as np

from climate.gcm import year_files as yf


class Var:
    def __init__(self, data, dimensions):
        self.data, self.dimensions = np.asarray(data, float), dimensions

    def __getitem__(self, index):
        return self.data[index]


class YearFile:
    """A stand-in for an open year file: its variables by name, each with its dimensions."""

    def __init__(self, **variables):
        self.variables = variables

    def __getitem__(self, name):
        return self.variables[name]


def year_file(with_mean=False):
    n = 6
    start = np.arange(n, dtype=float)[:, None] * 1e-3 + 1e-2      # the humidity at each interval's start, (time, lev)
    start[0] = 2e-4                                                # the restart's empty first record
    variables = dict(hus=Var(start, ('time', 'lev')), ta=Var(np.arange(n)[:, None] + 270.0, ('time', 'lev')),
                     lev=Var([0.5], ('lev',)))
    if with_mean:
        variables[yf.HUS_MEAN] = Var(start + 5e-4, ('time', 'lev'))
    return YearFile(**variables), start


def test_humidity_means_average_each_record_with_the_next():
    d, start = year_file()
    assert yf.records(d) == slice(1, -1)
    got = yf.humidity(d)
    assert got.shape == (4, 1)
    assert np.allclose(got[:, 0], 0.5 * (start[1:-1, 0] + start[2:, 0]))
    assert np.allclose(yf.interval_means(start), got)


def test_fields_take_the_same_records():
    d, start = year_file()
    f = yf.fields(d, ('hus', 'ta', 'lev'))
    assert f['ta'].shape == f['hus'].shape == (4, 1)
    assert np.allclose(f['ta'][:, 0], np.arange(1, 5) + 270.0)     # records 1 to n-2
    assert np.allclose(f['lev'], [0.5])                             # no time axis: whole


def test_a_run_that_writes_the_mean_keeps_every_record():
    d, start = year_file(with_mean=True)
    assert yf.records(d) == slice(None)
    f = yf.fields(d, ('hus', 'ta'))
    assert f['hus'].shape == f['ta'].shape == (6, 1)
    assert np.allclose(f['hus'], start + 5e-4)
