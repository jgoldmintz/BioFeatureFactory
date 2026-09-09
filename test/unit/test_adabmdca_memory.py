"""Parameter parsing and gauge regressions for bounded auxiliary memory."""

import tracemalloc

import numpy as np
import pytest

from biofeaturefactory.mutation_effects.adabmdca_pipeline import (
    apply_zero_sum_gauge,
    load_adabmdca_params,
)


@pytest.fixture
def params_file(tmp_path):
    params_path = tmp_path / "model_params.dat"
    params_path.write_text(
        "J 0 1 - A 1.250000e-01\n"
        "J 0 2 C - -1.750000e+00\n"
        "J 1 2 A C 2.500000e+00\n"
        "h 0 - -1.000000e+00\n"
        "h 0 A 2.500000e-01\n"
        "h 0 C 5.000000e-01\n"
        "h 1 - 0.000000e+00\n"
        "h 1 A 1.500000e+00\n"
        "h 1 C -5.000000e-01\n"
        "h 2 - 7.500000e-01\n"
        "h 2 A -1.250000e+00\n"
        "h 2 C 2.000000e+00\n"
    )
    return params_path


@pytest.mark.parametrize("dtype", [np.float16, np.float32, np.float64, "float64"])
def test_load_params_layout_values_and_dtype(params_file, dtype):
    fields, couplings = load_adabmdca_params(params_file, "-AC", dtype=dtype)
    expected_fields = np.array(
        [[-1, 0.25, 0.5], [0, 1.5, -0.5], [0.75, -1.25, 2]], dtype=dtype,
    )
    expected_couplings = np.zeros((3, 3, 3, 3), dtype=dtype)
    expected_couplings[0, 0, 1, 1] = 0.125
    expected_couplings[1, 1, 0, 0] = 0.125
    expected_couplings[0, 2, 2, 0] = -1.75
    expected_couplings[2, 0, 0, 2] = -1.75
    expected_couplings[1, 1, 2, 2] = 2.5
    expected_couplings[2, 2, 1, 1] = 2.5

    np.testing.assert_array_equal(fields, expected_fields)
    np.testing.assert_array_equal(couplings, expected_couplings)
    assert fields.dtype == np.dtype(dtype)
    assert couplings.dtype == np.dtype(dtype)
    assert couplings.strides == expected_couplings.transpose(0, 2, 1, 3).strides
    assert not couplings.flags.owndata


def test_loader_defaults_to_float32(params_file):
    fields, couplings = load_adabmdca_params(params_file, "-AC")
    assert fields.dtype == np.float32
    assert couplings.dtype == np.float32


def test_load_duplicate_reverse_and_diagonal_records(tmp_path):
    params_path = tmp_path / "duplicates.dat"
    params_path.write_text(
        "h 0 A 1\nh 0 A -2\nh 2 - 0\n"
        "J 0 1 A C 1\nJ 0 1 A C 2\nJ 1 0 C A 3\n"
        "J 1 1 A C 4\nJ 1 1 C A 5\nJ 2 2 C C 6\n"
    )
    fields, couplings = load_adabmdca_params(params_path, "-AC")

    assert fields[0, 1] == -2
    assert couplings[0, 1, 1, 2] == couplings[1, 2, 0, 1] == 3
    assert couplings[1, 1, 1, 2] == couplings[1, 2, 1, 1] == 5
    assert couplings[2, 2, 2, 2] == 6
    assert np.count_nonzero(couplings) == 5


@pytest.mark.parametrize("coupling_record", ["J 2 4 A C 0.25", "J 4 2 A C 0.25"])
def test_coupling_indices_extend_field_dimensions(tmp_path, coupling_record):
    params_path = tmp_path / "extended.dat"
    params_path.write_text(f"h 0 A 1\n{coupling_record}\n")
    fields, couplings = load_adabmdca_params(params_path, "AC")

    assert fields.shape == (5, 2)
    assert couplings.shape == (5, 2, 5, 2)
    assert np.count_nonzero(fields) == 1
    assert np.count_nonzero(couplings) == 2


def test_fields_only_and_ignored_lines(tmp_path):
    params_path = tmp_path / "fields.dat"
    params_path.write_text("\n# metadata\nignored record\nh 2 A 0.5 trailing text\n")
    fields, couplings = load_adabmdca_params(params_path, "AC")

    assert fields.shape == (3, 2)
    assert fields[2, 0] == 0.5
    assert not couplings.any()


@pytest.mark.parametrize("contents", ["", "\n# ignored\n", "J 0 1 A C 1\n"])
def test_missing_fields_error(tmp_path, contents):
    params_path = tmp_path / "missing_fields.dat"
    params_path.write_text(contents)

    with pytest.raises(RuntimeError, match="No 'h' entries found"):
        load_adabmdca_params(params_path, "AC")


@pytest.mark.parametrize(
    ("record", "error_type"),
    [
        ("h invalid A 1", ValueError),
        ("h 0 A invalid", ValueError),
        ("h 0 A", IndexError),
        ("h 0 Z 1", KeyError),
        ("J invalid 1 A C 1", ValueError),
        ("J 0 invalid A C 1", ValueError),
        ("J 0 1 A C invalid", ValueError),
        ("J 0 1 A C", IndexError),
        ("J 0 1 Z C 1", KeyError),
        ("J 0 1 A Z 1", KeyError),
        ("h -3 A 1", IndexError),
        ("J -3 0 A C 1", IndexError),
        ("J 0 -3 A C 1", IndexError),
    ],
)
def test_malformed_record_errors(tmp_path, record, error_type):
    params_path = tmp_path / "malformed.dat"
    params_path.write_text(f"h 1 A 0\n{record}\n")

    with pytest.raises(error_type):
        load_adabmdca_params(params_path, "AC")


def test_valid_negative_indices_retain_numpy_behavior(tmp_path):
    params_path = tmp_path / "negative.dat"
    params_path.write_text("h 1 A 0\nh -1 C 2\nJ -1 0 C A 3\n")
    fields, couplings = load_adabmdca_params(params_path, "AC")

    assert fields[1, 1] == 2
    assert couplings[1, 1, 0, 0] == couplings[0, 0, 1, 1] == 3


def test_loader_auxiliary_memory_does_not_scale_with_record_count(tmp_path):
    peaks = []
    for record_count in (100, 20000):
        params_path = tmp_path / f"repeated_{record_count}.dat"
        with params_path.open("w") as params_stream:
            for record_index in range(record_count):
                params_stream.write(f"h 0 A {record_index}\nJ 0 1 A C {record_index}\n")

        tracemalloc.start()
        try:
            fields, couplings = load_adabmdca_params(params_path, "AC")
            peaks.append(tracemalloc.get_traced_memory()[1])
        finally:
            tracemalloc.stop()

        assert fields[0, 0] == record_count - 1
        assert couplings[0, 0, 1, 1] == record_count - 1

    assert peaks[1] < 256 * 1024
    assert peaks[1] < peaks[0] + 128 * 1024


@pytest.mark.parametrize("dtype", [np.float16, np.float32, np.float64, np.int32])
@pytest.mark.parametrize("layout", ["contiguous", "loader", "strided"])
def test_chunked_gauge_matches_dense_and_preserves_inputs(dtype, layout):
    generator = np.random.default_rng(120)
    fields = generator.normal(size=(67, 3)).astype(dtype)
    couplings = generator.normal(size=(67, 3, 67, 3)).astype(dtype)
    if layout == "loader":
        couplings = couplings.transpose(0, 2, 1, 3).copy().transpose(0, 2, 1, 3)
    elif layout == "strided":
        couplings = couplings[::-1, :, ::-1, :]
    original_fields = fields.copy()
    original_couplings = couplings.copy()
    fields.flags.writeable = False
    couplings.flags.writeable = False
    expected_fields = fields - fields.mean(axis=1, keepdims=True)
    expected_couplings = (
        couplings
        - couplings.mean(axis=1, keepdims=True)
        - couplings.mean(axis=3, keepdims=True)
        + couplings.mean(axis=(1, 3), keepdims=True)
    )

    gauged_fields, gauged_couplings = apply_zero_sum_gauge(fields, couplings)

    np.testing.assert_array_equal(gauged_fields, expected_fields)
    np.testing.assert_array_equal(gauged_couplings, expected_couplings)
    np.testing.assert_array_equal(fields, original_fields)
    np.testing.assert_array_equal(couplings, original_couplings)
    assert gauged_fields.dtype == expected_fields.dtype
    assert gauged_couplings.dtype == expected_couplings.dtype
    assert not np.shares_memory(fields, gauged_fields)
    assert not np.shares_memory(couplings, gauged_couplings)


def test_gauge_accepts_empty_couplings():
    fields = np.empty((0, 3), dtype=np.float32)
    couplings = np.empty((0, 3, 0, 3), dtype=np.float32)
    gauged_fields, gauged_couplings = apply_zero_sum_gauge(fields, couplings)

    assert gauged_fields.shape == fields.shape
    assert gauged_couplings.shape == couplings.shape
    assert gauged_fields.dtype == fields.dtype
    assert gauged_couplings.dtype == couplings.dtype


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_gauge_peak_allocation_is_one_output_plus_bounded_scratch(dtype):
    fields = np.ones((160, 8), dtype=dtype)
    couplings = np.ones((160, 160, 8, 8), dtype=dtype).transpose(0, 2, 1, 3)
    tracemalloc.start()
    try:
        gauged_fields, gauged_couplings = apply_zero_sum_gauge(fields, couplings)
        peak_bytes = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()

    assert not gauged_fields.any()
    assert not gauged_couplings.any()
    output_bytes = gauged_fields.nbytes + gauged_couplings.nbytes
    assert output_bytes <= peak_bytes < output_bytes + 2 * 1024 * 1024
