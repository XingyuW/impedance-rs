def test_import():
    import pymodules.fformat
    import impedance_rs
    assert impedance_rs is not None
    assert pymodules.fformat is not None
