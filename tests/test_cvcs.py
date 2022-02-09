import pytest
from pytest import approx
import demins
from demins import ureg as ureg
from datetime import timedelta, datetime

@pytest.fixture
def make_vessel(density=3785):
    """
    Testing fixture to reproduce the spreadsheet provided by Wayne.
    tv = test vessel
    """
    tv = demins.Vessel()
    tv.anion_vol = 10 * ureg.cubic_foot
    tv.anion_specific_capacity = 1 * ureg.eq / ureg.L
    tv.cation_vol = 20 * ureg.cubic_foot
    tv.cation_specific_capacity = 2.4 * ureg.eq / ureg.L
    tv.flow = 120 * ureg.gallons / ureg.minute
    tv.start_date = datetime.now() - timedelta(days=305)
    boron = demins.Species(conc=1000, name='boron', units='ppm', equivalent_per_mole=1,
            molecular_weight=10.4)
    lithium = demins.Species(conc=3.5, name='lithium', units='ppm', equivalent_per_mole=1,
            molecular_weight=7)
    ammonia = demins.Species(conc=400, name='ammonia', units='ppb', equivalent_per_mole=1,
            molecular_weight=18.04, removal_efficiency=0.9)
    zinc = demins.Species(conc=10, name='zinc', units='ppb', equivalent_per_mole=2,
            molecular_weight=65.38, removal_efficiency=0.98)
    sodium = demins.Species(conc=0.5, name='sodium', units='ppb', equivalent_per_mole=1,
            molecular_weight=22.9898, removal_efficiency=0.98)
    magnesium = demins.Species(conc=0.5, name='magnesium', units='ppb', equivalent_per_mole=2,
            molecular_weight=24.305, removal_efficiency=0.98)
    calcium = demins.Species(conc=0.5, name='calcium', units='ppb', equivalent_per_mole=2,
            molecular_weight=40.08, removal_efficiency=0.98)
    cs137 = demins.Species(conc=4e-5, name='cs137', units='cs137', equivalent_per_mole=1,
            molecular_weight=137, removal_efficiency=0.98)
    co58 = demins.Species(conc=3e-4, name='co58', units='co58', equivalent_per_mole=2,
            molecular_weight=58, removal_efficiency=0.98)
    co60 = demins.Species(conc=4e-5, name='co60', units='co60', equivalent_per_mole=2,
            molecular_weight=60, removal_efficiency=0.98)
    fluoride = demins.Species(conc=0.5, name='fluoride', units='ppb', equivalent_per_mole=1,
            molecular_weight=19)
    chloride = demins.Species(conc=0.5, name='chloride', units='ppb', equivalent_per_mole=1,
            molecular_weight=34.453)
    sulfate = demins.Species(conc=0.5, name='sulfate', units='ppb', equivalent_per_mole=2,
            molecular_weight=96.06, removal_efficiency=0.98)
    i131 = demins.Species(conc=4e-5, name='i131', units='i131', equivalent_per_mole=1,
            molecular_weight=131, removal_efficiency=0.98)
    tv.anions = [boron, fluoride, chloride, sulfate, i131]
    tv.cations = [lithium, ammonia, zinc, sodium, magnesium, calcium, cs137, co58, co60]
    tv.initial_bed_capacity()
    for ion in tv.cations + tv.anions:
        ion.dot_e_calc(tv.flow, density)
        ion.uptake_calc()
    tv.vessel_throughput()
    return tv

def test_cation_uptake_rate(make_vessel):
    make_vessel.cation_uptake_rate_calc()
    assert make_vessel.cation_uptake_rate.magnitude == approx(9.24E-3, rel=1e-3)

def test_anion_uptake_rate(make_vessel):
    make_vessel.anion_uptake_rate_calc()
    assert make_vessel.anion_uptake_rate.magnitude == approx(4.63E-6, rel=1e-3)

def test_cation_capacity_init(make_vessel):
    assert make_vessel.cation_capacity_init.magnitude == approx(1359.2064, rel=1e-3)

def test_anion_capacity_init(make_vessel):
    assert make_vessel.anion_capacity_init.magnitude == approx(283.17, rel=1e-3)

def test_cation_capacity_now(make_vessel):
    assert make_vessel.cation_capacity_calc().magnitude == approx(-2698.62, rel=1e-3)

def test_anion_capacity_now(make_vessel):
    assert make_vessel.anion_capacity_calc().magnitude == approx(281.13, rel=1e-3)

def test_cation_percent_spent(make_vessel):
    make_vessel.cation_percent_exhausted()
    assert make_vessel.cation_pct_exh.magnitude == approx(298.54, rel=1e-3)

def test_anion_capacity_remain(make_vessel):
    make_vessel.anion_percent_exhausted()
    assert make_vessel.anion_pct_exh.magnitude == approx(0.7187, rel=1e-3)

def test_time_to_cation_capacity(make_vessel):
    assert make_vessel.time_to_cation_exhaust().magnitude == approx(-337294.14, rel=1e-3)

def test_time_to_anion(make_vessel):
    assert make_vessel.time_to_anion_exhaust().magnitude == approx(42337925, rel=1e-3)
