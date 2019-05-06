from neuronunit.tests.olfactory_bulb.publications import *
from neuronunit.tests.olfactory_bulb.tests import *


class RheobaseTestYu2015(Yu2015, RheobaseTest):
    """Tests for rheobase current (1 AP per 2 sec square injection), as described in Yu et. al. (2015)"""
    pass


class RheobaseTestBurtonUrban2014(BurtonUrban2014, RheobaseTest):
    """Tests for sag voltage, as described in Burton and Urban (2014)"""
    pass


class SagVoltageTestAngelo2012(Angelo2012, SagVoltageTest):
    """Tests for sag voltage, as described in Angelo 2012 """


class SagVoltageTestBurtonUrban2014(BurtonUrban2014, SagVoltageTest):
    """Tests for sag voltage, as described in Burton and Urban (2014)  """
    pass


class SagVoltageTestYu2015(Yu2015, SagVoltageTest):
    """Tests for sag voltage, as described in Yu et. al. (2015)"""
    pass


class SagVoltageTestHu2016(Hu2016, SagVoltageTest):
    """
    Tests for sag voltage, as described in Hu et. al. (2016)

    Paper does not specify what to do if min v is above injection steady
    state v -- here using same sag_window as Angelo et. al. (2012) = 100ms
    """
    pass


class AfterHyperpolarizationAmplitudeTestYu2015(Yu2015, AfterHyperpolarizationAmplitudeTest):
    """Tests for After-hyperpolarization (AHP) amplitude, as described in Yu et. al. (2015) """
    pass

class AfterHyperpolarizationTimeTestYu2015(Yu2015, AfterHyperpolarizationTimeTest):
    pass

class AfterHyperpolarizationAmplitudeTestBurtonUrban2014(BurtonUrban2014, AfterHyperpolarizationAmplitudeTest):
    pass

class AfterHyperpolarizationTimeTestBurtonUrban2014(BurtonUrban2014, AfterHyperpolarizationTimeTest):
    pass


class ReboundSpikingTestBurtonUrban2014(BurtonUrban2014,  ReboundSpikingTest):
    pass

class ReboundSpikingTestJohnsonDelaney2010(JohnsonDelaney2010,  ReboundSpikingTest):
    pass

class FISlopeTestBurtonUrban2014(BurtonUrban2014, FISlopeTest):
    pass

class ISICVTestBurtonUrban2014(BurtonUrban2014, ISICVTest):
    pass

class ISICVTestYu2015(Yu2015, ISICVTest):
    pass


class SpikeAccommodationTestBurtonUrban2014(BurtonUrban2014, SpikeAccommodationTest):
    pass

class SpikeAccommodationTestZibman2011(Zibman2011, SpikeAccommodationTest):
    pass


class SpikeAccommodationTimeConstantTestBurtonUrban2014(BurtonUrban2014, SpikeAccommodationTimeConstantTest):
    pass

class SpikeAccommodationTimeConstantTestZibman2011(Zibman2011, SpikeAccommodationTimeConstantTest):
    pass