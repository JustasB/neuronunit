from neuronunit.tests.olfactory_bulb.utilities import get_APs
from ..base import VmTest
from sciunit import capabilities as scap
from neuronunit import capabilities as ncap

class OlfactoryBulbCellTest(VmTest):
    pass


class OlfactoryBulbCellSpikeTest(OlfactoryBulbCellTest):

    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             ncap.SupportsSettingTemperature,
                             ncap.SupportsSettingStopTime,
                             ncap.ProducesActionPotentials,)

    def get_aps(self, voltage):
        return get_APs(voltage, self.ss_delay, self.threshold_method)