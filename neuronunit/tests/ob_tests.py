debug = True

from .base import np, pq, sciunit, ncap, VmTest, scores
from matplotlib import pyplot as plt

class SagVoltageTestAngelo2012(VmTest):
    """Tests for sag voltage, as described in Angelo et. al. (2012)"""

    ss_delay = 1000 * pq.ms
    current_duration = 1500 * pq.ms
    testing_voltage = -90 * pq.mV
    sag_window = 100 * pq.ms

    name = "Sag voltage test"

    description = ("A test of the sag voltage of a cell")

    score_type = scores.ZScore

    units = pq.mV

    def __init__(self, *args, **kwargs):
        super(SagVoltageTestAngelo2012, self).__init__(*args, **kwargs)

    def generate_prediction(self, model):

        # First, find the current that will bring the voltage to the testing voltage e.g. -90mV
        # Do this using a voltage clamp
        current = model.clamp_voltage(voltages= [self.testing_voltage,  0, 0]*pq.mV,
                                      durations=[self.current_duration, 0, 0]*pq.ms)

        inhibitory_current = np.median(current[-10:].magnitude)*pq.nA

        # Inject that current
        voltage = model.inject_square_current({"delay":self.ss_delay,
                                     "duration": self.current_duration,
                                     "amplitude": inhibitory_current})

        times = voltage.times

        # Find the minimum v within the ROI after SS and before end of injection
        t_roi = np.where((times > self.ss_delay) & (times <= self.ss_delay + self.current_duration))
        v_roi = voltage.base[t_roi]

        if debug:
            plt.plot(times[t_roi], v_roi)

        v_injection_ss = np.median(v_roi[-10:].magnitude)*pq.mV
        v_min_index = np.argmin(v_roi)
        v_min = v_roi[v_min_index]

        if debug:
            plt.plot([times[t_roi][-1]], [v_injection_ss], 'o', label="Current SS: "+str(v_injection_ss))
            plt.plot([times[t_roi][v_min_index]], [v_min], 'o', label="Min: "+str(v_min))

        if v_min < v_injection_ss:
            sag = v_injection_ss - v_min

        # If there is no sag within the window, use the voltage at specific time point e.g. 100ms after current onset
        else:
            t_roi = np.where((times > self.ss_delay + self.sag_window - 1*pq.ms) & (times < self.ss_delay + self.sag_window + 1*pq.ms))
            v_roi = voltage.base[t_roi]
            v_at_sag_window = np.median(v_roi.magnitude) * pq.mV

            sag = v_injection_ss - v_at_sag_window

            if debug:
                plt.plot([np.median(times[t_roi])], [v_at_sag_window], 'o', label="Alt min: "+str(v_at_sag_window))


        if debug:
            plt.ylim((-92, -87))
            plt.legend()
            plt.show()

        return sag


class SagVoltageTestBurtonUrban2014(SagVoltageTestAngelo2012):
    """Tests for sag voltage, as described in Burton and Urban (2014)"""
    current_duration = 2000 * pq.ms


class SagVoltageTestYu2015(SagVoltageTestBurtonUrban2014):
    """Tests for sag voltage, as described in Yu et. al. (2015)"""
    sag_window = 120 * pq.ms


class SagVoltageTestHu2016(SagVoltageTestAngelo2012):
    """
    Tests for sag voltage, as described in Hu et. al. (2016)

    Paper does not specify what to do if min v is above injection steady
    state v -- here using same sag_window as Angelo et. al. (2012)
    """
    current_duration = 1000 * pq.ms